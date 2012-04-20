// DML - Dependence Modeling Library
// Copyright (C) 2011-2012 Yasser González-Fernández <ygonzalezfernandez@gmail.com>

#include "config.h"

#include <glib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <igraph/igraph.h>

#include "src/dml.h"

// Macros to save pointers in the numeric attributes of igraph. Portability?

#define VAP(graph, n, v) ((void *) (long) (igraph_cattribute_VAN((graph), (n), (v))))
#define EAP(graph, n, e) ((void *) (long) (igraph_cattribute_EAN((graph), (n), (e))))
#define SETVAP(graph, n, vid, value) (igraph_cattribute_VAN_set((graph),(n),(vid),((double) (long) (value))))
#define SETEAP(graph, n, eid, value) (igraph_cattribute_EAN_set((graph),(n),(eid),((double) (long) (value))))

static double
cvine_calculate_weight(dml_vine_weight_t weight,
                       dml_measure_t *measure)
{
    double value;

    // The weight is maximized.

    switch (weight) {
    case DML_VINE_WEIGHT_TAU:
        value = fabs(dml_measure_tau_coef(measure));
        break;
    case DML_VINE_WEIGHT_CVM:
        value = dml_measure_cvm_stat(measure);
        break;
    default:
        value = 0;
        break;
    }

    return value;
}

static inline void
cvine_tree_cleanup(igraph_t *tree)
{
    igraph_integer_t a;
    gsl_vector *xa;

    for (a = 0; a < igraph_vcount(tree); a++) {
        xa = VAP(tree, "data", a);
        gsl_vector_free(xa);
    }
    DELVAS(tree);
}

static void
vine_fit_cvine(dml_vine_t *vine,
               const gsl_matrix *data,
               const dml_vine_weight_t weight,
               const dml_vine_trunc_t trunc,
               const dml_copula_indeptest_t indeptest,
               const double indeptest_level,
               const dml_copula_type_t *types,
               const size_t types_size,
               const dml_copula_select_t select,
               const double gof_level,
               const gsl_rng *rng)
{
    size_t m, n;
    igraph_t **trees;
    gsl_vector *x;
    size_t k;
    igraph_integer_t e; // Edge id.
    igraph_integer_t a, b; // Vertex id.
    igraph_integer_t root_vertex = -1; // Vertex id of the root of the tree.
    size_t root_index; // Variable index corresponding to the root of the tree.
    gsl_vector *xa, *xb; // Samples of 'a' and 'b', respectively.
    dml_measure_t *measure, ***measure_matrix;
    dml_copula_t *copula;
    double tree_weight, max_tree_weight;
    double tree_aic, copula_aic;
    gsl_vector_short *selected_roots;
    gsl_permutation **ranks;

    igraph_i_set_attribute_table(&igraph_cattribute_table);

    m = data->size1;
    n = data->size2;
    trees = g_malloc0_n(n - 1, sizeof(igraph_t *));
    selected_roots = gsl_vector_short_calloc(n);

    for (k = 0; k < n - 1; k++) { // Tree index.
        trees[k] = g_malloc(sizeof(igraph_t));
        igraph_empty(trees[k], n - k, IGRAPH_UNDIRECTED);

        // Assign original observations or computed pseudo-observations to the nodes.
        if (k == 0) {
            for (size_t i = 0; i < n; i++) { // Variable and node index.
                x = gsl_vector_alloc(m);
                gsl_matrix_get_col(x, data, i);
                SETVAP(trees[k], "data", i, x); // Observations.
                // Index of the original variable represented by this node.
                SETVAN(trees[k], "index", i, i);
            }
        } else {
            for (e = 0; e < igraph_ecount(trees[k - 1]); e++) {
                igraph_edge(trees[k - 1], e, &a, &b);
                copula = EAP(trees[k - 1], "copula", e);
                xa = VAP(trees[k - 1], "data", a);
                xb = VAP(trees[k - 1], "data", b);
                x = gsl_vector_alloc(m);
                if (a == root_vertex) {
                    dml_copula_h(copula, xb, xa, x);
                } else {
                    dml_copula_h(copula, xa, xb, x);
                }
                // The edge id of a tree are the vertex id of the next tree.
                SETVAP(trees[k], "data", e, x); // Pseudo-observations.
                SETVAN(trees[k], "index", e, EAN(trees[k - 1], "index", e));
            }
        }

        // Allocate the matrix for pairwise dependence measures.
        measure_matrix = g_malloc_n(n - k, sizeof(dml_measure_t **));
        for (size_t i = 0; i < n - k; i++) {
            measure_matrix[i] = g_malloc0_n(n - k, sizeof(dml_measure_t *));
        }

        // Select the root node of the tree.
        ranks = g_malloc0_n(n - k, sizeof(gsl_permutation *));
        max_tree_weight = GSL_NAN;
        for (a = 0; a < n - k; a++) {
            xa = VAP(trees[k], "data", a);
            tree_weight = 0;
            for (b = 0; b < n - k; b++) {
                if (a != b) {
                    xb = VAP(trees[k], "data", b);
                    if (measure_matrix[(size_t) a][(size_t) b] == NULL) {
                        measure_matrix[(size_t) a][(size_t) b] = dml_measure_alloc(xa, xb);
                        measure = measure_matrix[(size_t) a][(size_t) b];
                        measure_matrix[(size_t) b][(size_t) a] = measure;

                        // Pre-calculate the ranks of the variables.
                        if (ranks[(size_t) a] == NULL) {
                            ranks[(size_t) a] = dml_measure_x_rank(measure);
                        } else {
                            measure->x_rank = gsl_permutation_alloc(measure->x->size);
                            gsl_permutation_memcpy(measure->x_rank, ranks[(size_t) a]);
                        }
                        if (ranks[(size_t) b] == NULL) {
                            ranks[(size_t) b] = dml_measure_y_rank(measure);
                        } else {
                            measure->y_rank = gsl_permutation_alloc(measure->y->size);
                            gsl_permutation_memcpy(measure->y_rank, ranks[(size_t) b]);
                        }
                    }
                    tree_weight += cvine_calculate_weight(weight,
                            measure_matrix[(size_t) a][(size_t) b]);
                }
            }
            if (gsl_isnan(max_tree_weight) || tree_weight > max_tree_weight) {
                root_vertex = a;
                max_tree_weight = tree_weight;
            }
        }
        g_free(ranks);

        // Update the order of the variables and add edges from the root
        // vertex to the other vertexes.
        root_index = VAN(trees[k], "index", root_vertex);
        gsl_vector_short_set(selected_roots, root_index, 1);
        vine->order[k] = root_index;
        for (b = 0; b < igraph_vcount(trees[k]); b++) {
            if (root_vertex != b) {
                igraph_add_edge(trees[k], root_vertex, b);
                igraph_get_eid(trees[k], &e, root_vertex, b, 0);
                SETEAN(trees[k], "index", e, VAN(trees[k], "index", b));
            }
        }

        // Select the copulas in the tree.
        tree_aic = 0;
        for (e = 0; e < igraph_ecount(trees[k]); e++) {
            igraph_edge(trees[k], e, &a, &b);
            xa = VAP(trees[k], "data", a);
            xb = VAP(trees[k], "data", b);
            if (a == root_vertex) {
                copula = dml_copula_select(
                        xb, xa, measure_matrix[(size_t) b][(size_t) a],
                        indeptest, indeptest_level, types, types_size, select,
                        gof_level, rng);
            } else {
                copula = dml_copula_select(
                        xa, xb, measure_matrix[(size_t) a][(size_t) b],
                        indeptest, indeptest_level, types, types_size, select,
                        gof_level, rng);
            }
            SETEAP(trees[k], "copula", e, copula);

            // Get information for the truncation of the vine.
            if (trunc == DML_VINE_TRUNC_AIC) {
                dml_copula_aic(copula, xa, xb, &copula_aic);
                tree_aic += copula_aic;
            }
        }

        // Free the matrix for pairwise dependence measures and the memory
        // used by the attributes of the previous tree.
        if (k > 0) cvine_tree_cleanup(trees[k - 1]);
        for (size_t i = 0; i < n - k; i++) {
            for (size_t j = 0; j < i; j++) {
                dml_measure_free(measure_matrix[i][j]);
            }
            g_free(measure_matrix[i]);
        }
        g_free(measure_matrix);

        // Check if the vine should be truncated.
        if (trunc == DML_VINE_TRUNC_AIC && tree_aic >= 0) {
            // Undo the construction of the last tree.
            gsl_vector_short_set(selected_roots, root_index, 0);
            cvine_tree_cleanup(trees[k]);
            for (e = 0; e < igraph_ecount(trees[k]); e++) {
                copula = EAP(trees[k], "copula", e);
                dml_copula_free(copula);
            }
            igraph_destroy(trees[k]);
            g_free(trees[k]);
            trees[k] = NULL;
            break;
        }
    }

    // If the vine was completely estimated, cleanup the last tree.
    if (k == n - 1) {
        cvine_tree_cleanup(trees[k - 1]);
    }

    // Set the number of trees and complete the vector with the order of the variables.
    vine->trees = k;
    for (size_t i = 0; i < n; i++) {
        if (!gsl_vector_short_get(selected_roots, i)) {
            vine->order[k++] = i;
        }
    }

    // Initialize the matrix with the copulas and free the memory of the trees.
    for (k = 0; k < vine->trees; k++) {
        for (size_t i = 0; i < n - k - 1; i++) {
            for (e = 0; e < igraph_ecount(trees[k]); e++) {
                if (EAN(trees[k], "index", e) == vine->order[k + i + 1]) {
                    copula = EAP(trees[k], "copula", e);
                    vine->copulas[k][i] = copula;
                    break;
                }
            }
        }
        igraph_destroy(trees[k]);
        g_free(trees[k]);
    }
    g_free(trees);
    gsl_vector_short_free(selected_roots);
}

// Based on Algorithm 1 of Aas, K. and Czado, C. and Frigessi, A. and Bakken, H.
// Pair-Copula Constructions of Multiple Dependence. Insurance: Mathematics
// and Economics, 2009, Vol. 44, pp. 182-198.

static void
vine_ran_cvine(const dml_vine_t *vine,
               const gsl_rng *rng,
               gsl_matrix *data)
{
    size_t n;
    gsl_matrix *v;
    gsl_vector *w;
    gsl_vector *x, *y, *r;

    n = vine->dim;
    v = gsl_matrix_alloc(n, n);
    w = gsl_vector_alloc(n);
    x = gsl_vector_alloc(1);
    y = gsl_vector_alloc(1);
    r = gsl_vector_alloc(1);

    for (size_t s = 0; s < data->size1; s++) { // Loop over samples.
        for (size_t i = 0; i < n; i++) {
            gsl_vector_set(w, i, gsl_rng_uniform(rng));
        }
        gsl_matrix_set(data, s, vine->order[0], gsl_vector_get(w, 0));
        gsl_matrix_set(v, 0, 0, gsl_vector_get(w, 0));
        for (size_t i = 1; i < n; i++) { // Loop over the variables.
            gsl_matrix_set(v, i, 0, gsl_vector_get(w, i));
            for (size_t k = GSL_MIN(i - 1, vine->trees - 1); /* See break call. */; k--) {
                gsl_vector_set(x, 0, gsl_matrix_get(v, i, 0));
                gsl_vector_set(y, 0, gsl_matrix_get(v, k, k));
                dml_copula_hinv(vine->copulas[k][i - k - 1], x, y, r);
                gsl_matrix_set(v, i, 0, gsl_vector_get(r, 0));
                if (k == 0) break; // Avoid problems decrementing k if k is 0.
            }
            gsl_matrix_set(data, s, vine->order[i], gsl_matrix_get(v, i, 0));
            if (i == n - 1) {
                break;
            }
            for (size_t j = 0; j <= GSL_MIN(i - 1, vine->trees - 1); j++) {
                gsl_vector_set(x, 0, gsl_matrix_get(v, i, j));
                gsl_vector_set(y, 0, gsl_matrix_get(v, j, j));
                dml_copula_h(vine->copulas[j][i - j - 1], x, y, r);
                gsl_matrix_set(v, i, j + 1, gsl_vector_get(r, 0));
            }
        }
    }

    gsl_matrix_free(v);
    gsl_vector_free(w);
    gsl_vector_free(x);
    gsl_vector_free(y);
    gsl_vector_free(r);
}

static void
vine_free_cvine(dml_vine_t *vine)
{
    g_free(vine->order);
    if (vine->copulas != NULL) {
        for (size_t i = 0; i < vine->dim - 1; i++) {
            for (size_t j = 0; j < vine->dim - 1 - i; j++) {
                if (vine->copulas[i][j] != NULL) {
                    dml_copula_free(vine->copulas[i][j]);
                }
            }
            g_free(vine->copulas[i]);
        }
        g_free(vine->copulas);
    }
}

dml_vine_t *
dml_vine_alloc_cvine(const size_t dim)
{
    dml_vine_t *vine;

    vine = g_malloc(sizeof(dml_vine_t));
    vine->type = DML_VINE_CVINE;
    vine->dim = dim;
    vine->trees = 0;
    vine->order = g_malloc_n(dim, sizeof(size_t));
    vine->matrix = NULL; // Not used. C-vine represented by the order of the variables.
    // Upper triangular matrix with the parameters of the copulas.
    vine->copulas = g_malloc_n(dim - 1, sizeof(dml_copula_t **));
    for (size_t i = 0; i < dim - 1; i++) {
        vine->copulas[i] = g_malloc0_n(dim - 1 - i, sizeof(dml_copula_t *));
    }
    vine->fit = vine_fit_cvine;
    vine->ran = vine_ran_cvine;
    vine->free = vine_free_cvine;

    return vine;
}
