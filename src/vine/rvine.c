// DML - Dependence Modeling Library
// Copyright (C) 2011-2012 Yasser González-Fernández <ygonzalezfernandez@gmail.com>

#include "config.h"

#include <glib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_vector_short.h>
#include <igraph/igraph.h>

#include "src/dml.h"

// Macros to save pointers in the numeric attributes of igraph. Portability?
#define VAP(graph, n, v) ((void *) (long) (igraph_cattribute_VAN((graph), (n), (v))))
#define EAP(graph, n, e) ((void *) (long) (igraph_cattribute_EAN((graph), (n), (e))))
#define SETVAP(graph, n, vid, value) (igraph_cattribute_VAN_set((graph),(n),(vid),((double) (long) (value))))
#define SETEAP(graph, n, eid, value) (igraph_cattribute_EAN_set((graph),(n),(eid),((double) (long) (value))))

static void
rvine_set_weight(igraph_t *graph,
                 dml_vine_weight_t weight,
                 igraph_integer_t e,
                 const gsl_vector *x,
                 const gsl_vector *y)
{
    double value;
    dml_measure_t *measure = NULL;

    switch (weight) {
    case DML_VINE_WEIGHT_TAU:
        measure = dml_measure_alloc(x, y);
        value = 1 - fabs(dml_measure_tau_coef(measure));
        break;
    default:
        value = 0;
        break;
    }

    SETEAN(graph, "weight", e, value);
    SETEAP(graph, "tau", e, measure);
}

static inline void
rvine_tree_cleanup(igraph_t *tree)
{
    igraph_integer_t e, a;
    gsl_vector *xa;
    gsl_vector_short *Ue;
    dml_measure_t *measure;

    for (a = 0; a < igraph_vcount(tree); a++) {
        xa = VAP(tree, "data", a);
        gsl_vector_free(xa);
    }
    DELVAS(tree);

    for (e = 0; e < igraph_ecount(tree); e++) {
        Ue = EAP(tree, "Ue", e);
        gsl_vector_short_free(Ue);
        measure = EAP(tree, "tau", e);
        if (measure != NULL) dml_measure_free(measure);
    }
    DELEA(tree, "Ue");
    DELEA(tree, "weight");
    DELEA(tree, "tau");
}

static void
fit_rvine_trees(igraph_t **trees,
                const gsl_matrix *data,
                dml_vine_weight_t weight,
                dml_vine_truncation_t truncation,
                dml_copula_indeptest_t indeptest,
                double indeptest_level,
                const dml_copula_type_t *types,
                size_t types_size,
                dml_copula_selection_t selection)
{
    size_t m, n;
    igraph_t *graph;
    igraph_vector_t *graph_weight;
    dml_copula_t *copula;
    gsl_vector *x;
    igraph_integer_t e; // Edge id.
    igraph_integer_t a, b, aa, ab, ba, bb; // Vertex id.
    gsl_vector *xa, *xb; // Samples of 'a' and 'b', respectively.
    gsl_vector_short *Ue, *Ua, *Ub;
    size_t k;
    dml_measure_t *measure;
    double tree_aic, copula_aic;

    igraph_i_set_attribute_table(&igraph_cattribute_table);

    m = data->size1;
    n = data->size2;
    graph = g_malloc(sizeof(igraph_t));
    graph_weight = g_malloc(sizeof(igraph_vector_t));

    for (k = 0; k < n - 1; k++) { // Tree index.
        if (k == 0) {
            igraph_full(graph, n, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

            // Assign the observations to the nodes.
            for (size_t i = 0; i < n; i++) { // Variable and node index.
                x = gsl_vector_alloc(m);
                gsl_matrix_get_col(x, data, i);
                SETVAP(graph, "data", i, x);
            }

            for (e = 0; e < igraph_ecount(graph); e++) {
                igraph_edge(graph, e, &a, &b);

                // Calculate the weight of the edge.
                xa = VAP(graph, "data", a);
                xb = VAP(graph, "data", b);
                rvine_set_weight(graph, weight, e, xa, xb);

                // Variables "connected" by this edge.
                Ue = gsl_vector_short_calloc(n);
                gsl_vector_short_set(Ue, a, 1);
                gsl_vector_short_set(Ue, b, 1);
                SETEAP(graph, "Ue", e, Ue);

                // Conditioned set.
                SETEAN(graph, "Cea", e, a + 1);
                SETEAN(graph, "Ceb", e, b + 1);
            }
        } else {
            igraph_empty(graph, n - k, IGRAPH_UNDIRECTED);

            // Compute the pseudo-observations.
            for (e = 0; e < igraph_ecount(trees[k - 1]); e++) {
                igraph_edge(trees[k - 1], e, &a, &b);
                copula = EAP(trees[k - 1], "copula", e);
                xa = VAP(trees[k - 1], "data", a);
                xb = VAP(trees[k - 1], "data", b);
                x = gsl_vector_alloc(m);
                dml_copula_h(copula, xa, xb, x);
                // The edge id of a tree are the vertex id of the next tree.
                SETVAP(graph, "data", e, x);
            }

            // Adding all the "possible" edges.
            for (a = 0; a < igraph_vcount(graph) - 1; a++) {
                for (b = a + 1; b < igraph_vcount(graph); b++) {
                    // Checking the proximity condition.
                    igraph_edge(trees[k - 1], a, &aa, &ab);
                    igraph_edge(trees[k - 1], b, &ba, &bb);
                    if (aa == ba || aa == bb || ab == ba || ab == bb) {
                        igraph_add_edge(graph, a, b);
                        igraph_get_eid(graph, &e, a, b, 0);

                        // Calculate the weight of the edge.
                        xa = VAP(graph, "data", a);
                        xb = VAP(graph, "data", b);
                        rvine_set_weight(graph, weight, e, xa, xb);

                        // Variables "connected" by this edge and conditioned set.
                        Ua = EAP(trees[k - 1], "Ue", a);
                        Ub = EAP(trees[k - 1], "Ue", b);
                        Ue = gsl_vector_short_calloc(n);
                        for (size_t i = 0; i < n; i++) {
                            gsl_vector_short_set(Ue, i,
                                    gsl_vector_short_get(Ua, i)
                                            | gsl_vector_short_get(Ub, i));
                            if (gsl_vector_short_get(Ua, i)
                                    && !gsl_vector_short_get(Ub, i)) {
                                SETEAN(graph, "Cea", e, i + 1);
                            }
                            if (gsl_vector_short_get(Ub, i)
                                    && !gsl_vector_short_get(Ua, i)) {
                                SETEAN(graph, "Ceb", e, i + 1);
                            }
                        }
                        SETEAP(graph, "Ue", e, Ue);
                    }
                }
            }
        }

        // Compute the minimum weight spanning tree.
        trees[k] = g_malloc(sizeof(igraph_t));
        igraph_vector_init(graph_weight, igraph_ecount(graph));
        EANV(graph, "weight", graph_weight);
        igraph_minimum_spanning_tree_prim(graph, trees[k], graph_weight);
        igraph_vector_destroy(graph_weight);

        tree_aic = 0;
        for (e = 0; e < igraph_ecount(trees[k]); e++) {
            igraph_edge(trees[k], e, &a, &b);
            xa = VAP(trees[k], "data", a);
            xb = VAP(trees[k], "data", b);
            measure = EAP(trees[k], "tau", e);

            // Assign bivariate copulas to the edges.
            copula = dml_copula_select(xa, xb, measure, indeptest,
                                          indeptest_level, types, types_size,
                                          selection);
            SETEAP(trees[k], "copula", e, copula);

            // Get information for the truncation of the vine.
            if (truncation == DML_VINE_TRUNCATION_AIC) {
                dml_copula_aic(copula, xa, xb, &copula_aic);
                tree_aic += copula_aic;
            }
        }

        if (k > 0) rvine_tree_cleanup(trees[k - 1]);
        igraph_destroy(graph);

        // Check if the vine should be truncated.
        if (truncation == DML_VINE_TRUNCATION_AIC && tree_aic >= 0) {
            // Free the memory used for the last tree.
            rvine_tree_cleanup(trees[k]);
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

    // Cleanup the last tree if the vine was completely estimated.
    if (k == n - 1) {
        rvine_tree_cleanup(trees[k - 1]);
    }

    g_free(graph_weight);
    g_free(graph);
}

/* Compute an R-vine matrix from the trees. Based on Algorithm 3.2 of
 * Dissman, J. F. (2010). Statistical Inference for Regular Vines and
 * Application. Diploma thesis. University of Technology, Munich.
 */
static void
rvine_trees_to_vine(dml_vine_t *vine, igraph_t **trees)
{
    size_t n;
    size_t *order_inv;
    gsl_vector_short *B;
    igraph_integer_t Cea, Ceb;
    size_t x = 0, x_hat = 0, x_hat_hat = 0; // Initialized to avoid GCC warnings.
    dml_copula_t *copula = NULL; // Initialized to avoid GCC warnings.
    igraph_integer_t e; // Edge id.

    n = vine->dimension;
    order_inv = g_malloc0_n(n, sizeof(size_t));
    B = gsl_vector_short_calloc(n);

    vine->trees = n - 1;
    while (trees[vine->trees - 1] == NULL)
        vine->trees--;

    // for loop in line 2.
    for (size_t i = 0; i < n - 1; i++) {
        if (trees[n - i - 2] == NULL) {
            // Truncated. Get an unassigned variable from the last tree.
            for (e = 0; e < igraph_ecount(trees[vine->trees - 1]); e++) {
                x = EAN(trees[vine->trees - 1], "Cea", e);
                x_hat = EAN(trees[vine->trees - 1], "Ceb", e);
                if (!gsl_vector_short_get(B, x - 1)
                        && !gsl_vector_short_get(B, x_hat - 1)) {
                    // Mark the truncated entries with 0 and NULL.
                    x_hat = 0;
                    copula = NULL;
                    break;
                }
            }
        } else {
            for (e = 0; e < igraph_ecount(trees[n - i - 2]); e++) {
                x = EAN(trees[n - i - 2], "Cea", e);
                x_hat = EAN(trees[n - i - 2], "Ceb", e);
                if (!gsl_vector_short_get(B, x - 1)
                        && !gsl_vector_short_get(B, x_hat - 1)) {
                    copula = EAP(trees[n - i - 2], "copula", e);
                    break;
                }
            }
        }

        // Line 4.
        gsl_vector_short_set(B, x - 1, 1);
        vine->order[n - i - 1] = x - 1;
        order_inv[x - 1] = n - i;
        vine->matrix[i][i] = x;
        vine->matrix[i + 1][i] = x_hat;
        vine->copulas[i + 1][i] = copula;

        // for loop in line 5.
        for (size_t k = i + 2; k < n; k++) {
            if (trees[n - k - 1] != NULL) {
                for (e = 0; e < igraph_ecount(trees[n - k - 1]); e++) {
                    Cea = EAN(trees[n - k - 1], "Cea", e);
                    Ceb = EAN(trees[n - k - 1], "Ceb", e);
                    if (x == Cea) {
                        x_hat_hat = Ceb;
                        if (!gsl_vector_short_get(B, x_hat_hat - 1)) {
                            /* The pseudocode of the algorithm does not
                             * included this check. Invalid matrices
                             * were generated when xhathat is set to an
                             * index already assigned to a diagonal entry.
                             */
                            copula = EAP(trees[n - k - 1], "copula", e);
                            break;
                        }
                    } else if (x == Ceb) {
                        x_hat_hat = Cea;
                        if (!gsl_vector_short_get(B, x_hat_hat - 1)) {
                            // Ibdem to the previous comment.
                            copula = EAP(trees[n - k - 1], "copula", e);
                            break;
                        }
                    }
                }
                vine->matrix[k][i] = x_hat_hat;
                vine->copulas[k][i] = copula;
            }
        }
    }

    for (size_t i = 0; i < n; i++) {
        if (!gsl_vector_short_get(B, i)) {
            vine->matrix[n - 1][n - 1] = i + 1;
            vine->order[0] = i;
            order_inv[i] = 1;
            break;
        }
    }
    // Reorder the variables. The simulation algorithm assumes that the
    // diagonal entries of the R-vine matrix are ordered from n to 1.
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j <= i; j++) {
            if (vine->matrix[i][j] > 0) {
                vine->matrix[i][j] = order_inv[vine->matrix[i][j] - 1];
            }
        }
    }

    g_free(order_inv);
    gsl_vector_short_free(B);
}

/* Sequential R-vine estimation method. Based on Algorithm 4.1 of
 * Dissman, J. F. (2010). Statistical Inference for Regular Vines and
 * Application. Diploma thesis. University of Technology, Munich.
 */
static void
vine_fit_rvine(dml_vine_t *vine,
               const gsl_matrix *data,
               dml_vine_weight_t weight,
               dml_vine_truncation_t truncation,
               dml_copula_indeptest_t indeptest,
               double indeptest_level,
               const dml_copula_type_t *types,
               size_t types_size,
               dml_copula_selection_t selection)
{
    igraph_t **trees;

    trees = g_malloc0_n(data->size2 - 1, sizeof(igraph_t *));
    fit_rvine_trees(trees, data, weight, truncation, indeptest, indeptest_level,
                    types, types_size, selection);
    rvine_trees_to_vine(vine, trees);

    for (size_t i = 0; i < data->size2 - 1; i++) {
        if (trees[i] != NULL) {
            igraph_destroy(trees[i]);
            g_free(trees[i]);
        } else {
            break;
        }
    }
    g_free(trees);
}

/* Simulation of observations of an R-vine specification. Based on
 * Algorithm 5.2 of Dissman, J. F. (2010). Statistical Inference
 * for Regular Vines and Application. Diploma thesis. University
 * of Technology, Munich.
 */
static void
vine_ran_rvine(const dml_vine_t *vine, const gsl_rng *rng, gsl_matrix *data)
{
    size_t n, m;
    gsl_vector ***vdirect, ***vindirect;
    gsl_vector *z1 = NULL, *z2 = NULL, *hinv = NULL; // Initialized to avoid GCC warnings.
    size_t **M;

    n = vine->dimension;
    m = data->size1;
    vdirect = g_malloc_n(n, sizeof(gsl_vector **));
    vindirect = g_malloc_n(n, sizeof(gsl_vector **));
    for (size_t i = 0; i < n; i++) {
        vdirect[i] = g_malloc0_n(n, sizeof(gsl_vector *));
        vindirect[i] = g_malloc0_n(n, sizeof(gsl_vector *));
    }
    M = g_malloc_n(n, sizeof(size_t *));
    for (size_t i = 0; i < n; i++) {
        M[i] = g_malloc0_n(i + 1, sizeof(size_t));
    }

    // Line 4.
    for (size_t k = 0; k < n; k++) {
        vdirect[n - 1][k] = gsl_vector_alloc(m);
        for (size_t i = 0; i < m; i++) {
            gsl_vector_set(vdirect[n - 1][k], i, gsl_rng_uniform(rng));
        }
    }
    // Line 5.
    for (size_t k = 0; k < n; k++) {
        M[k][k] = vine->matrix[k][k];
        M[n - 1][k] = vine->matrix[n - 1][k];
        for (size_t i = n - 2; i > k; i--) {
            if (vine->matrix[i][k] > M[i + 1][k]) {
                M[i][k] = vine->matrix[i][k];
            } else {
                M[i][k] = M[i + 1][k];
            }
        }
    }
    // Line 6.
    gsl_matrix_set_col(data, vine->order[0], vdirect[n - 1][n - 1]);

    // for loop in line 7.
    for (size_t k = n - 2; /* See break call. */; k--) {
        // for loop in line 8.
        for (size_t i = k + 1; i < n; i++) {
            // Line 14.
            if (vine->matrix[i][k] != 0
                    && dml_copula_type(vine->copulas[i][k]) != DML_COPULA_INDEP) {
                if (M[i][k] == vine->matrix[i][k]) {
                    z2 = vdirect[i][n - M[i][k]];
                } else {
                    z2 = vindirect[i][n - M[i][k]];
                }
                hinv = gsl_vector_alloc(m);
                dml_copula_hinv(vine->copulas[i][k], vdirect[n - 1][k], z2,
                                hinv);
                gsl_vector_free(vdirect[n - 1][k]);
                vdirect[n - 1][k] = hinv;
            }
        }
        // Line 16.
        gsl_matrix_set_col(data, vine->order[n - k - 1], vdirect[n - 1][k]);

        if (k == 0) break; // Avoid problems decrementing the unsigned k if k is 0.

        // for loop in line 17.
        for (size_t i = n - 1; i > k; i--) {
            // Line 18.
            z1 = vdirect[i][k];
            // Line 19.
            if (vdirect[i - 1][k] == NULL) {
                vdirect[i - 1][k] = gsl_vector_alloc(m);
            }
            if (vine->matrix[i][k] == 0
                    || dml_copula_type(vine->copulas[i][k]) == DML_COPULA_INDEP) {
                // Vine truncated or Independence copula.
                gsl_vector_memcpy(vdirect[i - 1][k], z1);
            } else {
                dml_copula_h(vine->copulas[i][k], z1, z2, vdirect[i - 1][k]);
            }
            if (vindirect[i - 1][k] == NULL) {
                vindirect[i - 1][k] = gsl_vector_alloc(m);
            }
            gsl_vector_memcpy(vindirect[i - 1][k], vdirect[i - 1][k]);
        }
    }

    // Freeing memory.
    for (size_t i = 0; i < n; i++) {
        g_free(M[i]);
    }
    g_free(M);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (vdirect[i][j] != NULL) {
                gsl_vector_free(vdirect[i][j]);
            }
            if (vindirect[i][j] != NULL) {
                gsl_vector_free(vindirect[i][j]);
            }
        }
        g_free(vdirect[i]);
        g_free(vindirect[i]);
    }
    g_free(vdirect);
    g_free(vindirect);
}

static void
vine_free_rvine(dml_vine_t *vine)
{
    g_free(vine->order);
    if (vine->matrix != NULL) {
        for (size_t i = 0; i < vine->dimension; i++) {
            g_free(vine->matrix[i]);
        }
        g_free(vine->matrix);
    }
    if (vine->copulas != NULL) {
        for (size_t i = 1; i < vine->dimension; i++) {
            for (size_t j = 0; j < i; j++) {
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
dml_vine_alloc_rvine(const size_t dimension)
{
    dml_vine_t *vine;

    vine = g_malloc(sizeof(dml_vine_t));
    vine->type = DML_VINE_RVINE;
    vine->dimension = dimension;
    vine->trees = 0;
    vine->order = g_malloc_n(dimension, sizeof(size_t));
    // R-vine matrix.
    vine->matrix = g_malloc_n(dimension, sizeof(size_t *));
    for (size_t i = 0; i < dimension; i++) {
        vine->matrix[i] = g_malloc0_n(i + 1, sizeof(size_t));
    }
    // Lower triangular matrix with the copulas.
    vine->copulas = g_malloc_n(dimension, sizeof(dml_copula_t **));
    vine->copulas[0] = NULL;
    for (size_t i = 1; i < dimension; i++) {
        vine->copulas[i] = g_malloc0_n(i, sizeof(dml_copula_t *));
    }
    vine->fit = vine_fit_rvine;
    vine->ran = vine_ran_rvine;
    vine->free = vine_free_rvine;

    return vine;
}
