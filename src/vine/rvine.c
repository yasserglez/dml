/*
 * DML - Dependence Modeling Library
 * Copyright (C) 2011-2013 Yasser Gonzalez <contact@yassergonzalez.com>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"

#include <glib.h>
#include <gsl/gsl_math.h>
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
                 const gsl_vector *y,
                 const gsl_permutation *x_rank,
                 const gsl_permutation *y_rank)
{
    double value;
    dml_measure_t *measure;

    // The weight is minimized.

    measure = dml_measure_alloc(x, y);
    measure->x_rank = gsl_permutation_alloc(x_rank->size);
    gsl_permutation_memcpy(measure->x_rank, x_rank);
    measure->y_rank = gsl_permutation_alloc(y_rank->size);
    gsl_permutation_memcpy(measure->y_rank, y_rank);

    switch (weight) {
    case DML_VINE_WEIGHT_TAU:
        value = 1 - fabs(dml_measure_tau_coef(measure));
        break;
    case DML_VINE_WEIGHT_CVM:
        value = measure->x->size - dml_measure_cvm_stat(measure);
        break;
    default:
        value = 0;
        break;
    }

    SETEAN(graph, "weight", e, value);
    SETEAP(graph, "measure", e, measure);
}

static inline void
rvine_tree_cleanup(igraph_t *tree)
{
    igraph_integer_t e, a;

    for (a = 0; a < igraph_vcount(tree); a++) {
        if (VAP(tree, "h", a) != VAP(tree, "hrev", a)) {
            gsl_vector_free(VAP(tree, "h", a));
            gsl_vector_free(VAP(tree, "hrev", a));
            gsl_permutation_free(VAP(tree, "hrank", a));
            gsl_permutation_free(VAP(tree, "hrevrank", a));
        } else if (VAP(tree, "h", a) != NULL) {
            gsl_vector_free(VAP(tree, "h", a));
            gsl_permutation_free(VAP(tree, "hrank", a));
        } else {
            gsl_vector_free(VAP(tree, "hrev", a));
            gsl_permutation_free(VAP(tree, "hrevrank", a));
        }
    }
    DELVAS(tree);

    for (e = 0; e < igraph_ecount(tree); e++) {
        gsl_vector_short_free(EAP(tree, "Ue", e));
        dml_measure_free(EAP(tree, "measure", e));
    }
    DELEA(tree, "Ue");
    DELEA(tree, "weight");
    DELEA(tree, "measure");
}

// Select an R-vine (R-vine trees) to model the given sample. Based on
// Algorithm 4.1 of Dissman, J. F. (2010). Statistical Inference for Regular
// Vines and Application. Diploma thesis. University of Technology, Munich.

static void
fit_rvine_trees(igraph_t **trees,
                const gsl_matrix *data,
                const dml_vine_weight_t weight,
                const dml_vine_trunc_t trunc,
                const dml_copula_indeptest_t indeptest,
                const double indeptest_level,
                const dml_copula_type_t *types,
                const size_t types_size,
                const dml_copula_select_t select,
                const gsl_rng *rng)
{
    size_t m, n;
    igraph_t *graph;
    igraph_vector_t *graph_weight;
    dml_copula_t *copula;
    gsl_vector *x;
    igraph_integer_t e; // Edge id.
    igraph_integer_t a, aa, ab, b, ba, bb; // Vertex id.
    gsl_vector *u = NULL, *v = NULL;
    igraph_integer_t Cea, Ceb;
    gsl_vector_short *Ue, *Ua, *Ub;
    size_t k;
    dml_measure_t *measure;
    double tree_aic, copula_aic;
    gsl_permutation *perm, *rank, *u_rank = NULL, *v_rank = NULL;

    igraph_i_set_attribute_table(&igraph_cattribute_table);

    m = data->size1;
    n = data->size2;
    graph = g_malloc(sizeof(igraph_t));
    graph_weight = g_malloc(sizeof(igraph_vector_t));
    perm = gsl_permutation_alloc(m);

    for (k = 0; k < n - 1; k++) { // Tree index.
        if (k == 0) {
            igraph_full(graph, n, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

            // Assign the observations to the nodes.
            for (size_t i = 0; i < n; i++) { // Variable and node index.
                x = gsl_vector_alloc(m);
                gsl_matrix_get_col(x, data, i);

                // Results of the h-function of the copula assigned to the
                // edge that corresponds to this vertex in the previous tree.
                // h for the h-function with its arguments in order and
                // hrev for the h-function with its arguments reversed. In the
                // first tree both are equal to the observations of the
                // corresponding variable, in the rest of the trees they differ.
                SETVAP(graph, "h", i, x);
                SETVAP(graph, "hrev", i, x);
                gsl_sort_vector_index(perm, x);
                rank = gsl_permutation_alloc(m);
                gsl_permutation_inverse(rank, perm);
                // Ranks of the h and hrev vectors.
                SETVAP(graph, "hrank", i, rank);
                SETVAP(graph, "hrevrank", i, rank);
            }

            for (e = 0; e < igraph_ecount(graph); e++) {
                igraph_edge(graph, e, &a, &b);

                // Variables "connected" by this edge.
                Ue = gsl_vector_short_calloc(n);
                gsl_vector_short_set(Ue, a, 1);
                gsl_vector_short_set(Ue, b, 1);
                SETEAP(graph, "Ue", e, Ue);

                // Conditioned set.
                SETEAN(graph, "Cea", e, a + 1);
                SETEAN(graph, "Ceb", e, b + 1);
                Cea = EAN(graph, "Cea", e);
                Ceb = EAN(graph, "Ceb", e);

                // Calculate the weight of the edge.
                u = VAP(graph, "h", a);
                v = VAP(graph, "h", b);
                u_rank = VAP(graph, "hrank", a);
                v_rank = VAP(graph, "hrank", b);
                // The conditioned set is ordered to make the order of the
                // arguments in the bivariate copulas unique as suggested in
                // Czado, C. (2010) Pair-Copula Constructions of Multivariate
                // Copulas. In Jaworski, P. and Durante, F. and Hardle, W. K.
                // and Rychlik, T. (eds.) Copula Theory and Its Applications,
                // Springer-Verlag, 93-109.
                if (Cea < Ceb) {
                    rvine_set_weight(graph, weight, e, u, v, u_rank, v_rank);
                } else {
                    rvine_set_weight(graph, weight, e, v, u, v_rank, u_rank);
                }
            }
        } else {
            igraph_empty(graph, n - k, IGRAPH_UNDIRECTED);

            // Adding all "possible" edges.
            for (a = 0; a < igraph_vcount(graph) - 1; a++) {
                for (b = a + 1; b < igraph_vcount(graph); b++) {
                    igraph_edge(trees[k - 1], a, &aa, &ab);
                    igraph_edge(trees[k - 1], b, &ba, &bb);

                    // Checking the proximity condition.
                    if (aa == ba || aa == bb || ab == ba || ab == bb) {
                        igraph_add_edge(graph, a, b);
                        igraph_get_eid(graph, &e, a, b, IGRAPH_UNDIRECTED);

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

            // Compute pseudo-observations and edge weights.
            for (a = 0; a < igraph_vcount(graph); a++) {
                // See the comment in the code for the first tree.
                SETVAP(graph, "h", a, NULL);
                SETVAP(graph, "hrev", a, NULL);
                SETVAP(graph, "hrank", a, NULL);
                SETVAP(graph, "hrevrank", a, NULL);
            }
            for (e = 0; e < igraph_ecount(graph); e++) {
                igraph_edge(graph, e, &a, &b);
                Cea = EAN(graph, "Cea", e);
                Ceb = EAN(graph, "Ceb", e);

                // Assign u and u_rank.
                if ((Cea == EAN(trees[k - 1], "Cea", a)
                        && (EAN(trees[k - 1], "Cea", a)
                                < EAN(trees[k - 1], "Ceb", a)))
                        || (Cea != EAN(trees[k - 1], "Cea", a)
                                && (EAN(trees[k - 1], "Cea", a)
                                        > EAN(trees[k - 1], "Ceb", a)))) {
                    u = VAP(graph, "h", a);
                    if (u == NULL) {
                        copula = EAP(trees[k - 1], "copula", a);
                        measure = EAP(trees[k - 1], "measure", a);
                        u = gsl_vector_alloc(m);
                        dml_copula_h(copula, measure->x, measure->y, u);
                        SETVAP(graph, "h", a, u);
                        gsl_sort_vector_index(perm, u);
                        rank = gsl_permutation_alloc(m);
                        gsl_permutation_inverse(rank, perm);
                        SETVAP(graph, "hrank", a, rank);
                    }
                    u_rank = VAP(graph, "hrank", a);
                }
                if ((Cea == EAN(trees[k - 1], "Cea", a)
                        && (EAN(trees[k - 1], "Cea", a)
                                > EAN(trees[k - 1], "Ceb", a)))
                        || (Cea != EAN(trees[k - 1], "Cea", a)
                                && (EAN(trees[k - 1], "Cea", a)
                                        < EAN(trees[k - 1], "Ceb", a)))) {
                    u = VAP(graph, "hrev", a);
                    if (u == NULL) {
                        copula = EAP(trees[k - 1], "copula", a);
                        measure = EAP(trees[k - 1], "measure", a);
                        u = gsl_vector_alloc(m);
                        dml_copula_h(copula, measure->y, measure->x, u);
                        SETVAP(graph, "hrev", a, u);
                        gsl_sort_vector_index(perm, u);
                        rank = gsl_permutation_alloc(m);
                        gsl_permutation_inverse(rank, perm);
                        SETVAP(graph, "hrevrank", a, rank);
                    }
                    u_rank = VAP(graph, "hrevrank", a);
                }

                // Assign v and v_rank.
                if ((Ceb == EAN(trees[k - 1], "Cea", b)
                        && (EAN(trees[k - 1], "Cea", b)
                                < EAN(trees[k - 1], "Ceb", b)))
                        || (Ceb != EAN(trees[k - 1], "Cea", b)
                                && (EAN(trees[k - 1], "Cea", b)
                                        > EAN(trees[k - 1], "Ceb", b)))) {
                    v = VAP(graph, "h", b);
                    if (v == NULL) {
                        copula = EAP(trees[k - 1], "copula", b);
                        measure = EAP(trees[k - 1], "measure", b);
                        v = gsl_vector_alloc(m);
                        dml_copula_h(copula, measure->x, measure->y, v);
                        SETVAP(graph, "h", b, v);
                        gsl_sort_vector_index(perm, v);
                        rank = gsl_permutation_alloc(m);
                        gsl_permutation_inverse(rank, perm);
                        SETVAP(graph, "hrank", b, rank);
                    }
                    v_rank = VAP(graph, "hrank", b);

                }
                if ((Ceb == EAN(trees[k - 1], "Cea", b)
                        && (EAN(trees[k - 1], "Cea", b)
                                > EAN(trees[k - 1], "Ceb", b)))
                        || (Ceb != EAN(trees[k - 1], "Cea", b)
                                && (EAN(trees[k - 1], "Cea", b)
                                        < EAN(trees[k - 1], "Ceb", b)))) {
                    v = VAP(graph, "hrev", b);
                    if (v == NULL) {
                        copula = EAP(trees[k - 1], "copula", b);
                        measure = EAP(trees[k - 1], "measure", b);
                        v = gsl_vector_alloc(m);
                        dml_copula_h(copula, measure->y, measure->x, v);
                        SETVAP(graph, "hrev", b, v);
                        gsl_sort_vector_index(perm, v);
                        rank = gsl_permutation_alloc(m);
                        gsl_permutation_inverse(rank, perm);
                        SETVAP(graph, "hrevrank", b, rank);
                    }
                    v_rank = VAP(graph, "hrevrank", b);
                }

                // Set the weight of the edge. The arguments are ordered here.
                // The order determines the x and y fields of measure.
                if (Cea < Ceb) {
                    rvine_set_weight(graph, weight, e, u, v, u_rank, v_rank);
                } else {
                    rvine_set_weight(graph, weight, e, v, u, v_rank, u_rank);
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
            Cea = EAN(trees[k], "Cea", e);
            Ceb = EAN(trees[k], "Ceb", e);
            measure = EAP(trees[k], "measure", e);

            // Assign a bivariate copula to the edge.
            if (Cea < Ceb) {
                copula = dml_copula_select(measure->x, measure->y, measure,
                                           indeptest, indeptest_level, types,
                                           types_size, select, rng);
                // Get information for the truncation of the vine.
                if (trunc == DML_VINE_TRUNC_AIC) {
                    dml_copula_aic(copula, measure->x, measure->y, &copula_aic);
                    tree_aic += copula_aic;
                }
            } else {
                copula = dml_copula_select(measure->y, measure->x, measure,
                                           indeptest, indeptest_level, types,
                                           types_size, select, rng);
                // Get information for the truncation of the vine.
                if (trunc == DML_VINE_TRUNC_AIC) {
                    dml_copula_aic(copula, measure->y, measure->x, &copula_aic);
                    tree_aic += copula_aic;
                }
            }
            SETEAP(trees[k], "copula", e, copula);
        }

        igraph_destroy(graph);

        // Check if the vine should be truncated.
        if (trunc == DML_VINE_TRUNC_AIC && tree_aic >= 0) {
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

        if (k > 0) rvine_tree_cleanup(trees[k - 1]);
    }

    // Cleanup the last tree if the vine was completely estimated.
    // If the vine was truncated, the last tree will be freed in
    // the function vine_fit_rvine, because the rvine_trees_to_vine
    // function needs some attributes of its edges.
    if (k == n - 1) {
        rvine_tree_cleanup(trees[n - 2]);
    }

    g_free(graph_weight);
    g_free(graph);
    gsl_permutation_free(perm);
}

// Compute an R-vine matrix from the trees. Based on Algorithm 3.2 of
// Dissman, J. F. (2010). Statistical Inference for Regular Vines and
// Application. Diploma thesis. University of Technology, Munich.

static void
rvine_trees_to_vine(dml_vine_t *vine, igraph_t **trees)
{
    size_t n = vine->dim;
    size_t *order_inv;
    gsl_vector_short *B;
    igraph_integer_t Cea, Ceb;
    size_t x = 0, x_hat = 0, x_hat_hat = 0; // Initialized to avoid GCC warnings.
    dml_copula_t *copula = NULL; // Initialized to avoid GCC warnings.
    igraph_integer_t e; // Edge id.
    igraph_integer_t a, b, aa, ab, ba, bb; // Vertex id.
    igraph_t **last_trees = NULL;
    igraph_t *graph = NULL;
    gsl_vector_short *Ue, *Ua, *Ub;

    // Set the number of trees of the vines.
    vine->trees = n - 1;
    while (trees[vine->trees - 1] == NULL) vine->trees--;

    // Nothing to do for vines without trees.
    if (vine->trees == 0) return;

    // Selecting a structure for the trees that were truncated.
    // Is this really necessary? Think a better solution.
    if (vine->trees != n - 1) {
        igraph_i_set_attribute_table(&igraph_cattribute_table);
        last_trees = g_malloc_n(n - 1 - vine->trees, sizeof(igraph_t *));
        graph = g_malloc(sizeof(igraph_t));

        for (size_t k = vine->trees; k < n - 1; k++) { // Tree index.
            igraph_empty(graph, n - k, IGRAPH_UNDIRECTED);

            // Adding all "possible" edges.
            for (a = 0; a < igraph_vcount(graph) - 1; a++) {
                for (b = a + 1; b < igraph_vcount(graph); b++) {
                    // Checking the proximity condition.
                    igraph_edge(k <= vine->trees ? trees[k - 1] : last_trees[k - 1 - vine->trees], a, &aa, &ab);
                    igraph_edge(k <= vine->trees ? trees[k - 1] : last_trees[k - 1 - vine->trees], b, &ba, &bb);
                    if (aa == ba || aa == bb || ab == ba || ab == bb) {
                        igraph_add_edge(graph, a, b);
                        igraph_get_eid(graph, &e, a, b, IGRAPH_UNDIRECTED);

                        // Variables "connected" by this edge and conditioned set.
                        Ua = EAP(k <= vine->trees ? trees[k - 1] : last_trees[k - 1 - vine->trees], "Ue", a);
                        Ub = EAP(k <= vine->trees ? trees[k - 1] : last_trees[k - 1 - vine->trees], "Ue", b);
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

            // Compute the minimum weight spanning tree.
            last_trees[k - vine->trees] = g_malloc(sizeof(igraph_t));
            igraph_minimum_spanning_tree_unweighted(graph, last_trees[k - vine->trees]);

            igraph_destroy(graph);
        }
    }

    order_inv = g_malloc0_n(n, sizeof(size_t));
    B = gsl_vector_short_calloc(n);

    // for loop in line 2.
    for (size_t i = 0; i < n - 1; i++) {
        if (trees[n - i - 2] == NULL) {
            for (e = 0; e < igraph_ecount(last_trees[n - i - 2 - vine->trees]); e++) {
                x = EAN(last_trees[n - i - 2 - vine->trees], "Cea", e);
                x_hat = EAN(last_trees[n - i - 2 - vine->trees], "Ceb", e);
                if (!gsl_vector_short_get(B, x - 1)
                        && !gsl_vector_short_get(B, x_hat - 1)) {
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
                            // The pseudocode of the algorithm does not included
                            // this check. Invalid matrices were generated when
                            // x_hat_hat is set to an index already assigned
                            // to a diagonal entry.
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

    if (vine->trees != n - 1) {
        for (size_t i = 0; i < n - 1 - vine->trees; i++) {
            for (e = 0; e < igraph_ecount(last_trees[i]); e++) {
                Ue = EAP(last_trees[i], "Ue", e);
                gsl_vector_short_free(Ue);
            }
            DELEA(last_trees[i], "Ue");
            igraph_destroy(last_trees[i]);
            g_free(last_trees[i]);
        }
        g_free(last_trees);
        g_free(graph);
    }
    g_free(order_inv);
    gsl_vector_short_free(B);
}

// Sequential R-vine estimation method. Based on Algorithm 4.1 of
// Dissman, J. F. (2010). Statistical Inference for Regular Vines and
// Application. Diploma thesis. University of Technology, Munich.

static void
vine_fit_rvine(dml_vine_t *vine,
               const gsl_matrix *data,
               const dml_vine_weight_t weight,
               const dml_vine_trunc_t trunc,
               const dml_copula_indeptest_t indeptest,
               const double indeptest_level,
               const dml_copula_type_t *types,
               const size_t types_size,
               const dml_copula_select_t select,
               const gsl_rng *rng)
{
    igraph_t **trees;

    trees = g_malloc0_n(data->size2 - 1, sizeof(igraph_t *));
    fit_rvine_trees(trees, data, weight, trunc, indeptest, indeptest_level,
                    types, types_size, select, rng);
    rvine_trees_to_vine(vine, trees);

    // If the vine was truncated, free the memory of the last tree.
    if (vine->trees > 0 && vine->trees != data->size2 - 1) {
        rvine_tree_cleanup(trees[vine->trees - 1]);
    }
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

// Simulation of observations of an R-vine specification. Based on
// Algorithm 5.2 of Dissman, J. F. (2010). Statistical Inference
// for Regular Vines and Application. Diploma thesis. University
// of Technology, Munich.

static void
vine_ran_rvine(const dml_vine_t *vine, const gsl_rng *rng, gsl_matrix *data)
{
    size_t n, m;
    gsl_vector ***vdirect, ***vindirect;
    gsl_vector *z1 = NULL, *z2 = NULL, *hinv = NULL; // Initialized to avoid GCC warnings.
    size_t **M;

    n = vine->dim;
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
                dml_copula_hinv(vine->copulas[i][k], vdirect[n - 1][k], z2, hinv);
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
                // Vine truncated or independence copula.
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
        for (size_t i = 0; i < vine->dim; i++) {
            g_free(vine->matrix[i]);
        }
        g_free(vine->matrix);
    }
    if (vine->copulas != NULL) {
        for (size_t i = 1; i < vine->dim; i++) {
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
dml_vine_alloc_rvine(const size_t dim)
{
    dml_vine_t *vine;

    vine = g_malloc(sizeof(dml_vine_t));
    vine->type = DML_VINE_RVINE;
    vine->dim = dim;
    vine->trees = 0;
    vine->order = g_malloc_n(dim, sizeof(size_t));
    // R-vine matrix.
    vine->matrix = g_malloc_n(dim, sizeof(size_t *));
    for (size_t i = 0; i < dim; i++) {
        vine->matrix[i] = g_malloc0_n(i + 1, sizeof(size_t));
    }
    // Lower triangular matrix with the copulas.
    vine->copulas = g_malloc_n(dim, sizeof(dml_copula_t **));
    vine->copulas[0] = NULL;
    for (size_t i = 1; i < dim; i++) {
        vine->copulas[i] = g_malloc0_n(i, sizeof(dml_copula_t *));
    }
    vine->fit = vine_fit_rvine;
    vine->ran = vine_ran_rvine;
    vine->free = vine_free_rvine;

    return vine;
}
