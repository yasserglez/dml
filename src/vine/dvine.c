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
#include <gsl/gsl_vector_ulong.h>
#include <gsl/gsl_vector_short.h>

#include "src/dml.h"

static void
dvine_select_order(dml_vine_t *vine,
                   const gsl_matrix *data,
                   dml_vine_weight_t weight,
                   dml_measure_t ***measure_matrix,
                   const gsl_rng *rng)
{
    int n = (int) vine->dim;
    double **weight_matrix;
    size_t selected_node;
    gsl_vector_ulong *tour;
    gsl_vector_short *in_tour;
    double selected_cost, current_cost;
    size_t current_pos = 0, selected_pos = 0; // Initialized to avoid GCC warnings.
    double dik, dkj, dij;
    size_t i, j;
    int cut_index;

    // Compute the weights. The weights are minimized.
    weight_matrix = g_malloc_n(n, sizeof(double *));
    for (size_t i = 0; i < n; i++) {
        weight_matrix[i] = g_malloc_n(n, sizeof(double));
        for (size_t j = 0; j < i; j++) {
            switch (weight) {
            case DML_VINE_WEIGHT_TAU:
                weight_matrix[i][j] = 1
                        - fabs(dml_measure_tau_coef(measure_matrix[i][j]));
                break;
            case DML_VINE_WEIGHT_CVM:
                weight_matrix[i][j] = measure_matrix[i][j]->x->size
                        - dml_measure_cvm_stat(measure_matrix[i][j]);
                break;
            default:
                weight_matrix[i][j] = 0;
                break;
            }
            weight_matrix[j][i] = weight_matrix[i][j];
        }
    }

    // Compute an approximate solution for the TSP instance using the
    // Cheapest insertion heuristic. The dummy node has index n.
    tour = gsl_vector_ulong_alloc(n + 1);
    in_tour = gsl_vector_short_alloc(n + 1);
    gsl_vector_short_set_all(in_tour, FALSE);

    for (size_t node_count = 0; node_count < n + 1; node_count++) {
        // Select the node to be inserted.
        if (node_count == 0) {
            // Initial node (randomly selected).
            selected_node = floor((n + 1) * gsl_rng_uniform(rng));
            gsl_vector_ulong_set(tour, node_count, selected_node);
            gsl_vector_short_set(in_tour, selected_node, TRUE);
        } else {
            selected_cost = GSL_DBL_MAX;

            // Rest of the nodes. Choose the nodes with the minimal insertion cost.
            for (size_t k = 0; k < n + 1; k++) {
                if (!gsl_vector_short_get(in_tour, k)) {
                    if (node_count == 1) {
                        i = gsl_vector_ulong_get(tour, 0);
                        current_cost = (i == n || k == n) ? 0: weight_matrix[i][k];
                        current_pos = 0;
                    } else {
                        current_cost = GSL_DBL_MAX;
                        for (size_t pos = 0; pos < node_count - 1; pos++) {
                            i = gsl_vector_ulong_get(tour, pos);
                            j = gsl_vector_ulong_get(tour, pos + 1);
                            dik = (i == n || k == n) ? 0: weight_matrix[i][k];
                            dkj = (k == n || j == n) ? 0: weight_matrix[k][j];
                            dij = (i == n || j == n) ? 0: weight_matrix[i][j];

                            if (dik + dkj + dij < current_cost) {
                                current_cost = dik + dkj + dij;
                                current_pos = pos;
                            }
                        }
                    }

                    // Check the last position.
                    i = gsl_vector_ulong_get(tour, node_count - 1);
                    j = gsl_vector_ulong_get(tour, 0);
                    dik = (i == n || k == n) ? 0: weight_matrix[i][k];
                    dkj = (k == n || j == n) ? 0: weight_matrix[k][j];
                    dij = (i == n || j == n) ? 0: weight_matrix[i][j];

                    if (dik + dkj + dij < current_cost) {
                        current_cost = dik + dkj + dij;
                        current_pos = node_count - 1;
                    }

                    if (current_cost < selected_cost) {
                        selected_node = k;
                        selected_cost = current_cost;
                        selected_pos = current_pos;
                    }

                }
            }

            // Add the selected node to the tour.
            for (size_t pos = selected_pos; pos < node_count; pos++) {
                i = gsl_vector_ulong_get(tour, pos + 1);
                if (pos == selected_pos) {
                    gsl_vector_ulong_set(tour, pos + 1, selected_node);
                } else {
                    gsl_vector_ulong_set(tour, pos + 1, j);
                }
                j = i;
            }
            gsl_vector_short_set(in_tour, selected_node, TRUE);
        }
    }

    // Cut the tour at the dummy node.
    cut_index = -1;
    for (int i = 0; i < n + 1; i++) {
        if (cut_index >= 0) {
            vine->order[i - cut_index - 1] = gsl_vector_ulong_get(tour, i);
        } else if (gsl_vector_ulong_get(tour, i) == n) {
            cut_index = i;
        }
    }
    for (int i = 0; i < cut_index; i++) {
        vine->order[n - cut_index + i] = gsl_vector_ulong_get(tour, i);
    }

    gsl_vector_ulong_free(tour);
    gsl_vector_short_free(in_tour);
    for (size_t i = 0; i < n; i++) {
        g_free(weight_matrix[i]);
    }
    g_free(weight_matrix);
}

// Based on Algorithm 4 of Aas, K. and Czado, C. and Frigessi, A. and Bakken, H.
// Pair-Copula Constructions of Multiple Dependence. Insurance: Mathematics
// and Economics, 2009, Vol. 44, pp. 182-198.

static void
vine_fit_dvine(dml_vine_t *vine,
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
    dml_measure_t ***measure_matrix;
    dml_measure_t *measure;
    gsl_vector ***v;
    dml_copula_t *copula;
    double tree_aic, copula_aic;
    gsl_vector_view *x_view;
    gsl_permutation **ranks;

    m = data->size1;
    n = data->size2;

    // Allocate a matrix with pairwise dependence measures.
    measure_matrix = g_malloc_n(n, sizeof(dml_measure_t **));
    for (size_t i = 0; i < n; i++) {
        measure_matrix[i] = g_malloc0_n(n, sizeof(dml_measure_t *));
    }
    x_view = g_malloc_n(n, sizeof(gsl_vector_view));

    // Calculate pairwise dependence measures of the original variables.
    ranks = g_malloc0_n(n, sizeof(gsl_permutation *));
    x_view[0] = gsl_matrix_column((gsl_matrix *) data, 0);
    for (size_t i = 1; i < n; i++) {
        x_view[i] = gsl_matrix_column((gsl_matrix *) data, i);
        for (size_t j = 0; j < i; j++) {
            x_view[j] = gsl_matrix_column((gsl_matrix *) data, j);
            measure_matrix[i][j] = dml_measure_alloc(&x_view[i].vector, &x_view[j].vector);
            measure_matrix[j][i] = measure_matrix[i][j];

            // Pre-calculate the ranks of the variables.
            if (ranks[i] == NULL) {
                ranks[i] = dml_measure_x_rank(measure_matrix[i][j]);
            } else {
                measure_matrix[i][j]->x_rank = gsl_permutation_alloc(
                        measure_matrix[i][j]->x->size);
                gsl_permutation_memcpy(measure_matrix[i][j]->x_rank, ranks[i]);
            }
            if (ranks[j] == NULL) {
                ranks[j] = dml_measure_y_rank(measure_matrix[i][j]);
            } else {
                measure_matrix[i][j]->y_rank = gsl_permutation_alloc(
                        measure_matrix[i][j]->y->size);
                gsl_permutation_memcpy(measure_matrix[i][j]->y_rank, ranks[j]);
            }
        }
    }
    g_free(ranks);

    // Select the order of the variables in the vine. The order of the variables
    // determines the structure of the first tree (and therefore the rest of
    // the vine).
    if (n > 2) {
        dvine_select_order(vine, data, weight, measure_matrix, rng);
    } else {
        vine->order[0] = 0;
        vine->order[1] = 1;
    }

    // This allocates more memory than required. In the future, the code should
    // be modified to use 0-based indexes for the v vector.
    v = g_malloc_n(n, sizeof(gsl_vector **));
    for (size_t i = 0; i <= n - 1; i++) {
        v[i] = g_malloc0_n(GSL_MAX(n, 2*n-4) + 1, sizeof(gsl_vector *));
    }

    // Estimation of the first tree of the vine.
    for (size_t i = 1; i <= n; i++) {
        v[0][i] = gsl_vector_alloc(m);
        gsl_matrix_get_col(v[0][i], data, vine->order[i-1]);
    }

    // Selection of the copulas in the first tree.
    tree_aic = 0;
    for (size_t i = 1; i <= n - 1; i++) {
        measure = measure_matrix[vine->order[i-1]][vine->order[i+1-1]];
        copula = dml_copula_select(v[0][i], v[0][i+1], measure, indeptest,
                                   indeptest_level, types, types_size,
                                   select, rng);
        vine->copulas[0][i-1] = copula;
        // Get information for the truncation of the vine.
        if (trunc == DML_VINE_TRUNC_AIC) {
            dml_copula_aic(copula, v[0][i], v[0][i+1], &copula_aic);
            tree_aic += copula_aic;
        }
    }

    // Free the matrix with pairwise dependence measures.
    for (size_t i = 1; i < n; i++) {
        for (size_t j = 0; j < i; j++) {
            dml_measure_free(measure_matrix[i][j]);
        }
        g_free(measure_matrix[i]);
    }
    g_free(measure_matrix);
    g_free(x_view);

    // Check if the vine should be truncated.
    if (trunc == DML_VINE_TRUNC_AIC && tree_aic >= 0) {
        // Undo the construction of the first tree, free memory and
        // return a vine without trees.
        for (size_t i = 0; i < n - 1; i++) {
            dml_copula_free(vine->copulas[0][i]);
            vine->copulas[0][i] = NULL;
        }
        for (size_t i = 1; i <= n; i++) {
            gsl_vector_free(v[0][i]);
            g_free(v[i-1]);
        }
        g_free(v);
        vine->trees = 0;
        return;
    }

    // Is it a 2-dimensional D-vine?
    if (n == 2) {
        for (size_t i = 1; i <= n; i++) {
            gsl_vector_free(v[0][i]);
            g_free(v[i-1]);
        }
        g_free(v);
        vine->trees = n - 1;
        return;
    }

    // Compute observations for the second tree.
    v[1][1] = gsl_vector_alloc(m);
    dml_copula_h(vine->copulas[0][0], v[0][1], v[0][2], v[1][1]);
    for (size_t k = 1; k <= n - 3; k++) {
        if (v[1][2*k] == NULL)
            v[1][2*k] = gsl_vector_alloc(m);
        dml_copula_h(vine->copulas[0][k], v[0][k+2], v[0][k+1], v[1][2*k]);
        if (v[1][2*k+1] == NULL)
            v[1][2*k+1] = gsl_vector_alloc(m);
        dml_copula_h(vine->copulas[0][k], v[0][k+1], v[0][k+2], v[1][2*k+1]);
    }

    if (v[1][2*n-4] == NULL)
        v[1][2*n-4] = gsl_vector_alloc(m);
    dml_copula_h(vine->copulas[0][n-2], v[0][n], v[0][n-1], v[1][2*n-4]);

    // Estimation of the rest of the trees.
    for (size_t j = 2; j <= n - 1; j++) { // Tree index.
        // Select the copulas of the j-th tree.
        tree_aic = 0;
        for (size_t i = 1; i <= n - j; i++) { // Edges of the tree.
            copula = dml_copula_select(v[j-1][2*i-1], v[j-1][2*i], NULL,
                                       indeptest, indeptest_level, types,
                                       types_size, select, rng);
            vine->copulas[j-1][i-1] = copula;
            // Get information for the truncation of the vine.
            if (trunc == DML_VINE_TRUNC_AIC) {
                dml_copula_aic(copula, v[j-1][2*i-1], v[j-1][2*i], &copula_aic);
                tree_aic += copula_aic;
            }
        }
        // Check if the vine should be truncated.
        if (trunc == DML_VINE_TRUNC_AIC && tree_aic >= 0) {
            for (size_t i = 1; i <= n - j; i++) {
                dml_copula_free(vine->copulas[j-1][i-1]);
                vine->copulas[j-1][i-1] = NULL;
            }
            vine->trees = j - 1;
            break;
        }

        // Last tree of the vine.
        if (j == n - 1) {
            vine->trees = n - 1;
            break;
        }

        // Compute observations for the next tree.
        if (v[j][1] == NULL) v[j][1] = gsl_vector_alloc(m);
        dml_copula_h(vine->copulas[j-1][0], v[j-1][1], v[j-1][2], v[j][1]);
        if (n > 4) {
            for (size_t i = 1; i <= n - j - 2; i++) {
                if (v[j][2*i] == NULL)
                    v[j][2*i] = gsl_vector_alloc(m);
                dml_copula_h(vine->copulas[j-1][i], v[j-1][2*i+2],
                             v[j-1][2*i+1], v[j][2*i]);

                if (v[j][2*i+1] == NULL)
                    v[j][2*i+1] = gsl_vector_alloc(m);
                dml_copula_h(vine->copulas[j-1][i], v[j-1][2*i+1],
                             v[j-1][2*i+2], v[j][2*i+1]);
            }
        }
        if (v[j][2*n-2*j-2] == NULL)
            v[j][2*n-2*j-2] = gsl_vector_alloc(m);
        dml_copula_h(vine->copulas[j-1][n-j-1], v[j-1][2*n-2*j],
                     v[j-1][2*n-2*j-1], v[j][2*n-2*j-2]);
    }

    // Freeing memory.
    for (size_t i = 0; i <= n - 1; i++) {
        for (size_t j = 0; j < 2*n-4+1; j++) {
            if (v[i][j] != NULL) {
                gsl_vector_free(v[i][j]);
            }
        }
        g_free(v[i]);
    }
    g_free(v);
}

// Based on Algorithm 2 of Aas, K. and Czado, C. and Frigessi, A. and Bakken, H.
// Pair-Copula Constructions of Multiple Dependence. Insurance: Mathematics
// and Economics, 2009, Vol. 44, pp. 182-198.

static void
vine_ran_dvine(const dml_vine_t *vine,
               const gsl_rng *rng,
               gsl_matrix *data)
{
    size_t n;
    gsl_matrix *v;
    gsl_vector *w;
    gsl_vector *x, *y, *r;

    n = vine->dim;
    v = gsl_matrix_alloc(n+1, GSL_MAX(n, 2*n-4) + 1);
    w = gsl_vector_alloc(n);
    x = gsl_vector_alloc(1);
    y = gsl_vector_alloc(1);
    r = gsl_vector_alloc(1);

    for (size_t s = 0; s < data->size1; s++) { // Loop over samples.
        for (size_t i = 0; i < n; i++) {
            gsl_vector_set(w, i, gsl_rng_uniform(rng));
        }

        gsl_matrix_set(v, 1, 1, gsl_vector_get(w, 0));
        gsl_matrix_set(data, s, vine->order[0], gsl_vector_get(w, 0));

        gsl_vector_set(x, 0, gsl_vector_get(w, 1));
        gsl_vector_set(y, 0, gsl_matrix_get(v, 1, 1));
        dml_copula_hinv(vine->copulas[0][0], x, y, r);
        gsl_matrix_set(v, 2, 1, gsl_vector_get(r, 0));
        gsl_matrix_set(data, s, vine->order[1], gsl_vector_get(r, 0));

        gsl_vector_set(x, 0, gsl_matrix_get(v, 1, 1));
        gsl_vector_set(y, 0, gsl_matrix_get(v, 2, 1));
        dml_copula_h(vine->copulas[0][0], x, y, r);
        gsl_matrix_set(v, 2, 2, gsl_vector_get(w, 1));

        for (size_t i = 3; i <= n; i++) { // Loop over the rest of the variables.
            gsl_matrix_set(v, i, 1, gsl_vector_get(w, i-1));

            if (vine->trees >= 2) {
                for (size_t k = GSL_MIN(vine->trees, i - 1); k >= 2; k--) {
                    gsl_vector_set(x, 0, gsl_matrix_get(v, i, 1));
                    gsl_vector_set(y, 0, gsl_matrix_get(v, i-1, 2*k-2));
                    dml_copula_hinv(vine->copulas[k-1][i-k-1], x, y, r);
                    gsl_matrix_set(v, i, 1, gsl_vector_get(r, 0));
                }
            }

            gsl_vector_set(x, 0, gsl_matrix_get(v, i, 1));
            gsl_vector_set(y, 0, gsl_matrix_get(v, i-1, 1));
            dml_copula_hinv(vine->copulas[0][i-2], x, y, r);
            gsl_matrix_set(v, i, 1, gsl_vector_get(r, 0));
            gsl_matrix_set(data, s, vine->order[i-1], gsl_vector_get(r, 0));

            if (i == n) break;

            if (vine->trees >= 2) {
                gsl_vector_set(x, 0, gsl_matrix_get(v, i-1, 1));
                gsl_vector_set(y, 0, gsl_matrix_get(v, i, 1));
                dml_copula_h(vine->copulas[0][i-2], x, y, r);
                gsl_matrix_set(v, i, 2, gsl_vector_get(r, 0));
            }

            if (vine->trees >= 3) {
                gsl_vector_set(x, 0, gsl_matrix_get(v, i, 1));
                gsl_vector_set(y, 0, gsl_matrix_get(v, i-1, 1));
                dml_copula_h(vine->copulas[0][i-2], x, y, r);
                gsl_matrix_set(v, i, 3, gsl_vector_get(r, 0));
            }

            if (vine->trees >= 3 && i > 3) {
                for (size_t j = 2; j <= GSL_MIN(vine->trees - 1, i - 2); j++) {
                    gsl_vector_set(x, 0, gsl_matrix_get(v, i-1, 2*j-2));
                    gsl_vector_set(y, 0, gsl_matrix_get(v, i, 2*j-1));
                    dml_copula_h(vine->copulas[j-1][i-j-1], x, y, r);
                    gsl_matrix_set(v, i, 2*j, gsl_vector_get(r, 0));

                    gsl_vector_set(x, 0, gsl_matrix_get(v, i, 2*j-1));
                    gsl_vector_set(y, 0, gsl_matrix_get(v, i-1, 2*j-2));
                    dml_copula_h(vine->copulas[j-1][i-j-1], x, y, r);
                    gsl_matrix_set(v, i, 2*j+1, gsl_vector_get(r, 0));
                }
            }

            if (vine->trees >= i) {
                gsl_vector_set(x, 0, gsl_matrix_get(v, i-1, 2*i-4));
                gsl_vector_set(y, 0, gsl_matrix_get(v, i, 2*i-3));
                dml_copula_h(vine->copulas[i-2][0], x, y, r);
                gsl_matrix_set(v, i, 2*i-2, gsl_vector_get(r, 0));
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
vine_free_dvine(dml_vine_t *vine)
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
dml_vine_alloc_dvine(const size_t dim)
{
    dml_vine_t *vine;

    vine = g_malloc(sizeof(dml_vine_t));
    vine->type = DML_VINE_DVINE;
    vine->dim = dim;
    vine->trees = 0;
    vine->order = g_malloc_n(dim, sizeof(size_t));
    vine->matrix = NULL; // Not used. D-vine represented by the order of the variables.
    // Upper triangular matrix with the parameters of the copulas.
    vine->copulas = g_malloc_n(dim - 1, sizeof(dml_copula_t **));
    for (size_t i = 0; i < dim - 1; i++) {
        vine->copulas[i] = g_malloc0_n(dim - 1 - i, sizeof(dml_copula_t *));
    }
    vine->fit = vine_fit_dvine;
    vine->ran = vine_ran_dvine;
    vine->free = vine_free_dvine;

    return vine;
}
