/* DML - Dependence Modeling Library
 * Copyright (C) 2011 Yasser González-Fernández <ygonzalezfernandez@gmail.com>
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
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "config.h"
#include "src/dml.h"
#include "tests/common.h"

void
test_cvine_alloc()
{
    dml_vine_t *vine;

    vine = dml_vine_alloc(DML_VINE_CVINE, 2);
    g_assert(dml_vine_type(vine) == DML_VINE_CVINE);
    dml_vine_free(vine);
}

void
test_cvine_ran_2d()
{
    gsl_rng *rng;
    dml_vine_t *vine;
    gsl_matrix *data;
    size_t m = 100;
    gsl_vector_view x, y;
    double corr = 0.9, vine_corr;

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, g_test_rand_int());

    vine = dml_vine_alloc(DML_VINE_CVINE, 2);
    vine->trees = 1;
    vine->order[0] = 0;
    vine->order[1] = 1;
    vine->copulas[0][0] = dml_copula_alloc_normal(corr);
    data = gsl_matrix_alloc(m, 2);
    dml_vine_ran(vine, rng, data);
    x = gsl_matrix_column(data, 0);
    y = gsl_matrix_column(data, 1);
    vine_corr = gsl_stats_correlation(x.vector.data, x.vector.stride,
                                      y.vector.data, y.vector.stride, m);
    g_assert(fabs(corr - vine_corr) < 0.1);

    dml_vine_free(vine);
    gsl_matrix_free(data);
}

void
test_cvine_fit_2d()
{
    gsl_rng *rng;
    dml_copula_t *copula;
    gsl_matrix *data;
    size_t m = 100;
    gsl_vector_view x, y;
    dml_vine_t *vine;
    dml_copula_type_t types[] = { DML_COPULA_NORMAL };
    size_t types_size = 1;
    double corr = 0.9, vine_corr;

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, g_test_rand_int());

    data = gsl_matrix_alloc(m, 2);
    x = gsl_matrix_column(data, 0);
    y = gsl_matrix_column(data, 1);

    copula = dml_copula_alloc_normal(corr);
    dml_copula_ran(copula, rng, &x.vector, &y.vector);
    vine = dml_vine_alloc(DML_VINE_CVINE, 2);
    dml_vine_fit(vine, data, DML_VINE_WEIGHT_TAU, DML_VINE_TRUNCATION_NONE,
                 DML_COPULA_INDEPTEST_NONE, 1.0, &types[0], types_size,
                 DML_COPULA_SELECTION_AIC);
    vine_corr = ((double *) vine->copulas[0][0]->data)[0];
    g_assert(fabs(corr - vine_corr) < 0.1);

    gsl_matrix_free(data);
    dml_copula_free(copula);
    dml_vine_free(vine);
}

void
test_cvine_ran_fit_20d_normal_trunc()
{
    size_t n = 20, m = 2500;
    gsl_rng *rng;
    dml_copula_type_t types[] = { DML_COPULA_NORMAL };
    size_t types_size = 1;
    dml_vine_t *vine, *fitted;
    gsl_matrix *vine_data, *fitted_data;
    gsl_vector_view xi, xj;
    double corr;

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, g_test_rand_int());

    vine = dml_vine_alloc(DML_VINE_CVINE, n);
    vine->trees = 1;
    for (size_t i = 0; i < n; i++) {
        vine->order[i] = i;
        if (i < n - 1) {
            if (i % 2 == 0) {
                vine->copulas[0][i] = dml_copula_alloc_normal(0.5);
            } else {
                vine->copulas[0][i] = dml_copula_alloc_indep();
            }
        }
    }
    vine_data = gsl_matrix_alloc(m, n);
    dml_vine_ran(vine, rng, vine_data);
    fitted = dml_vine_alloc(DML_VINE_CVINE, n);
    dml_vine_fit(fitted, vine_data, DML_VINE_WEIGHT_TAU,
                 DML_VINE_TRUNCATION_AIC, DML_COPULA_INDEPTEST_NONE, 1.0,
                 &types[0], types_size, DML_COPULA_SELECTION_AIC);
    fitted_data = gsl_matrix_alloc(m, n);
    dml_vine_ran(fitted, rng, fitted_data);
    xi = gsl_matrix_column(fitted_data, 0);
    for (size_t i = 1; i < n - 1; i++) {
        xj = gsl_matrix_column(fitted_data, i);
        corr = gsl_stats_correlation(xi.vector.data, xi.vector.stride,
                                     xj.vector.data, xj.vector.stride, m);
        if (i % 2 != 0) {
            g_assert(fabs(corr - 0.5) < 0.1);
        } else {
            g_assert(fabs(corr) < 0.1);
        }
    }

    dml_vine_free(vine);
    dml_vine_free(fitted);
    gsl_matrix_free(vine_data);
    gsl_matrix_free(fitted_data);
    gsl_rng_free(rng);
}

void
test_cvine_ran_fit_5d_normal_indep()
{
    size_t n = 5, m = 1000;
    gsl_rng *rng;
    dml_copula_type_t types[] = { DML_COPULA_NORMAL };
    size_t types_size = 1;
    dml_vine_t *vine, *fitted;
    gsl_matrix *vine_data, *fitted_data;
    gsl_vector_view xi, xj;
    double corr;

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, g_test_rand_int());

    // Initialize the original C-vine.
    vine = dml_vine_alloc(DML_VINE_CVINE, n);
    vine->trees = n - 1;
    vine->order[0] = 0;
    vine->order[1] = 1;
    vine->order[2] = 2;
    vine->order[3] = 3;
    vine->order[4] = 4;
    vine->copulas[0][0] = dml_copula_alloc_normal(0.75);
    vine->copulas[0][1] = dml_copula_alloc_normal(-0.75);
    vine->copulas[0][2] = dml_copula_alloc_normal(0.95);
    vine->copulas[0][3] = dml_copula_alloc_normal(-0.95);
    vine->copulas[1][0] = dml_copula_alloc_indep();
    vine->copulas[1][1] = dml_copula_alloc_indep();
    vine->copulas[1][2] = dml_copula_alloc_indep();
    vine->copulas[2][0] = dml_copula_alloc_indep();
    vine->copulas[2][1] = dml_copula_alloc_indep();
    vine->copulas[3][0] = dml_copula_alloc_indep();

    // Simulate the original vine and estimate a new vine from the sample.
    // Then, sample the fitted vine.
    vine_data = gsl_matrix_alloc(m, n);
    dml_vine_ran(vine, rng, vine_data);
    fitted = dml_vine_alloc(DML_VINE_CVINE, n);
    dml_vine_fit(fitted, vine_data, DML_VINE_WEIGHT_TAU,
                 DML_VINE_TRUNCATION_NONE, DML_COPULA_INDEPTEST_TAU, 0.05,
                 &types[0], types_size, DML_COPULA_SELECTION_AIC);
    fitted_data = gsl_matrix_alloc(m, n);
    dml_vine_ran(fitted, rng, fitted_data);

    // Check the correlations of the first tree.
    xi = gsl_matrix_column(fitted_data, 0);
    xj = gsl_matrix_column(fitted_data, 1);
    corr = gsl_stats_correlation(xi.vector.data, xi.vector.stride,
                                 xj.vector.data, xj.vector.stride, m);
    g_assert(fabs(corr - 0.75) < 0.1);
    xj = gsl_matrix_column(fitted_data, 2);
    corr = gsl_stats_correlation(xi.vector.data, xi.vector.stride,
                                 xj.vector.data, xj.vector.stride, m);
    g_assert(fabs(corr + 0.75) < 0.1);
    xj = gsl_matrix_column(fitted_data, 3);
    corr = gsl_stats_correlation(xi.vector.data, xi.vector.stride,
                                 xj.vector.data, xj.vector.stride, m);
    g_assert(fabs(corr - 0.95) < 0.1);
    xj = gsl_matrix_column(fitted_data, 4);
    corr = gsl_stats_correlation(xi.vector.data, xi.vector.stride,
                                 xj.vector.data, xj.vector.stride, m);
    g_assert(fabs(corr + 0.95) < 0.1);

    dml_vine_free(vine);
    dml_vine_free(fitted);
    gsl_matrix_free(vine_data);
    gsl_matrix_free(fitted_data);
    gsl_rng_free(rng);
}

void
test_dvine_alloc()
{
    dml_vine_t *vine;

    vine = dml_vine_alloc(DML_VINE_DVINE, 2);
    g_assert(dml_vine_type(vine) == DML_VINE_DVINE);
    dml_vine_free(vine);
}

void
test_rvine_alloc()
{
    dml_vine_t *vine;

    vine = dml_vine_alloc(DML_VINE_RVINE, 2);
    g_assert(dml_vine_type(vine) == DML_VINE_RVINE);
    dml_vine_free(vine);
}

void
test_rvine_ran_fit_3d_normal()
{
    size_t n = 3, m = 1000;
    gsl_rng *rng;
    dml_copula_type_t types[] = { DML_COPULA_NORMAL };
    size_t types_size = 1;
    dml_vine_t *vine, *fitted;
    gsl_matrix *vine_data, *fitted_data;
    gsl_vector_view vine_xi, vine_xj, fitted_xi, fitted_xj;
    double vine_corr, fitted_corr;

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, g_test_rand_int());

    // Initialize the original vine.
    vine = dml_vine_alloc(DML_VINE_RVINE, n);
    vine->trees = n - 1;
    vine->matrix[0][0] = 3;
    vine->matrix[1][0] = 1;
    vine->matrix[1][1] = 2;
    vine->matrix[2][0] = 2;
    vine->matrix[2][1] = 1;
    vine->matrix[2][2] = 1;
    vine->copulas[1][0] = dml_copula_alloc_normal(0.25);
    vine->copulas[2][0] = dml_copula_alloc_normal(0.9);
    vine->copulas[2][1] = dml_copula_alloc_normal(-0.75);
    vine->order[0] = 0;
    vine->order[1] = 1;
    vine->order[2] = 2;

    // Simulate the original vine and estimate a new vine from the sample.
    vine_data = gsl_matrix_alloc(m, n);
    dml_vine_ran(vine, rng, vine_data);
    fitted = dml_vine_alloc(DML_VINE_RVINE, n);
    dml_vine_fit(fitted, vine_data, DML_VINE_WEIGHT_TAU, DML_VINE_TRUNCATION_NONE,
             DML_COPULA_INDEPTEST_NONE, 0, &types[0], types_size,
             DML_COPULA_SELECTION_AIC);

    // Sample the fitted vine and check correlations.
    fitted_data = gsl_matrix_alloc(m, n);
    dml_vine_ran(fitted, rng, fitted_data);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < i; j++) {
            vine_xi = gsl_matrix_column(vine_data, i);
            vine_xj = gsl_matrix_column(vine_data, j);
            fitted_xi = gsl_matrix_column(fitted_data, i);
            fitted_xj = gsl_matrix_column(fitted_data, j);
            vine_corr = gsl_stats_correlation(vine_xi.vector.data,
                                              vine_xi.vector.stride,
                                              vine_xj.vector.data,
                                              vine_xj.vector.stride, m);
            fitted_corr = gsl_stats_correlation(fitted_xi.vector.data,
                                                fitted_xi.vector.stride,
                                                fitted_xj.vector.data,
                                                fitted_xj.vector.stride, m);
            g_assert(fabs(vine_corr - fitted_corr) < 0.1);
        }
    }

    dml_vine_free(vine);
    dml_vine_free(fitted);
    gsl_matrix_free(vine_data);
    gsl_matrix_free(fitted_data);
    gsl_rng_free(rng);
}

void
test_rvine_ran_fit_7d_normal()
{
    size_t n = 7, m = 5000;
    gsl_rng *rng;
    dml_copula_type_t types[] = { DML_COPULA_NORMAL };
    size_t types_size = 1;
    dml_vine_t *vine, *fitted;
    gsl_matrix *vine_data, *fitted_data;
    gsl_vector_view vine_xi, vine_xj, fitted_xi, fitted_xj;
    double vine_corr, fitted_corr;

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, g_test_rand_int());

    // Initialize the original vine.
    vine = dml_vine_alloc(DML_VINE_RVINE, n);
    vine->trees = n - 1;
    vine->matrix[0][0] = 7;
    vine->matrix[1][0] = 4;
    vine->matrix[1][1] = 6;
    vine->matrix[2][0] = 3;
    vine->matrix[2][1] = 4;
    vine->matrix[2][2] = 5;
    vine->matrix[3][0] = 6;
    vine->matrix[3][1] = 3;
    vine->matrix[3][2] = 4;
    vine->matrix[3][3] = 4;
    vine->matrix[4][0] = 5;
    vine->matrix[4][1] = 5;
    vine->matrix[4][2] = 3;
    vine->matrix[4][3] = 2;
    vine->matrix[4][4] = 3;
    vine->matrix[5][0] = 2;
    vine->matrix[5][1] = 1;
    vine->matrix[5][2] = 1;
    vine->matrix[5][3] = 1;
    vine->matrix[5][4] = 2;
    vine->matrix[5][5] = 2;
    vine->matrix[6][0] = 1;
    vine->matrix[6][1] = 2;
    vine->matrix[6][2] = 2;
    vine->matrix[6][3] = 3;
    vine->matrix[6][4] = 1;
    vine->matrix[6][5] = 1;
    vine->matrix[6][6] = 1;
    vine->copulas[1][0] = dml_copula_alloc_normal(0.2);
    vine->copulas[2][0] = dml_copula_alloc_normal(0.1);
    vine->copulas[2][1] = dml_copula_alloc_normal(-0.1);
    vine->copulas[3][0] = dml_copula_alloc_normal(0.1);
    vine->copulas[3][1] = dml_copula_alloc_normal(0.2);
    vine->copulas[3][2] = dml_copula_alloc_normal(-0.1);
    vine->copulas[4][0] = dml_copula_alloc_normal(0.2);
    vine->copulas[4][1] = dml_copula_alloc_normal(0.2);
    vine->copulas[4][2] = dml_copula_alloc_normal(-0.3);
    vine->copulas[4][3] = dml_copula_alloc_normal(-0.2);
    vine->copulas[5][0] = dml_copula_alloc_normal(0.2);
    vine->copulas[5][1] = dml_copula_alloc_normal(0.2);
    vine->copulas[5][2] = dml_copula_alloc_normal(-0.1);
    vine->copulas[5][3] = dml_copula_alloc_normal(-0.1);
    vine->copulas[5][4] = dml_copula_alloc_normal(0.1);
    vine->copulas[6][0] = dml_copula_alloc_normal(0.8);
    vine->copulas[6][1] = dml_copula_alloc_normal(0.3);
    vine->copulas[6][2] = dml_copula_alloc_normal(-0.7);
    vine->copulas[6][3] = dml_copula_alloc_normal(-0.6);
    vine->copulas[6][4] = dml_copula_alloc_normal(0.6);
    vine->copulas[6][5] = dml_copula_alloc_normal(0.9);
    vine->order[4] = 0;
    vine->order[1] = 1;
    vine->order[0] = 2;
    vine->order[6] = 3;
    vine->order[5] = 4;
    vine->order[2] = 5;
    vine->order[3] = 6;

    // Simulate the original vine and estimate a new vine from the sample.
    vine_data = gsl_matrix_alloc(m, n);
    dml_vine_ran(vine, rng, vine_data);
    fitted = dml_vine_alloc(DML_VINE_RVINE, n);
    dml_vine_fit(fitted, vine_data, DML_VINE_WEIGHT_TAU, DML_VINE_TRUNCATION_NONE,
             DML_COPULA_INDEPTEST_NONE, 0, &types[0], types_size,
             DML_COPULA_SELECTION_AIC);

    // Sample the fitted vine and check correlations.
    fitted_data = gsl_matrix_alloc(m, n);
    dml_vine_ran(fitted, rng, fitted_data);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < i; j++) {
            vine_xi = gsl_matrix_column(vine_data, i);
            vine_xj = gsl_matrix_column(vine_data, j);
            fitted_xi = gsl_matrix_column(fitted_data, i);
            fitted_xj = gsl_matrix_column(fitted_data, j);
            vine_corr = gsl_stats_correlation(vine_xi.vector.data,
                                              vine_xi.vector.stride,
                                              vine_xj.vector.data,
                                              vine_xj.vector.stride, m);
            fitted_corr = gsl_stats_correlation(fitted_xi.vector.data,
                                                fitted_xi.vector.stride,
                                                fitted_xj.vector.data,
                                                fitted_xj.vector.stride, m);
            g_assert(fabs(vine_corr - fitted_corr) < 0.25);
        }
    }

    dml_vine_free(vine);
    dml_vine_free(fitted);
    gsl_matrix_free(vine_data);
    gsl_matrix_free(fitted_data);
    gsl_rng_free(rng);
}

void test_rvine_ran_fit_9d_normal_trunc()
{
    size_t n = 7, m = 5000;
    gsl_rng *rng;
    dml_copula_type_t types[] = { DML_COPULA_NORMAL };
    size_t types_size = 1;
    dml_vine_t *vine, *fitted;
    gsl_matrix *vine_data, *fitted_data;
    gsl_vector_view vine_xi, vine_xj, fitted_xi, fitted_xj;
    double vine_corr, fitted_corr;

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, g_test_rand_int());

    // Initialize the original vine.
    vine = dml_vine_alloc(DML_VINE_RVINE, n);
    vine->trees = n - 1;
    vine->matrix[0][0] = 7;
    vine->matrix[1][0] = 4;
    vine->matrix[1][1] = 6;
    vine->matrix[2][0] = 3;
    vine->matrix[2][1] = 4;
    vine->matrix[2][2] = 5;
    vine->matrix[3][0] = 6;
    vine->matrix[3][1] = 3;
    vine->matrix[3][2] = 4;
    vine->matrix[3][3] = 4;
    vine->matrix[4][0] = 5;
    vine->matrix[4][1] = 5;
    vine->matrix[4][2] = 3;
    vine->matrix[4][3] = 2;
    vine->matrix[4][4] = 3;
    vine->matrix[5][0] = 2;
    vine->matrix[5][1] = 1;
    vine->matrix[5][2] = 1;
    vine->matrix[5][3] = 1;
    vine->matrix[5][4] = 2;
    vine->matrix[5][5] = 2;
    vine->matrix[6][0] = 1;
    vine->matrix[6][1] = 2;
    vine->matrix[6][2] = 2;
    vine->matrix[6][3] = 3;
    vine->matrix[6][4] = 1;
    vine->matrix[6][5] = 1;
    vine->matrix[6][6] = 1;
    vine->copulas[1][0] = dml_copula_alloc_normal(0.2);
    vine->copulas[2][0] = dml_copula_alloc_normal(0.1);
    vine->copulas[2][1] = dml_copula_alloc_normal(-0.1);
    vine->copulas[3][0] = dml_copula_alloc_normal(0.1);
    vine->copulas[3][1] = dml_copula_alloc_normal(0.2);
    vine->copulas[3][2] = dml_copula_alloc_normal(-0.1);
    vine->copulas[4][0] = dml_copula_alloc_normal(0.2);
    vine->copulas[4][1] = dml_copula_alloc_normal(0.2);
    vine->copulas[4][2] = dml_copula_alloc_normal(-0.3);
    vine->copulas[4][3] = dml_copula_alloc_normal(-0.2);
    vine->copulas[5][0] = dml_copula_alloc_normal(0.2);
    vine->copulas[5][1] = dml_copula_alloc_normal(0.2);
    vine->copulas[5][2] = dml_copula_alloc_normal(-0.1);
    vine->copulas[5][3] = dml_copula_alloc_normal(-0.1);
    vine->copulas[5][4] = dml_copula_alloc_normal(0.1);
    vine->copulas[6][0] = dml_copula_alloc_normal(0.8);
    vine->copulas[6][1] = dml_copula_alloc_normal(0.3);
    vine->copulas[6][2] = dml_copula_alloc_normal(-0.7);
    vine->copulas[6][3] = dml_copula_alloc_normal(-0.6);
    vine->copulas[6][4] = dml_copula_alloc_normal(0.6);
    vine->copulas[6][5] = dml_copula_alloc_normal(0.9);
    vine->order[4] = 0;
    vine->order[1] = 1;
    vine->order[0] = 2;
    vine->order[6] = 3;
    vine->order[5] = 4;
    vine->order[2] = 5;
    vine->order[3] = 6;

    // Simulate the original vine and estimate a new vine from the sample.
    vine_data = gsl_matrix_alloc(m, n);
    dml_vine_ran(vine, rng, vine_data);
    fitted = dml_vine_alloc(DML_VINE_RVINE, n);
    dml_vine_fit(fitted, vine_data, DML_VINE_WEIGHT_TAU, DML_VINE_TRUNCATION_AIC,
             DML_COPULA_INDEPTEST_NONE, 0, &types[0], types_size,
             DML_COPULA_SELECTION_AIC);

    // Sample the fitted vine and check correlations.
    fitted_data = gsl_matrix_alloc(m, n);
    dml_vine_ran(fitted, rng, fitted_data);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < i; j++) {
            vine_xi = gsl_matrix_column(vine_data, i);
            vine_xj = gsl_matrix_column(vine_data, j);
            fitted_xi = gsl_matrix_column(fitted_data, i);
            fitted_xj = gsl_matrix_column(fitted_data, j);
            vine_corr = gsl_stats_correlation(vine_xi.vector.data,
                                              vine_xi.vector.stride,
                                              vine_xj.vector.data,
                                              vine_xj.vector.stride, m);
            fitted_corr = gsl_stats_correlation(fitted_xi.vector.data,
                                                fitted_xi.vector.stride,
                                                fitted_xj.vector.data,
                                                fitted_xj.vector.stride, m);
            g_assert(fabs(vine_corr - fitted_corr) < 0.25);
        }
    }

    dml_vine_free(vine);
    dml_vine_free(fitted);
    gsl_matrix_free(vine_data);
    gsl_matrix_free(fitted_data);
    gsl_rng_free(rng);
}
