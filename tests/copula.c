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

static double u[] = {
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
};

static double v[] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
    0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
    0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
    0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
    0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
    0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7,
    0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
    0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
};

static const double normal_rho[] = { -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4,
    -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

static const double clayton_theta_positive[] = { 0, 0.5, 1, 1.5, 2, 2.5, 3,
    3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50 };

static const double clayton_theta_negative[] = { 0, -0.5, -1, -1.5, -2, -2.5,
    -3, -3.5, -4, -4.5, -5, -6, -7, -8, -9, -10, -20, -30, -40, -50 };


#define TEST_DML_COPULA_ALLOC(c, C) \
    void \
    test_##c##_alloc() \
    { \
        dml_copula_t *copula; \
        copula = dml_copula_alloc(DML_COPULA_##C); \
        g_assert(dml_copula_type(copula) == DML_COPULA_##C); \
        dml_copula_free(copula); \
    }

#define TEST_DML_COPULA_FIT(c) \
    void \
    test_##c##_fit() \
    { \
        size_t n; \
        gsl_rng *rng; \
        gsl_vector *u, *v; \
        gsl_vector *pdf_copula, *pdf_fitted; \
        dml_copula_t *copula, *fitted; \
        n = 1000; \
        rng = gsl_rng_alloc(gsl_rng_taus); \
        gsl_rng_set(rng, g_test_rand_int()); \
        u = gsl_vector_alloc(n); \
        v = gsl_vector_alloc(n); \
        pdf_copula = gsl_vector_alloc(n); \
        pdf_fitted = gsl_vector_alloc(n); \
        copula = dml_copula_alloc_##c(); \
        fitted = dml_copula_alloc(dml_copula_type(copula)); \
        dml_copula_ran(copula, rng, u, v); \
        dml_copula_fit(fitted, u, v, NULL); \
        dml_copula_pdf(copula, u, v, pdf_copula); \
        dml_copula_pdf(fitted, u, v, pdf_fitted); \
        test_vectors_equal(pdf_copula, pdf_fitted, 0.25, 0.75); \
        dml_copula_free(copula); \
        dml_copula_free(fitted); \
        gsl_rng_free(rng); \
        gsl_vector_free(u); \
        gsl_vector_free(v); \
        gsl_vector_free(pdf_copula); \
        gsl_vector_free(pdf_fitted); \
    }

#define TEST_DML_COPULA_FIT_P(c, P, p) \
    void \
    test_##c##_fit() \
    { \
        size_t n; \
        gsl_rng *rng; \
        gsl_vector *u, *v; \
        gsl_vector *pdf_copula, *pdf_fitted; \
        dml_copula_t *copula, *fitted; \
        n = 1000; \
        rng = gsl_rng_alloc(gsl_rng_taus); \
        gsl_rng_set(rng, g_test_rand_int()); \
        u = gsl_vector_alloc(n); \
        v = gsl_vector_alloc(n); \
        pdf_copula = gsl_vector_alloc(n); \
        pdf_fitted = gsl_vector_alloc(n); \
        for (int i = 0; i < p; i++) { \
            copula = dml_copula_alloc_##c(P[i]); \
            fitted = dml_copula_alloc(dml_copula_type(copula)); \
            dml_copula_ran(copula, rng, u, v); \
            dml_copula_fit(fitted, u, v, NULL); \
            dml_copula_pdf(copula, u, v, pdf_copula); \
            dml_copula_pdf(fitted, u, v, pdf_fitted); \
            test_vectors_equal(pdf_copula, pdf_fitted, 0.25, 0.75); \
            dml_copula_free(copula); \
            dml_copula_free(fitted); \
        } \
        gsl_rng_free(rng); \
        gsl_vector_free(u); \
        gsl_vector_free(v); \
        gsl_vector_free(pdf_copula); \
        gsl_vector_free(pdf_fitted); \
    }

#define TEST_DML_COPULA_FUNC(c, f) \
    void \
    test_##c##_##f() \
    { \
        FILE *f; \
        char *path; \
        gsl_vector_view u_view, v_view; \
        gsl_vector *dat, *result; \
        dml_copula_t *copula; \
        dat = gsl_vector_alloc(121); \
        result = gsl_vector_alloc(121); \
        path = g_build_filename("tests", "data", "copula_" #c "_" #f ".dat", NULL); \
        f = fopen(path, "r"); \
        gsl_vector_fscanf(f, dat); \
        copula = dml_copula_alloc_##c(); \
        u_view = gsl_vector_view_array(u, 121); \
        v_view = gsl_vector_view_array(v, 121); \
        dml_copula_##f(copula, &u_view.vector, &v_view.vector, result); \
        test_vectors_equal(dat, result, 0.25, 1.0); \
        dml_copula_free(copula); \
        gsl_vector_free(dat); \
        gsl_vector_free(result); \
        g_free(path); \
        fclose(f); \
    }

#define TEST_DML_COPULA_FUNC_P(c, f, P, p) \
    void \
    test_##c##_##f() \
    { \
        FILE *f; \
        char *path; \
        gsl_vector_view u_view, v_view; \
        gsl_vector *dat, *result; \
        dml_copula_t *copula; \
        dat = gsl_vector_alloc(121); \
        result = gsl_vector_alloc(121); \
        path = g_build_filename("tests", "data", "copula_" #c "_" #f ".dat", NULL); \
        f = fopen(path, "r"); \
        for (int i = 0; i < p; i++) { \
            gsl_vector_fscanf(f, dat); \
            copula = dml_copula_alloc_##c(P[i]); \
            u_view = gsl_vector_view_array(u, 121); \
            v_view = gsl_vector_view_array(v, 121); \
            dml_copula_##f(copula, &u_view.vector, &v_view.vector, result); \
            test_vectors_equal(dat, result, 0.25, 1.0); \
            dml_copula_free(copula); \
        } \
        gsl_vector_free(dat); \
        gsl_vector_free(result); \
        g_free(path); \
        fclose(f); \
    }

TEST_DML_COPULA_ALLOC(indep, INDEP)
TEST_DML_COPULA_FIT(indep)
TEST_DML_COPULA_FUNC(indep, pdf)
TEST_DML_COPULA_FUNC(indep, cdf)
TEST_DML_COPULA_FUNC(indep, h)
TEST_DML_COPULA_FUNC(indep, hinv)

TEST_DML_COPULA_ALLOC(normal, NORMAL)
TEST_DML_COPULA_FIT_P(normal, normal_rho, 21)
TEST_DML_COPULA_FUNC_P(normal, pdf, normal_rho, 21)
TEST_DML_COPULA_FUNC_P(normal, cdf, normal_rho, 21)
TEST_DML_COPULA_FUNC_P(normal, h, normal_rho, 21)
TEST_DML_COPULA_FUNC_P(normal, hinv, normal_rho, 21)

TEST_DML_COPULA_ALLOC(clayton, CLAYTON)
TEST_DML_COPULA_FIT_P(clayton, clayton_theta_positive, 20)
TEST_DML_COPULA_FUNC_P(clayton, pdf, clayton_theta_positive, 20)
TEST_DML_COPULA_FUNC_P(clayton, cdf, clayton_theta_positive, 20)
TEST_DML_COPULA_FUNC_P(clayton, h, clayton_theta_positive, 20)
TEST_DML_COPULA_FUNC_P(clayton, hinv, clayton_theta_positive, 20)

TEST_DML_COPULA_ALLOC(rclayton90, RCLAYTON90)
TEST_DML_COPULA_FIT_P(rclayton90, clayton_theta_negative, 20)
TEST_DML_COPULA_FUNC_P(rclayton90, pdf, clayton_theta_negative, 20)
TEST_DML_COPULA_FUNC_P(rclayton90, cdf, clayton_theta_negative, 20)
TEST_DML_COPULA_FUNC_P(rclayton90, h, clayton_theta_negative, 20)
TEST_DML_COPULA_FUNC_P(rclayton90, hinv, clayton_theta_negative, 20)

TEST_DML_COPULA_ALLOC(rclayton180, RCLAYTON180)
TEST_DML_COPULA_FIT_P(rclayton180, clayton_theta_positive, 20)
TEST_DML_COPULA_FUNC_P(rclayton180, pdf, clayton_theta_positive, 20)
TEST_DML_COPULA_FUNC_P(rclayton180, cdf, clayton_theta_positive, 20)
TEST_DML_COPULA_FUNC_P(rclayton180, h, clayton_theta_positive, 20)
TEST_DML_COPULA_FUNC_P(rclayton180, hinv, clayton_theta_positive, 20)

TEST_DML_COPULA_ALLOC(rclayton270, RCLAYTON270)
TEST_DML_COPULA_FIT_P(rclayton270, clayton_theta_negative, 20)
TEST_DML_COPULA_FUNC_P(rclayton270, pdf, clayton_theta_negative, 20)
TEST_DML_COPULA_FUNC_P(rclayton270, cdf, clayton_theta_negative, 20)
TEST_DML_COPULA_FUNC_P(rclayton270, h, clayton_theta_negative, 20)
TEST_DML_COPULA_FUNC_P(rclayton270, hinv, clayton_theta_negative, 20)

void
test_copula_select_indeptest_none()
{
    size_t m = 1000;
    gsl_rng *rng;
    gsl_vector *u, *v;
    dml_copula_t *copula, *selected;
    dml_copula_type_t types[] = { DML_COPULA_NORMAL };
    size_t types_size = 1;

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, g_test_rand_int());

    u = gsl_vector_alloc(m);
    v = gsl_vector_alloc(m);

    // Independence.
    copula = dml_copula_alloc_indep();
    dml_copula_ran(copula, rng, u, v);
    selected = dml_copula_select(u, v, NULL, DML_COPULA_INDEPTEST_NONE, 1,
                                    &types[0], types_size,
                                    DML_COPULA_SELECTION_AIC);
    g_assert(DML_COPULA_NORMAL == dml_copula_type(selected));
    dml_copula_free(selected);
    dml_copula_free(copula);

    // Normal.
    copula = dml_copula_alloc_normal(0.75);
    dml_copula_ran(copula, rng, u, v);
    selected = dml_copula_select(u, v, NULL, DML_COPULA_INDEPTEST_NONE, 1,
                                    &types[0], types_size,
                                    DML_COPULA_SELECTION_AIC);
    g_assert(DML_COPULA_NORMAL == dml_copula_type(selected));
    dml_copula_free(selected);
    dml_copula_free(copula);

    gsl_vector_free(u);
    gsl_vector_free(v);
}

void
test_copula_select_indeptest_tau()
{
    size_t m = 1000;
    gsl_rng *rng;
    gsl_vector *u, *v;
    dml_copula_t *copula, *selected;
    dml_copula_type_t types[] = { DML_COPULA_NORMAL };
    size_t types_size = 1;

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, g_test_rand_int());

    u = gsl_vector_alloc(m);
    v = gsl_vector_alloc(m);

    // Independence.
    copula = dml_copula_alloc_indep();
    dml_copula_ran(copula, rng, u, v);
    selected = dml_copula_select(u, v, NULL, DML_COPULA_INDEPTEST_TAU, 0.01,
                                    &types[0], types_size,
                                    DML_COPULA_SELECTION_AIC);
    g_assert(dml_copula_type(copula) == dml_copula_type(selected));
    dml_copula_free(selected);
    dml_copula_free(copula);

    // Normal.
    copula = dml_copula_alloc_normal(0.75);
    dml_copula_ran(copula, rng, u, v);
    selected = dml_copula_select(u, v, NULL, DML_COPULA_INDEPTEST_TAU, 0.01,
                                    &types[0], types_size,
                                    DML_COPULA_SELECTION_AIC);
    g_assert(dml_copula_type(copula) == dml_copula_type(selected));
    dml_copula_free(selected);
    dml_copula_free(copula);

    gsl_vector_free(u);
    gsl_vector_free(v);
}

void
test_copula_select_aic()
{
    size_t m = 1000;
    gsl_rng *rng;
    gsl_vector *u, *v;
    dml_copula_t *copula, *selected;
    dml_copula_t *copulas[] = { dml_copula_alloc_normal(0.75),
        dml_copula_alloc_clayton(10), dml_copula_alloc_rclayton90(-10),
        dml_copula_alloc_rclayton180(10), dml_copula_alloc_rclayton270(-10) };
    dml_copula_type_t types[] = { DML_COPULA_NORMAL, DML_COPULA_CLAYTON,
        DML_COPULA_RCLAYTON90, DML_COPULA_RCLAYTON180, DML_COPULA_RCLAYTON270, };
    size_t num_copulas = 5;

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, g_test_rand_int());

    u = gsl_vector_alloc(m);
    v = gsl_vector_alloc(m);

    for (size_t i = 0; i < num_copulas; i++) {
        copula = copulas[i];
        dml_copula_ran(copula, rng, u, v);
        selected = dml_copula_select(u, v, NULL, DML_COPULA_INDEPTEST_NONE, 1,
                                        &types[0], num_copulas,
                                        DML_COPULA_SELECTION_AIC);
        g_assert(dml_copula_type(copula) == dml_copula_type(selected));
        dml_copula_free(selected);
        dml_copula_free(copula);
    }

    gsl_vector_free(u);
    gsl_vector_free(v);
}
