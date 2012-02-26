// DML - Dependence Modeling Library
// Copyright (C) 2011-2012 Yasser González-Fernández <ygonzalezfernandez@gmail.com>

#include "config.h"

#include "src/dml.h"
#include "tests/common.h"

void
test_measure_alloc()
{
    dml_measure_t *measure;

    measure = dml_measure_alloc(NULL, NULL);
    dml_measure_free(measure);
}

// Test using the data set presented in the following reference:
// Genest, C. and Favre, A. C. (2007) Everything you always want to known
// about copula modeling but were afraid to ask. Journal of Hydrologic
// Engineering, 12, 347-68.

void
test_measure_tau_small()
{
    dml_measure_t *measure;
    gsl_vector_view x_view, y_view;
    double x[] = { -2.224, -1.538, -0.807, 0.024, 0.052, 1.324 };
    double y[] = { 0.431, 1.035, 0.586, 1.465, 1.115, -0.847 };

    x_view = gsl_vector_view_array(x, 6);
    y_view = gsl_vector_view_array(y, 6);

    measure = dml_measure_alloc(&x_view.vector, &y_view.vector);
    g_assert(fabs(dml_measure_tau_coef(measure) - 0.06) <= 0.01);
    g_assert(fabs(dml_measure_tau_pvalue(measure) - 0.85) <= 0.01);

    dml_measure_free(measure);
}

void
test_measure_tau_large()
{
    FILE *f;
    char *path;
    gsl_matrix *data;
    gsl_vector_view x, y;
    dml_measure_t *measure;

    data = gsl_matrix_alloc(1000, 2);
    path = g_build_filename("tests", "data", "tau_large.dat", NULL);
    f = fopen(path, "r");
    gsl_matrix_fscanf(f, data);

    x = gsl_matrix_column(data, 0);
    y = gsl_matrix_column(data, 1);
    measure = dml_measure_alloc(&x.vector, &y.vector);
    g_assert(fabs(dml_measure_tau_coef(measure) - 0.03) <= 0.01);
    g_assert(fabs(dml_measure_tau_pvalue(measure) - 0.09) <= 0.01);

    dml_measure_free(measure);
    gsl_matrix_free(data);
    fclose(f);
    g_free(path);
}

void
test_measure_cvm_normal()
{
    size_t n = 150;
    gsl_rng *rng;
    gsl_vector *u, *v;
    dml_measure_t *measure;
    dml_copula_t *copula;

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, g_test_rand_int());

    u = gsl_vector_alloc(n);
    v = gsl_vector_alloc(n);

    copula = dml_copula_alloc_normal(0.75);
    dml_copula_ran(copula, rng, u, v);
    measure = dml_measure_alloc(u, v);
    g_assert(dml_measure_cvm_pvalue(measure, rng) < 0.01);

    dml_measure_free(measure);
    dml_copula_free(copula);
    gsl_vector_free(u);
    gsl_vector_free(v);
    gsl_rng_free(rng);
}

void
test_measure_cvm_indep()
{
    size_t n = 350;
    gsl_rng *rng;
    gsl_vector *u, *v;
    dml_measure_t *measure;
    dml_copula_t *copula;

    rng = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(rng, g_test_rand_int());

    u = gsl_vector_alloc(n);
    v = gsl_vector_alloc(n);

    copula = dml_copula_alloc_indep();
    dml_copula_ran(copula, rng, u, v);
    measure = dml_measure_alloc(u, v);
    g_assert(dml_measure_cvm_pvalue(measure, rng) > 0.01);

    dml_measure_free(measure);
    dml_copula_free(copula);
    gsl_vector_free(u);
    gsl_vector_free(v);
    gsl_rng_free(rng);
}
