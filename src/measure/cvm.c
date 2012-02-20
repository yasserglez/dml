// DML - Dependence Modeling Library
// Copyright (C) 2011-2012 Yasser González-Fernández <ygonzalezfernandez@gmail.com>

#include "config.h"

#include <stdbool.h>

#include <glib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

#include "src/dml.h"

// Global variables for the approximate independent realizations of the test
// statistics under mutual independence. This memory is not freed!

static size_t stats_sample_size = 0;
static size_t stats_count = 100;
static double *stats = NULL;

// Auxiliar functions from the copula R package.

static void
J_u(int n, int p, double *R, double *J)
{
    int i, j, l, m;

    m = 0;
    for (j = 0; j < p; j++)
        for (l = 0; l < n; l++)
            for (i = 0; i < n; i++)
                J[m++] = 1.0 - fmax(R[i + n * j], R[l + n * j]) / n;

}

static void
K_array(int n, int p, double *J, double *K)
{
    int i, j, l, m, n2 = n * n;

    m = 0;
    for (j = 0; j < p; j++)
        for (i = 0; i < n; i++) {
            K[m] = 0.0;
            for (l = 0; l < n; l++)
                K[m] += J[n2 * j + n * l + i];
            K[m++] /= (double) n;
        }
}

static void
L_array(int n, int p, double *K, double *L)
{
    int i, j;

    for (j = 0; j < p; j++) {
        L[j] = 0.0;
        for (i = 0; i < n; i++)
            L[j] += K[n * j + i];
        L[j] /= (double) n;
    }
}

static double
I_n(int n, int p, double *J, double *K, double *L)
{
    int i, j, l, n2 = n * n;
    double In, sum, prod, part1, part2, part3;

    /* First term */
    sum = 0.0;
    for (i = 0; i < n; i++)
        for (l = 0; l < n; l++) {
            prod = 1.0;
            for (j = 0; j < p; j++)
                prod *= J[n2 * j + n * l + i];
            sum += prod;
        }
    part1 = sum / (double) (n);

    /* Second term */
    sum = 0.0;
    for (i = 0; i < n; i++) {
        prod = 1.0;
        for (j = 0; j < p; j++)
            prod *= K[j * n + i]; /* K(i, j) */
        sum += prod;
    }
    part2 = 2.0 * sum;

    /* Third term */
    prod = 1.0;
    for (j = 0; j < p; j++)
        prod *= L[j];
    part3 = prod * (double) n;

    In = part1 - part2 + part3;

    return In;
}

static void
compute_cvm_stat(dml_measure_t *measure)
{
    size_t n;
    double *R, *J, *K, *L;
    gsl_permutation *x_rank, *y_rank;
    double stat;

    n = measure->x->size;

    R = g_malloc0_n(n * 2, sizeof(double));
    J = g_malloc0_n(n * n * 2, sizeof(double));
    K = g_malloc0_n(n * 2, sizeof(double));
    L = g_malloc0_n(2, sizeof(double));

    // Compute the ranks of the data.
    x_rank = dml_measure_x_rank(measure);
    y_rank = dml_measure_y_rank(measure);
    for (size_t i = 0; i < n; i++) {
        R[0 * n + i] = x_rank->data[i] + 1;
        R[1 * n + i] = y_rank->data[i] + 1;
    }

    // Compute arrays J, K, L.
    J_u(n, 2, R, J);
    K_array(n, 2, J, K);
    L_array(n, 2, K, L);

    // Compute the value of the global statistic.
    stat = I_n(n, 2, J, K, L);

    measure->cvm_stat = stat;

    g_free(R);
    g_free(J);
    g_free(K);
    g_free(L);
}

static void
compute_cvm_pvalue(dml_measure_t *measure, const gsl_rng *rng)
{
    size_t k, count;
    int i, j, index;
    double *R, *J, *K, *L;
    double r, stat, pvalue;

    if (stats_sample_size != measure->x->size) {
        // Update the realizations of the test statistics. This code is based
        // on the simulate_empirical_copula function of the copula R package.
        // It generates approximate independent realizations of the test
        // statistics under mutual independence. The original code of the
        // function was modified, since only the bivariate case is used here.

        // Update the sample size.
        stats_sample_size = measure->x->size;

        if (stats == NULL) {
            // First execution. This memory is not freed!
            stats = g_malloc0_n(stats_count, sizeof(double));
        }

        R = g_malloc0_n(stats_sample_size * 2, sizeof(double));
        J = g_malloc0_n(stats_sample_size * stats_sample_size * 2, sizeof(double));
        K = g_malloc0_n(stats_sample_size * 2, sizeof(double));
        L = g_malloc0_n(2, sizeof(double));

        for (k = 0; k < stats_count; k++) {
            /* Generate data */
            for (j = 0; j < 2; j++) {
                for (i = 0; i < stats_sample_size; i++)
                    R[i + stats_sample_size * j] = i + 1;

                /* Permutation = Random ranks in column j */
                for (i = stats_sample_size - 1; i >= 0; i--) {
                    r = R[j * stats_sample_size + i];
                    index = (int) ((i + 1) * gsl_rng_uniform(rng));
                    R[j * stats_sample_size + i] = R[j * stats_sample_size + index];
                    R[j * stats_sample_size + index] = r;
                }
            }

            /* Compute arrays J, K, L */
            J_u(stats_sample_size, 2, R, J);
            K_array(stats_sample_size, 2, J, K);
            L_array(stats_sample_size, 2, K, L);

            /* Global statistic under independence */
            stats[k] = I_n(stats_sample_size, 2, J, K, L);
        }

        g_free(R);
        g_free(J);
        g_free(K);
        g_free(L);
    }

    // Compute the corresponding p-value.
    stat = dml_measure_cvm_stat(measure);
    count = 0;
    for (k = 0; k < stats_count; k++)
        if (stats[k] >= stat)
            count++;
    pvalue = (double) (count + 0.5) / (stats_count + 1.0);

    measure->cvm_pvalue = pvalue;
}

double
dml_measure_cvm_stat(dml_measure_t *measure)
{
    if (gsl_isnan(measure->cvm_stat)) {
        compute_cvm_stat(measure);
    }

    return measure->cvm_stat;
}

double
dml_measure_cvm_pvalue(dml_measure_t *measure, const gsl_rng *rng)
{
    if (gsl_isnan(measure->cvm_pvalue)) {
        compute_cvm_pvalue(measure, rng);
    }

    return measure->cvm_pvalue;
}
