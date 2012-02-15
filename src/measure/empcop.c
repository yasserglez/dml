// DML - Dependence Modeling Library
// Copyright (C) 2011-2012 Yasser González-Fernández <ygonzalezfernandez@gmail.com>

#include "config.h"

#include <stdbool.h>

#include <glib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>

#include "src/dml.h"

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
compute_empcop_cvm_stat(dml_measure_t *measure)
{
    measure->empcop_cvm_stat = 0.72806;
}

static void
compute_empcop_cvm_pvalue(dml_measure_t *measure,
                          size_t num_stats,
                          double *stats)
{
    size_t k, count;
    double stat, pvalue;

    stat = dml_measure_empcop_cvm_stat(measure);
    count = 0;
    for (k = 0; k < num_stats; k++)
        if (stats[k] >= stat)
            count++;
    pvalue = (double) (count + 0.5) / (num_stats + 1.0);

    measure->empcop_cvm_pvalue = pvalue;
}

// Based on the simulate_empirical_copula function of the copula R package.
// This function generates approximate independent realizations of the test
// statistics under mutual independence. The original code of the function
// was modified, since only the bivariate case is used here. The parameters
// are: n (sample size), N (number of repetitions), rng (GSL random number
// generator) and stats (values of the global statistic under independence).
void
dml_measure_empcop_cvm_sim(size_t n,
                           const gsl_rng *rng,
                           size_t num_stats,
                           double *stats)
{
    int i, j, k, index;
    double *R = g_malloc0_n(n * 2, sizeof(double));
    double *J = g_malloc0_n(n * n * 2, sizeof(double));
    double *K = g_malloc0_n(n * 2, sizeof(double));
    double *L = g_malloc0_n(2, sizeof(double));
    double r;

    /* N repetitions */
    for (k = 0; k < num_stats; k++) {
        /* Generate data */
        for (j = 0; j < 2; j++) {
            for (i = 0; i < n; i++)
                R[i + n * j] = i + 1;

            /* Permutation = Random ranks in column j */
            for (i = n - 1; i >= 0; i--) {
                r = R[j * n + i];
                index = (int) ((i + 1) * gsl_rng_uniform(rng));
                R[j * n + i] = R[j * n + index];
                R[j * n + index] = r;
            }
        }

        /* Compute arrays J, K, L */
        J_u(n, 2, R, J);
        K_array(n, 2, J, K);
        L_array(n, 2, K, L);

        /* Global statistic under independence */
        stats[k] = I_n(n, 2, J, K, L);
    }

    g_free(R);
    g_free(J);
    g_free(K);
    g_free(L);
}

double
dml_measure_empcop_cvm_stat(dml_measure_t *measure)
{
    if (gsl_isnan(measure->empcop_cvm_stat)) {
        compute_empcop_cvm_stat(measure);
    }

    return measure->empcop_cvm_stat;
}

double
dml_measure_empcop_cvm_pvalue(dml_measure_t *measure,
                              size_t num_stats,
                              double *stats)
{
    if (gsl_isnan(measure->empcop_cvm_pvalue)) {
        compute_empcop_cvm_pvalue(measure, num_stats, stats);
    }

    return measure->empcop_cvm_pvalue;
}
