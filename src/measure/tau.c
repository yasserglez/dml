// DML - Dependence Modeling Library
// Copyright (C) 2011-2012 Yasser González-Fernández <ygonzalezfernandez@gmail.com>

#include "config.h"

#include <stdbool.h>

#include <glib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cdf.h>

#include "src/dml.h"

/* Based on the ktau function of the CDVine R package. This function
 * computes Kendall's tau rank correlation coefficient in O(n log n) by means
 * of the algorithm presented in Knight, W. R. (1966). A computer method for
 * calculating Kendall's tau with ungrouped data. Journal of the American
 * Statistical Association 61, 436-439.
 */
static void
compute_tau_coef(dml_measure_tau_t *tau)
{
    double *x, *y, *x_aux, *y_aux, *tmp;
    size_t i, j, k, l, m, n;
    size_t i_end, j_end;
    bool i_flag, j_flag, x_flag;
    double score = 0, denom = 0;
    size_t t = 0, u = 0, v = 0;

    n = tau->x->size;
    x = g_malloc_n(n, sizeof(double));
    y = g_malloc_n(n, sizeof(double));
    for (i = 0; i < n; i++) {
        x[i] = gsl_vector_get(tau->x, i);
        y[i] = gsl_vector_get(tau->y, i);
    }
    x_aux = g_malloc_n(n, sizeof(double));
    y_aux = g_malloc_n(n, sizeof(double));

    /* 1.1. Sort x and y in x order. */
    k = 1;
    do {
        l = 0;
        do {
            i = l;
            j = (i + k) < n ? (i + k) : n;
            i_end = j;
            j_end = (j + k) < n ? (j + k) : n;
            do {
                i_flag = (i < i_end);
                j_flag = (j < j_end);
                x_flag = ((x[i] > x[j]) | ((x[i] == x[j]) & (y[i] > y[j])));
                if ((i_flag & !j_flag) | (i_flag & j_flag & !x_flag)) {
                    x_aux[l] = x[i]; y_aux[l] = y[i];
                    i++; l++;
                }
                if ((!i_flag & j_flag) | (i_flag & j_flag & x_flag)) {
                    x_aux[l] = x[j]; y_aux[l] = y[j];
                    j++; l++;
                }
            } while (i_flag | j_flag);
        } while (l < n);

        tmp = x; x = x_aux; x_aux = tmp;
        tmp = y; y = y_aux; y_aux = tmp;
        k *= 2;
    } while (k < n);

    /* 1.2. Count pairs of tied x's in t. */
    j = 1;
    m = 1;
    for (i = 1; i < n; i++)
        if (x[i] == x[i - 1]) {
            j++;
            if (y[i] == y[i - 1]) m++;
        } else if (j > 1) {
            t += j * (j - 1) / 2;
            if (m > 1) v += m * (m - 1) / 2;
            j = 1;
            m = 1;
        }
    t += j * (j - 1) / 2;
    v += m * (m - 1) / 2;

    /* 2.1. Sort y again and count exchanges in score. */
    k = 1;
    do {
        l = 0;
        do {
            i = l;
            j = (i + k) < n ? (i + k) : n;
            i_end = j;
            j_end = (j + k) < n ? (j + k) : n;
            do {
                i_flag = (i < i_end);
                j_flag = (j < j_end);
                x_flag = (y[i] > y[j]);
                if ((i_flag & !j_flag) | (i_flag & j_flag & !x_flag)) {
                    x_aux[l] = x[i]; y_aux[l] = y[i];
                    i++; l++;
                }
                if ((!i_flag & j_flag) | (i_flag & j_flag & x_flag)) {
                    x_aux[l] = x[j]; y_aux[l] = y[j];
                    score += i_end - i;
                    j++; l++;
                }
            } while (i_flag | j_flag);
        } while (l < n);

        tmp = x; x = x_aux; x_aux = tmp;
        tmp = y; y = y_aux; y_aux = tmp;
        k *= 2;
    } while (k < n);

    /* 2.2. Count pairs of tied y's in u. */
    j = 1;
    for (i = 1; i < n; i++) {
        if (y[i] == y[i - 1]) {
            j++;
        } else if (j > 1) {
            u += j * (j - 1) / 2;
            j = 1;
        }
    }
    u += j * (j - 1) / 2;

    /* 3. Calculate Kendall's score and denominator. */
    denom = 0.5 * n * (n - 1);
    score = denom - (2.0 * score + t + u - v);
    denom = sqrt((denom - t) * (denom - u));
    tau->coef = score / denom;

    g_free(x); g_free(y);
    g_free(x_aux); g_free(y_aux);
}

static void
compute_tau_pvalue(dml_measure_tau_t *tau)
{
    size_t n;
    double x, coef;

    n = tau->x->size;
    coef = dml_measure_tau_coef(tau);
    x = sqrt((9 * n * (n - 1)) / (2 * (2 * n + 5))) * fabs(coef);
    tau->pvalue = 2 * gsl_cdf_ugaussian_Q(x);
}

dml_measure_tau_t *
dml_measure_tau_alloc(const gsl_vector *x, const gsl_vector *y)
{
    dml_measure_tau_t *tau;

    tau = g_malloc(sizeof(dml_measure_tau_t));
    tau->x = x;
    tau->y = y;
    tau->coef = GSL_NAN;
    tau->pvalue = GSL_NAN;

    return tau;
}

inline double
dml_measure_tau_coef(dml_measure_tau_t *tau)
{
    if (gsl_isnan(tau->coef)) {
        compute_tau_coef(tau);
    }

    return tau->coef;
}

inline double
dml_measure_tau_pvalue(dml_measure_tau_t *tau)
{
    if (gsl_isnan(tau->pvalue)) {
        compute_tau_pvalue(tau);
    }

    return tau->pvalue;
}

inline void
dml_measure_tau_free(dml_measure_tau_t *tau)
{
    g_free(tau);
}
