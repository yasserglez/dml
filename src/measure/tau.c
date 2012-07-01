/*
 * DML - Dependence Modeling Library
 * Copyright (C) 2011-2012 Yasser González Fernández <ygonzalezfernandez@gmail.com>
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

#include <stdbool.h>

#include <glib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

#include "src/dml.h"

// Based on the ktau function of the CDVine R package. This function
// computes Kendall's tau rank correlation coefficient in O(n log n) by means
// of the algorithm presented in Knight, W. R. (1966). A computer method for
// calculating Kendall's tau with ungrouped data. Journal of the American
// Statistical Association 61, 436-439.

static void
compute_tau_coef(dml_measure_t *measure)
{
    double *x, *y, *x_aux, *y_aux, *tmp;
    size_t i, j, k, l, m, n;
    size_t i_end, j_end;
    bool i_flag, j_flag, x_flag;
    double score = 0, denom = 0;
    size_t t = 0, u = 0, v = 0;
    gsl_permutation *x_rank;

    n = measure->x->size;
    x = g_malloc_n(n, sizeof(double));
    y = g_malloc_n(n, sizeof(double));
    x_aux = g_malloc_n(n, sizeof(double));
    y_aux = g_malloc_n(n, sizeof(double));

    // 1.1. Sort x and y in x order. This step differs with the original
    // algorithm to reuse the ranks of the data if already computed.
    x_rank = dml_measure_x_rank(measure);
    for (i = 0; i < n; i++) {
        x[x_rank->data[i]] = gsl_vector_get(measure->x, i);
        y[x_rank->data[i]] = gsl_vector_get(measure->y, i);
    }

    // 1.2. Count pairs of tied x's in t.
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

    // 2.1. Sort y again and count exchanges in score.
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

    // 2.2. Count pairs of tied y's in u.
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

    // 3. Calculate Kendall's score and denominator.
    denom = 0.5 * n * (n - 1);
    score = denom - (2.0 * score + t + u - v);
    denom = sqrt((denom - t) * (denom - u));
    measure->tau_coef = score / denom;

    g_free(x);
    g_free(y);
    g_free(x_aux);
    g_free(y_aux);
}

static void
compute_tau_pvalue(dml_measure_t *measure)
{
    size_t n;
    double x, coef;

    n = measure->x->size;
    coef = dml_measure_tau_coef(measure);
    x = sqrt((9 * n * (n - 1)) / (2 * (2 * n + 5))) * fabs(coef);
    measure->tau_pvalue = 2 * gsl_cdf_ugaussian_Q(x);
}

double
dml_measure_tau_coef(dml_measure_t *measure)
{
    if (gsl_isnan(measure->tau_coef)) {
        compute_tau_coef(measure);
    }

    return measure->tau_coef;
}

double
dml_measure_tau_pvalue(dml_measure_t *measure)
{
    if (gsl_isnan(measure->tau_pvalue)) {
        compute_tau_pvalue(measure);
    }

    return measure->tau_pvalue;
}
