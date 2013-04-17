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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

#include "src/dml.h"

inline dml_copula_t *
dml_copula_alloc(const dml_copula_type_t type)
{
    dml_copula_t *copula;

    switch (type) {
    case DML_COPULA_INDEP:
        copula = dml_copula_alloc_indep();
        break;
    case DML_COPULA_NORMAL:
        copula = dml_copula_alloc_normal(0);
        break;
    case DML_COPULA_CLAYTON:
        copula = dml_copula_alloc_clayton(0);
        break;
    case DML_COPULA_RCLAYTON90:
        copula = dml_copula_alloc_rclayton90(0);
        break;
    case DML_COPULA_RCLAYTON180:
        copula = dml_copula_alloc_rclayton180(0);
        break;
    case DML_COPULA_RCLAYTON270:
        copula = dml_copula_alloc_rclayton270(0);
        break;
    default:
        copula = NULL;
        break;
    }

    return copula;
}

inline dml_copula_type_t
dml_copula_type(const dml_copula_t *copula)
{
    return copula->type;
}

inline void
dml_copula_fit(dml_copula_t *copula,
               const gsl_vector *u,
               const gsl_vector *v,
               dml_measure_t *measure)
{
    if (copula->fit != NULL) {
        copula->fit(copula, u, v, measure);
    }
}

void
dml_copula_ran(const dml_copula_t *copula,
               const gsl_rng *rng,
               gsl_vector *u,
               gsl_vector *v)
{
    gsl_vector *tmp;

    tmp = gsl_vector_alloc(u->size);
    for (size_t i = 0; i < u->size; i++) {
        gsl_vector_set(v, i, gsl_rng_uniform(rng));
        gsl_vector_set(tmp, i, gsl_rng_uniform(rng));
    }
    dml_copula_hinv(copula, tmp, v, u);

    gsl_vector_free(tmp);
}

inline void
dml_copula_pdf(const dml_copula_t *copula,
               const gsl_vector *u,
               const gsl_vector *v,
               gsl_vector *pdf)
{
    copula->pdf(copula, u, v, pdf);
}

inline void
dml_copula_cdf(const dml_copula_t *copula,
               const gsl_vector *u,
               const gsl_vector *v,
               gsl_vector *cdf)
{
    copula->cdf(copula, u, v, cdf);
}

inline void
dml_copula_h(const dml_copula_t *copula,
             const gsl_vector *u,
             const gsl_vector *v,
             gsl_vector *h)
{
    copula->h(copula, u, v, h);
}

inline void
dml_copula_hinv(const dml_copula_t *copula,
                const gsl_vector *u,
                const gsl_vector *v,
                gsl_vector *hinv)
{
    copula->hinv(copula, u, v, hinv);
}

inline void
dml_copula_aic(const dml_copula_t *copula,
               const gsl_vector *u,
               const gsl_vector *v,
               double *aic)
{
    copula->aic(copula, u, v, aic);
}

// Two functions from the copula R package that implement the normal
// copula goodness-of-fit test based on a Cramer-von Mises statistic.
// Only the value of the statistic is calculated here and it is used
// in the copula selection procedure of the vines.

static double
empcop(int n, int p, double *U, double *V, int m, int k)
{
    int i, j, ind;
    double ec = 0.0;

    for (i = 0; i < n; i++) {
        ind = 1;
        for (j = 0; j < p; j++)
            ind *= (U[i + n * j] <= V[k + m * j]);
        ec += (double) ind;
    }
    return ec / (double) n;
}

static void
cramer_vonMises_2(int *p,
                  double *U,
                  int *n,
                  double *V,
                  int *m,
                  double *Ctheta,
                  double *stat)
{
    int i;
    double s = 0.0, diff;

    for (i = 0; i < *m; i++) {
        diff = empcop(*n, *p, U, V, *m, i) - Ctheta[i];
        s += diff * diff;
    }

    *stat = s * (*n) / (*m);
}

void
dml_copula_cvm_stat(const dml_copula_t *copula,
                    const gsl_vector *u,
                    const gsl_vector *v,
                    dml_measure_t *measure,
                    const gsl_rng *rng,
                    double *cvm_stat)
{
    int n, p = 2;
    gsl_vector *cdf;
    double *U, *Ctheta;
    double stat;
    gsl_vector *u_pseudo, *v_pseudo;
    bool measure_dealloc;
    gsl_permutation *u_rank, *v_rank;

    if (measure == NULL) {
        measure = dml_measure_alloc(u, v);
        measure_dealloc = true;
    } else {
        measure_dealloc = false;
    }

    n = (int) u->size;
    u_pseudo = gsl_vector_alloc(n);
    v_pseudo = gsl_vector_alloc(n);
    cdf = gsl_vector_alloc(n);
    U = g_malloc_n(2 * n, sizeof(double));
    Ctheta = g_malloc_n(n, sizeof(double));

    // Compute the ranks of the data.
    u_rank = dml_measure_x_rank(measure);
    v_rank = dml_measure_y_rank(measure);
    for (size_t i = 0; i < n; i++) {
        U[i + n*0] = (double) (u_rank->data[i] + 1.0) / (n + 1.0);
        gsl_vector_set(u_pseudo, i, U[i + n*0]);
        U[i + n*1] = (double) (v_rank->data[i] + 1.0) / (n + 1.0);
        gsl_vector_set(v_pseudo, i, U[i + n*1]);
    }

    dml_copula_cdf(copula, u_pseudo, v_pseudo, cdf);
    for (size_t i = 0; i < n; i++) {
        Ctheta[i] = gsl_vector_get(cdf, i);
    }

    cramer_vonMises_2(&p, U, &n, U, &n, Ctheta, &stat);

    *cvm_stat = stat;

    gsl_vector_free(cdf);
    gsl_vector_free(u_pseudo);
    gsl_vector_free(v_pseudo);
    g_free(U);
    g_free(Ctheta);
    if (measure_dealloc) dml_measure_free(measure);
}

inline void
dml_copula_free(dml_copula_t *copula)
{
    if (copula->free != NULL) {
        copula->free(copula);
    }
    g_free(copula);
}
