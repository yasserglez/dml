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

#include <glib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

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
               dml_measure_tau_t *tau)
{
    if (copula->fit != NULL) {
        copula->fit(copula, u, v, tau);
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

inline void
dml_copula_free(dml_copula_t *copula)
{
    if (copula->free != NULL) {
        copula->free(copula);
    }
    g_free(copula);
}
