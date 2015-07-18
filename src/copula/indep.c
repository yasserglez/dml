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
#include <gsl/gsl_vector.h>

#include "src/dml.h"

static void
copula_pdf_indep(const dml_copula_t *copula,
                 const gsl_vector *u,
                 const gsl_vector *v,
                 gsl_vector *pdf)
{
    gsl_vector_set_all(pdf, 1);
}

static void
copula_cdf_indep(const dml_copula_t *copula,
                 const gsl_vector *u,
                 const gsl_vector *v,
                 gsl_vector *cdf)
{
    gsl_vector_memcpy(cdf, u);
    gsl_vector_mul(cdf, v);
}

static void
copula_h_indep(const dml_copula_t *copula,
               const gsl_vector *u,
               const gsl_vector *v,
               gsl_vector *h)
{
    gsl_vector_memcpy(h, u);
}

static void
copula_hinv_indep(const dml_copula_t *copula,
                  const gsl_vector *u,
                  const gsl_vector *v,
                  gsl_vector *hinv)
{
    gsl_vector_memcpy(hinv, u);
}

static void
copula_aic_indep(const dml_copula_t *copula,
                 const gsl_vector *u,
                 const gsl_vector *v,
                 double *aic)
{
    *aic = 0;
}

dml_copula_t *
dml_copula_alloc_indep()
{
    dml_copula_t *copula;

    copula = g_malloc(sizeof(dml_copula_t));
    copula->type = DML_COPULA_INDEP;
    copula->fit = NULL; // Disabled.
    copula->pdf = copula_pdf_indep;
    copula->cdf = copula_cdf_indep;
    copula->h = copula_h_indep;
    copula->hinv = copula_hinv_indep;
    copula->aic = copula_aic_indep;
    copula->free = NULL; // Disabled.
    copula->data = NULL; // Disabled.

    return copula;
}
