// DML - Dependence Modeling Library
// Copyright (C) 2011-2012 Yasser González-Fernández <ygonzalezfernandez@gmail.com>

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
