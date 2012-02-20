// DML - Dependence Modeling Library
// Copyright (C) 2011-2012 Yasser González-Fernández <ygonzalezfernandez@gmail.com>

#include "config.h"

#include <glib.h>
#include <gsl/gsl_math.h>

#include "src/dml.h"

dml_measure_t *
dml_measure_alloc(const gsl_vector *x, const gsl_vector *y)
{
    dml_measure_t *measure;

    measure = g_malloc(sizeof(dml_measure_t));
    measure->x = x;
    measure->y = y;
    measure->tau_coef = GSL_NAN;
    measure->tau_pvalue = GSL_NAN;
    measure->cvm_stat = GSL_NAN;
    measure->cvm_pvalue = GSL_NAN;

    return measure;
}

inline void
dml_measure_free(dml_measure_t *measure)
{
    g_free(measure);
}
