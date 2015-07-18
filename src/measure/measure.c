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
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

#include "src/dml.h"

dml_measure_t *
dml_measure_alloc(const gsl_vector *x, const gsl_vector *y)
{
    dml_measure_t *measure;

    measure = g_malloc(sizeof(dml_measure_t));
    measure->x = x;
    measure->y = y;
    measure->x_rank = NULL;
    measure->y_rank = NULL;
    measure->tau_coef = GSL_NAN;
    measure->tau_pvalue = GSL_NAN;
    measure->cvm_stat = GSL_NAN;
    measure->cvm_pvalue = GSL_NAN;

    return measure;
}

inline void
dml_measure_free(dml_measure_t *measure)
{
    if (measure->x_rank != NULL) {
        gsl_permutation_free(measure->x_rank);
    }
    if (measure->y_rank != NULL) {
        gsl_permutation_free(measure->y_rank);
    }
    g_free(measure);
}
