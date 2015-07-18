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
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

#include "src/dml.h"

gsl_permutation *
dml_measure_x_rank(dml_measure_t *measure)
{
    gsl_permutation *x_perm;

    if (measure->x_rank == NULL) {
        x_perm = gsl_permutation_alloc(measure->x->size);
        measure->x_rank = gsl_permutation_alloc(measure->x->size);
        gsl_sort_vector_index(x_perm, measure->x);
        gsl_permutation_inverse(measure->x_rank, x_perm);
        gsl_permutation_free(x_perm);
    }

    return measure->x_rank;
}

gsl_permutation *
dml_measure_y_rank(dml_measure_t *measure)
{
    gsl_permutation *y_perm;

    if (measure->y_rank == NULL) {
        y_perm = gsl_permutation_alloc(measure->y->size);
        measure->y_rank = gsl_permutation_alloc(measure->y->size);
        gsl_sort_vector_index(y_perm, measure->y);
        gsl_permutation_inverse(measure->y_rank, y_perm);
        gsl_permutation_free(y_perm);
    }

    return measure->y_rank;
}
