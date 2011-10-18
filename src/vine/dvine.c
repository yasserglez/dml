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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>

static void
vine_fit_dvine(dml_vine_t *vine,
               const gsl_matrix *data,
               dml_vine_weight_t weight,
               dml_vine_truncation_t truncation,
               dml_copula_indeptest_t indeptest,
               double indeptest_level,
               const dml_copula_type_t *types,
               size_t types_size,
               dml_copula_selection_t selection)
{
}

static void
vine_ran_dvine(const dml_vine_t *vine,
               const gsl_rng *rng,
               gsl_matrix *data)
{
}

static void
vine_free_dvine(dml_vine_t *vine)
{
}

dml_vine_t *
dml_vine_alloc_dvine(const size_t dimension)
{
    dml_vine_t *vine;

    vine = g_malloc(sizeof(dml_vine_t));
    vine->type = DML_VINE_DVINE;
    vine->fit = vine_fit_dvine;
    vine->ran = vine_ran_dvine;
    vine->free = vine_free_dvine;

    return vine;
}
