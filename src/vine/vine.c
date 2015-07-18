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
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>

#include "src/dml.h"

inline dml_vine_t *
dml_vine_alloc(const dml_vine_type_t type, const size_t dim)
{
    dml_vine_t *vine;

    switch (type) {
    case DML_VINE_CVINE:
        vine = dml_vine_alloc_cvine(dim);
        break;
    case DML_VINE_DVINE:
        vine = dml_vine_alloc_dvine(dim);
        break;
    case DML_VINE_RVINE:
        vine = dml_vine_alloc_rvine(dim);
        break;
    default:
        vine = NULL;
        break;
    }

    return vine;
}

inline dml_vine_type_t
dml_vine_type(const dml_vine_t *vine)
{
    return vine->type;
}

inline void
dml_vine_fit(dml_vine_t *vine,
             const gsl_matrix *data,
             const dml_vine_weight_t weight,
             const dml_vine_trunc_t trunc,
             const dml_copula_indeptest_t indeptest,
             const double indeptest_level,
             const dml_copula_type_t *types,
             const size_t types_size,
             const dml_copula_select_t select,
             const gsl_rng *rng)
{
    vine->fit(vine, data, weight, trunc, indeptest, indeptest_level, types,
              types_size, select, rng);
}

inline void
dml_vine_ran(const dml_vine_t *vine, const gsl_rng *rng, gsl_matrix *data)
{
    if (vine->trees > 0) {
        vine->ran(vine, rng, data);
    } else {
        // Vine without trees. Independence.
        for (size_t i = 0; i < data->size1; i++) {
            for (size_t j = 0; j < data->size2; j++) {
                gsl_matrix_set(data, i, j, gsl_rng_uniform(rng));
            }
        }
    }
}

inline void
dml_vine_free(dml_vine_t *vine)
{
    if (vine->free != NULL) {
        vine->free(vine);
    }
    g_free(vine);
}
