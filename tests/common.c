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

#include "src/dml.h"
#include "tests/common.h"

void
test_vectors_equal(const gsl_vector *x,
                   const gsl_vector *y,
                   double err_eps,
                   double rate_eps)
{
    size_t k = 0, n = x->size;
    double xi, yi;
    double err, rate;

    for (size_t i = 0; i < n; i++) {
        xi = gsl_vector_get(x, i);
        yi = gsl_vector_get(y, i);
        if (fabs(xi) < 0.01 || fabs(yi) < 0.01) {
            err = fabs(xi - yi);
        } else {
            err = fabs((yi - xi) / xi);
        }
        k += (err <= err_eps);
    }
    rate = k / (double) n;

    g_assert(rate >= rate_eps);
}
