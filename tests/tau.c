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
#include "tests/common.h"

/* Test using the data set presented in the following reference:
 * Genest, C. and Favre, A. C. (2007) Everything you always want to known
 * about copula modeling but were afraid to ask. Journal of Hydrologic
 * Engineering, 12, 347-68.
 */
void
test_measure_tau_small()
{
    dml_measure_tau_t *tau;
    gsl_vector_view x_view, y_view;
    double x[] = { -2.224, -1.538, -0.807, 0.024, 0.052, 1.324 };
    double y[] = { 0.431, 1.035, 0.586, 1.465, 1.115, -0.847 };

    x_view = gsl_vector_view_array(x, 6);
    y_view = gsl_vector_view_array(y, 6);

    tau = dml_measure_tau_alloc(&x_view.vector, &y_view.vector);
    g_assert(fabs(dml_measure_tau_coef(tau) - 0.06) <= 0.01);
    g_assert(fabs(dml_measure_tau_pvalue(tau) - 0.85) <= 0.01);

    dml_measure_tau_free(tau);
}

void
test_measure_tau_large()
{
    FILE *f;
    char *path;
    gsl_matrix *data;
    gsl_vector_view x, y;
    dml_measure_tau_t *tau;

    data = gsl_matrix_alloc(1000, 2);
    path = g_build_filename("tests", "data", "tau_large.dat", NULL);
    f = fopen(path, "r");
    gsl_matrix_fscanf(f, data);

    x = gsl_matrix_column(data, 0);
    y = gsl_matrix_column(data, 1);
    tau = dml_measure_tau_alloc(&x.vector, &y.vector);
    g_assert(fabs(dml_measure_tau_coef(tau) - 0.03) <= 0.01);
    g_assert(fabs(dml_measure_tau_pvalue(tau) - 0.09) <= 0.01);

    dml_measure_tau_free(tau);
    gsl_matrix_free(data);
    fclose(f);
    g_free(path);
}
