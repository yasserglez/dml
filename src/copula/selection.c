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

#include <stdbool.h>

#include <glib.h>
#include <gsl/gsl_vector.h>

dml_copula_t *
dml_copula_selection(const gsl_vector *u,
                     const gsl_vector *v,
                     dml_measure_tau_t *tau,
                     const dml_copula_indeptest_t indeptest,
                     const double indeptest_level,
                     const dml_copula_type_t *types,
                     const size_t types_size,
                     const dml_copula_selection_t selection)
{
    dml_copula_t *selected, *candidate;
    double selected_fit = 0, candidate_fit = 0; // Initialized to avoid GCC warnings.
    bool tau_allocated;

    selected = NULL;

    if (tau == NULL) {
        tau = dml_measure_tau_alloc(u, v);
        tau_allocated = true;
    } else {
        tau_allocated = false;
    }

    // Independence tests.
    if (indeptest == DML_COPULA_INDEPTEST_TAU) {
        double pvalue = dml_measure_tau_pvalue(tau);
        if (pvalue >= indeptest_level) {
            // The null hypothesis is not rejected.
            selected = dml_copula_alloc(DML_COPULA_INDEP);
        }
    }

    // Goodness-of-fit.
    if (selected == NULL && types_size > 0) {
        if (types_size == 1) {
            selected = dml_copula_alloc(types[0]);
            dml_copula_fit(selected, u, v, tau);
        } else {
            for (size_t t = 0; t < types_size; t++) {
                if (dml_measure_tau_coef(tau) > 0
                        && (types[t] == DML_COPULA_RCLAYTON90
                                || types[t] == DML_COPULA_RCLAYTON270)) {
                    // Ignore the copula if it does not represent positive dependence.
                    continue;
                }
                if (dml_measure_tau_coef(tau) < 0
                        && (types[t] == DML_COPULA_CLAYTON
                                || types[t] == DML_COPULA_RCLAYTON180)) {
                    // Ignore the copula if it does not represent negative dependence.
                    continue;
                }

                candidate = dml_copula_alloc(types[t]);
                dml_copula_fit(candidate, u, v, tau);

                switch (selection) {
                case DML_COPULA_SELECTION_AIC:
                    dml_copula_aic(candidate, u, v, &candidate_fit);
                    break;
                default:
                    break;
                }
                if (selected == NULL) {
                    selected = candidate;
                    selected_fit = candidate_fit;
                } else if (candidate_fit < selected_fit) {
                    dml_copula_free(selected);
                    selected = candidate;
                    selected_fit = candidate_fit;
                } else {
                    dml_copula_free(candidate);
                }
            }
        }
    }

    if (tau_allocated) dml_measure_tau_free(tau);

    return selected;
}
