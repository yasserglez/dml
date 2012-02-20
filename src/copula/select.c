// DML - Dependence Modeling Library
// Copyright (C) 2011-2012 Yasser González-Fernández <ygonzalezfernandez@gmail.com>

#include "config.h"

#include <stdbool.h>

#include <glib.h>
#include <gsl/gsl_vector.h>

#include "src/dml.h"

dml_copula_t *
dml_copula_select(const gsl_vector *u,
                  const gsl_vector *v,
                  dml_measure_t *measure,
                  const dml_copula_indeptest_t indeptest,
                  const double indeptest_level,
                  const dml_copula_type_t *types,
                  const size_t types_size,
                  const dml_copula_select_t select,
                  const double gof_level,
                  const gsl_rng *rng)
{
    dml_copula_t *selected, *candidate;
    double selected_fit = 0, candidate_fit = 0; // Initialized to avoid GCC warnings.
    double indeptest_pvalue = -1;
    bool measure_dealloc;

    selected = NULL;

    if (measure == NULL) {
        measure = dml_measure_alloc(u, v);
        measure_dealloc = true;
    } else {
        measure_dealloc = false;
    }

    // Independence tests.
    if (indeptest == DML_COPULA_INDEPTEST_TAU) {
        indeptest_pvalue = dml_measure_tau_pvalue(measure);
    } else if (indeptest == DML_COPULA_INDEPTEST_CVM) {
        indeptest_pvalue = dml_measure_cvm_pvalue(measure, rng);
    }
    if (indeptest_pvalue >= indeptest_level) {
        // The null hypothesis was not rejected.
        selected = dml_copula_alloc(DML_COPULA_INDEP);
    }

    // Goodness-of-fit.
    if (selected == NULL && types_size > 0) {
        if (types_size == 1) {
            selected = dml_copula_alloc(types[0]);
            dml_copula_fit(selected, u, v, measure);

            if (select == DML_COPULA_SELECT_GOF) {
                dml_copula_gof(selected, u, v, measure, rng, &selected_fit);
                if (selected_fit < gof_level) {
                    dml_copula_free(selected);
                    selected = NULL;
                }
            }
        } else {
            for (size_t t = 0; t < types_size; t++) {
                if (dml_measure_tau_coef(measure) > 0
                        && (types[t] == DML_COPULA_RCLAYTON90
                                || types[t] == DML_COPULA_RCLAYTON270)) {
                    // Ignore the copula if it does not represent positive dependence.
                    continue;
                }
                if (dml_measure_tau_coef(measure) < 0
                        && (types[t] == DML_COPULA_CLAYTON
                                || types[t] == DML_COPULA_RCLAYTON180)) {
                    // Ignore the copula if it does not represent negative dependence.
                    continue;
                }

                candidate = dml_copula_alloc(types[t]);
                dml_copula_fit(candidate, u, v, measure);

                switch (select) {
                case DML_COPULA_SELECT_AIC:
                    dml_copula_aic(candidate, u, v, &candidate_fit);
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
                    break;
                case DML_COPULA_SELECT_GOF:
                    dml_copula_gof(candidate, u, v, measure, rng, &candidate_fit);
                    if (candidate_fit >= gof_level) {
                        // The candidate copula was not rejected.
                        if (selected == NULL) {
                            selected = candidate;
                            selected_fit = candidate_fit;
                        } else if (candidate_fit > selected_fit) {
                            dml_copula_free(selected);
                            selected = candidate;
                            selected_fit = candidate_fit;
                        } else {
                            dml_copula_free(candidate);
                        }
                    } else {
                        // The candidate copula was rejected.
                        dml_copula_free(candidate);
                    }
                    break;
                default:
                    break;
                }
            }
        }
    }

    // Non-parametric estimation.
    if (selected == NULL) {
        // TODO: Kernel-smoothed empirical copula.
    }

    if (measure_dealloc) dml_measure_free(measure);

    return selected;
}
