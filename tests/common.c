// DML - Dependence Modeling Library
// Copyright (C) 2011-2012 Yasser González-Fernández <ygonzalezfernandez@gmail.com>

#include "config.h"

#include "src/dml.h"
#include "tests/common.h"

void
test_vectors_equal(const gsl_vector *x,
                   const gsl_vector *y,
                   double err_eps,
                   double rate_eps)
{
    size_t k, n;
    double xi, yi;
    double err, rate;

    k = 0;
    n = x->size;
    for (size_t i = 0; i < n; i++) {
        xi = gsl_vector_get(x, i);
        yi = gsl_vector_get(y, i);

        if (fabs(xi) < 0.01 || fabs(yi) < 0.01) {
            // Absolute error.
            err = fabs(xi - yi);
        } else {
            // Relative error.
            err = fabs((yi - xi) / xi);
        }
        k += (err <= err_eps);
    }
    rate = k / (double) n;
    g_assert(rate >= rate_eps);
}
