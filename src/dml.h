// DML - Dependence Modeling Library
// Copyright (C) 2011-2012 Yasser González-Fernández <ygonzalezfernandez@gmail.com>

#ifndef DML_H_
#define DML_H_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct dml_measure_tau_s {
    const gsl_vector *x; // Observations of the first variable.
    const gsl_vector *y; // Observations of the second variable.
    double coef; // Kendall's tau rank correlation coefficient.
    double pvalue; // p-value of the independence based on Kendall's tau.
} dml_measure_tau_t;


typedef enum dml_copula_selection_e {
    DML_COPULA_SELECTION_AIC, // Akaike Information Criterion.
} dml_copula_selection_t;

typedef enum dml_copula_indeptest_e {
    DML_COPULA_INDEPTEST_NONE, // Disabled.
    DML_COPULA_INDEPTEST_TAU, // Test based on Kendall's tau.
} dml_copula_indeptest_t;

typedef enum dml_copula_type_e {
    DML_COPULA_INDEP, // Independence copula.
    DML_COPULA_NORMAL, // Normal copula.
    DML_COPULA_CLAYTON, // Clayton copula.
    DML_COPULA_RCLAYTON90, // Clayton copula rotated 90 degrees.
    DML_COPULA_RCLAYTON180, // Clayton copula rotated 180 degrees.
    DML_COPULA_RCLAYTON270, // Clayton copula rotated 270 degrees.
} dml_copula_type_t;

typedef struct dml_copula_s {
    dml_copula_type_t type;
    void *data;

    void (*fit)(struct dml_copula_s *copula,
                const gsl_vector *u,
                const gsl_vector *v,
                dml_measure_tau_t *tau);
    void (*pdf)(const struct dml_copula_s *copula,
                const gsl_vector *u,
                const gsl_vector *v,
                gsl_vector *pdf);
    void (*cdf)(const struct dml_copula_s *copula,
                const gsl_vector *u,
                const gsl_vector *v,
                gsl_vector *cdf);
    void (*h)(const struct dml_copula_s *copula,
              const gsl_vector *u,
              const gsl_vector *v,
              gsl_vector *h);
    void (*hinv)(const struct dml_copula_s *copula,
                 const gsl_vector *u,
                 const gsl_vector *v,
                 gsl_vector *hinv);
    void (*aic)(const struct dml_copula_s *copula,
                const gsl_vector *u,
                const gsl_vector *v,
                double *aic);
    void (*free)(struct dml_copula_s *copula);
} dml_copula_t;


typedef enum dml_vine_weight_e {
    DML_VINE_WEIGHT_TAU, // Absolute value of Kendall's tau.
} dml_vine_weight_t;

typedef enum dml_vine_truncation_e {
    DML_VINE_TRUNCATION_NONE, // Disabled.
    DML_VINE_TRUNCATION_AIC, // Akaike Information Criterion.
} dml_vine_truncation_t;

typedef enum dml_vine_type_e {
    DML_VINE_CVINE, // Canonical vine.
    DML_VINE_DVINE, // D-vine.
    DML_VINE_RVINE, // Regular vine.
} dml_vine_type_t;

typedef struct dml_vine_s {
    dml_vine_type_t type;
    size_t dimension;
    size_t trees;
    size_t *order; // order[i] = j means that the variable i-th variable of the vine
                   // represents the j-th variable. Both i and j are zero-based.
    size_t **matrix; // R-vine matrix.
    dml_copula_t ***copulas; // Matrix with the parameters of the copulas.

    void (*fit)(struct dml_vine_s *vine,
                const gsl_matrix *data,
                const dml_vine_weight_t weight,
                const dml_vine_truncation_t truncation,
                const dml_copula_indeptest_t indeptest,
                const double indeptest_level,
                const dml_copula_type_t *types,
                const size_t types_size,
                const dml_copula_selection_t selection);
    void (*ran)(const struct dml_vine_s *vine,
                const gsl_rng *rng,
                gsl_matrix *data);
    void (*free)(struct dml_vine_s *vine);

} dml_vine_t;


dml_measure_tau_t *
dml_measure_tau_alloc(const gsl_vector *x, const gsl_vector *y);

double
dml_measure_tau_coef(dml_measure_tau_t *tau);

double
dml_measure_tau_pvalue(dml_measure_tau_t *tau);

void
dml_measure_tau_free(dml_measure_tau_t *tau);


dml_copula_t *
dml_copula_alloc(const dml_copula_type_t type);

dml_copula_t *
dml_copula_alloc_indep();

dml_copula_t *
dml_copula_alloc_normal(const double rho);

dml_copula_t *
dml_copula_alloc_clayton(const double theta);

dml_copula_t *
dml_copula_alloc_rclayton90(const double theta);

dml_copula_t *
dml_copula_alloc_rclayton180(const double theta);

dml_copula_t *
dml_copula_alloc_rclayton270(const double theta);

dml_copula_type_t
dml_copula_type(const dml_copula_t *copula);

dml_copula_t *
dml_copula_select(const gsl_vector *u,
                  const gsl_vector *v,
                  dml_measure_tau_t *tau,
                  const dml_copula_indeptest_t indeptest,
                  const double indeptest_level,
                  const dml_copula_type_t *types,
                  const size_t types_size,
                  const dml_copula_selection_t selection);

void
dml_copula_fit(dml_copula_t *copula,
               const gsl_vector *u,
               const gsl_vector *v,
               dml_measure_tau_t *tau);

void
dml_copula_ran(const dml_copula_t *copula,
               const gsl_rng *rng,
               gsl_vector *u,
               gsl_vector *v);

void
dml_copula_pdf(const dml_copula_t *copula,
               const gsl_vector *u,
               const gsl_vector *v,
               gsl_vector *pdf);

void
dml_copula_cdf(const dml_copula_t *copula,
               const gsl_vector *u,
               const gsl_vector *v,
               gsl_vector *cdf);

void
dml_copula_h(const dml_copula_t *copula,
             const gsl_vector *u,
             const gsl_vector *v,
             gsl_vector *h);

void
dml_copula_hinv(const dml_copula_t *copula,
                const gsl_vector *u,
                const gsl_vector *v,
                gsl_vector *hinv);

void
dml_copula_aic(const dml_copula_t *copula,
               const gsl_vector *u,
               const gsl_vector *v,
               double *aic);

void
dml_copula_free(dml_copula_t *copula);


dml_vine_t *
dml_vine_alloc(const dml_vine_type_t type, const size_t dimension);

dml_vine_t *
dml_vine_alloc_cvine(const size_t dimension);

dml_vine_t *
dml_vine_alloc_dvine(const size_t dimension);

dml_vine_t *
dml_vine_alloc_rvine(const size_t dimension);

dml_vine_type_t
dml_vine_type(const dml_vine_t *vine);

void
dml_vine_fit(dml_vine_t *vine,
             const gsl_matrix *data,
             const dml_vine_weight_t weight,
             const dml_vine_truncation_t truncation,
             const dml_copula_indeptest_t indeptest,
             const double indeptest_level,
             const dml_copula_type_t *types,
             const size_t types_size,
             const dml_copula_selection_t selection);

void
dml_vine_ran(const dml_vine_t *vine, const gsl_rng *rng, gsl_matrix *data);

void
dml_vine_free(dml_vine_t *vine);

#endif /* DML_H_ */
