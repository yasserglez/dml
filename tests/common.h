// DML - Dependence Modeling Library
// Copyright (C) 2011-2012 Yasser González-Fernández <ygonzalezfernandez@gmail.com>

#ifndef TEST_H_
#define TEST_H_

#include <glib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>

void
test_vectors_equal(const gsl_vector *x,
                   const gsl_vector *y,
                   double err_eps,
                   double rate_eps);


void test_measure_alloc();
void test_measure_tau_small();
void test_measure_tau_large();
void test_measure_cvm_normal();
void test_measure_cvm_indep();


void test_indep_alloc();
void test_indep_fit();
void test_indep_pdf();
void test_indep_cdf();
void test_indep_h();
void test_indep_hinv();

void test_normal_alloc();
void test_normal_fit();
void test_normal_pdf();
void test_normal_cdf();
void test_normal_h();
void test_normal_hinv();
void test_normal_gof_normal();
void test_normal_gof_clayton();

void test_clayton_alloc();
void test_clayton_fit();
void test_clayton_pdf();
void test_clayton_cdf();
void test_clayton_h();
void test_clayton_hinv();

void test_rclayton90_alloc();
void test_rclayton90_fit();
void test_rclayton90_pdf();
void test_rclayton90_cdf();
void test_rclayton90_h();
void test_rclayton90_hinv();

void test_rclayton180_alloc();
void test_rclayton180_fit();
void test_rclayton180_pdf();
void test_rclayton180_cdf();
void test_rclayton180_h();
void test_rclayton180_hinv();

void test_rclayton270_alloc();
void test_rclayton270_fit();
void test_rclayton270_pdf();
void test_rclayton270_cdf();
void test_rclayton270_h();
void test_rclayton270_hinv();

void test_copula_select_indeptest_none();
void test_copula_select_indeptest_tau();
void test_copula_select_indeptest_cvm();
void test_copula_select_aic();


void test_cvine_alloc();
void test_cvine_ran_2d();
void test_cvine_fit_2d();
void test_cvine_ran_fit_tau_5d_normal_indep();
void test_cvine_ran_fit_cvm_5d_normal_indep();
void test_cvine_ran_fit_20d_normal_trunc();
void test_cvine_bugfix1();

void test_dvine_alloc();
void test_dvine_ran_2d();
void test_dvine_fit_2d();
void test_dvine_ran_fit_tau_5d_normal_indep();
void test_dvine_ran_fit_cvm_5d_normal_indep();
void test_dvine_ran_fit_20d_normal_trunc();

void test_rvine_alloc();
void test_rvine_ran_2d();
void test_rvine_fit_2d();
void test_rvine_ran_fit_3d_normal();
void test_rvine_ran_fit_tau_7d_normal();
void test_rvine_ran_fit_cvm_7d_normal();
void test_rvine_ran_fit_9d_normal_trunc();
void test_rvine_bugfix1();

#endif /* TEST_H_ */
