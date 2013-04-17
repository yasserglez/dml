/*
 * DML - Dependence Modeling Library
 * Copyright (C) 2011-2013 Yasser Gonzalez-Fernandez <ygonzalezfernandez@gmail.com>
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

#include "tests/common.h"

int
main(int argc, char **argv)
{
    g_test_init(&argc, &argv, NULL);

    g_test_add_func("/measure/alloc", test_measure_alloc);
    g_test_add_func("/measure/tau/small", test_measure_tau_small);
    g_test_add_func("/measure/tau/large", test_measure_tau_large);
    g_test_add_func("/measure/cvm/normal", test_measure_cvm_normal);
    g_test_add_func("/measure/cvm/indep", test_measure_cvm_indep);

    g_test_add_func("/copula/indep/alloc", test_indep_alloc);
    g_test_add_func("/copula/indep/fit", test_indep_fit);
    g_test_add_func("/copula/indep/pdf", test_indep_pdf);
    g_test_add_func("/copula/indep/cdf", test_indep_cdf);
    g_test_add_func("/copula/indep/h", test_indep_h);
    g_test_add_func("/copula/indep/hinv", test_indep_hinv);
    g_test_add_func("/copula/normal/alloc", test_normal_alloc);
    g_test_add_func("/copula/normal/fit", test_normal_fit);
    g_test_add_func("/copula/normal/pdf", test_normal_pdf);
    g_test_add_func("/copula/normal/cdf", test_normal_cdf);
    g_test_add_func("/copula/normal/h", test_normal_h);
    g_test_add_func("/copula/normal/hinv", test_normal_hinv);
    g_test_add_func("/copula/clayton/alloc", test_clayton_alloc);
    g_test_add_func("/copula/clayton/fit", test_clayton_fit);
    g_test_add_func("/copula/clayton/pdf", test_clayton_pdf);
    g_test_add_func("/copula/clayton/cdf", test_clayton_cdf);
    g_test_add_func("/copula/clayton/h", test_clayton_h);
    g_test_add_func("/copula/clayton/hinv", test_clayton_hinv);
    g_test_add_func("/copula/rclayton90/alloc", test_rclayton90_alloc);
    g_test_add_func("/copula/rclayton90/fit", test_rclayton90_fit);
    g_test_add_func("/copula/rclayton90/pdf", test_rclayton90_pdf);
    g_test_add_func("/copula/rclayton90/cdf", test_rclayton90_cdf);
    g_test_add_func("/copula/rclayton90/h", test_rclayton90_h);
    g_test_add_func("/copula/rclayton90/hinv", test_rclayton90_hinv);
    g_test_add_func("/copula/rclayton180/alloc", test_rclayton180_alloc);
    g_test_add_func("/copula/rclayton180/fit", test_rclayton180_fit);
    g_test_add_func("/copula/rclayton180/pdf", test_rclayton180_pdf);
    g_test_add_func("/copula/rclayton180/cdf", test_rclayton180_cdf);
    g_test_add_func("/copula/rclayton180/h", test_rclayton180_h);
    g_test_add_func("/copula/rclayton180/hinv", test_rclayton180_hinv);
    g_test_add_func("/copula/rclayton270/alloc", test_rclayton270_alloc);
    g_test_add_func("/copula/rclayton270/fit", test_rclayton270_fit);
    g_test_add_func("/copula/rclayton270/pdf", test_rclayton270_pdf);
    g_test_add_func("/copula/rclayton270/cdf", test_rclayton270_cdf);
    g_test_add_func("/copula/rclayton270/h", test_rclayton270_h);
    g_test_add_func("/copula/rclayton270/hinv", test_rclayton270_hinv);
    g_test_add_func("/copula/select/indeptest_none", test_copula_select_indeptest_none);
    g_test_add_func("/copula/select/indeptest_tau", test_copula_select_indeptest_tau);
    g_test_add_func("/copula/select/indeptest_cvm", test_copula_select_indeptest_cvm);
    g_test_add_func("/copula/select/aic", test_copula_select_aic);

    g_test_add_func("/vine/cvine/alloc", test_cvine_alloc);
    g_test_add_func("/vine/cvine/ran_2d", test_cvine_ran_2d);
    g_test_add_func("/vine/cvine/fit_2d", test_cvine_fit_2d);
    g_test_add_func("/vine/cvine/ran_fit_tau_5d_normal_indep", test_cvine_ran_fit_tau_5d_normal_indep);
    g_test_add_func("/vine/cvine/ran_fit_cvm_5d_normal_indep", test_cvine_ran_fit_cvm_5d_normal_indep);
    g_test_add_func("/vine/cvine/ran_fit_20d_normal_trunc", test_cvine_ran_fit_20d_normal_trunc);
    g_test_add_func("/vine/cvine/bugfix1", test_cvine_bugfix1);
    g_test_add_func("/vine/dvine/alloc", test_dvine_alloc);
    g_test_add_func("/vine/dvine/ran_2d", test_dvine_ran_2d);
    g_test_add_func("/vine/dvine/fit_2d", test_dvine_fit_2d);
    g_test_add_func("/vine/dvine/ran_fit_tau_5d_normal_indep", test_dvine_ran_fit_tau_5d_normal_indep);
    g_test_add_func("/vine/dvine/ran_fit_cvm_5d_normal_indep", test_dvine_ran_fit_cvm_5d_normal_indep);//    g_test_add_func("/vine/dvine/ran_fit_20d_normal_trunc", test_dvine_ran_fit_20d_normal_trunc);
    g_test_add_func("/vine/rvine/alloc", test_rvine_alloc);
    g_test_add_func("/vine/rvine/ran_2d", test_rvine_ran_2d);
    g_test_add_func("/vine/rvine/fit_2d", test_rvine_fit_2d);
    g_test_add_func("/vine/rvine/ran_fit_3d_normal", test_rvine_ran_fit_3d_normal);
    g_test_add_func("/vine/rvine/ran_fit_tau_7d_normal", test_rvine_ran_fit_tau_7d_normal);
    g_test_add_func("/vine/rvine/ran_fit_cvm_7d_normal", test_rvine_ran_fit_cvm_7d_normal);
    g_test_add_func("/vine/rvine/ran_fit_9d_normal_trunc", test_rvine_ran_fit_9d_normal_trunc);
    g_test_add_func("/vine/rvine/bugfix1", test_rvine_bugfix1);
    g_test_add_func("/vine/rvine/bugfix2", test_rvine_bugfix2);
    g_test_add_func("/vine/rvine/bugfix3", test_rvine_bugfix3);

    return g_test_run();
}
