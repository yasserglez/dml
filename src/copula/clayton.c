// DML - Dependence Modeling Library
// Copyright (C) 2011-2012 Yasser González-Fernández <ygonzalezfernandez@gmail.com>

#include "config.h"

#include <glib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>

#include "src/dml.h"

static void
copula_fit_rclayton(dml_copula_t *copula,
                    const gsl_vector *u,
                    const gsl_vector *v,
                    dml_measure_t *measure)
{
    double *params;
    double theta;
    double coef;

    if (measure == NULL) {
        measure = dml_measure_alloc(u, v);
        coef = dml_measure_tau_coef(measure);
        dml_measure_free(measure);
    } else {
    	coef = dml_measure_tau_coef(measure);
    }
    if (dml_copula_type(copula) == DML_COPULA_CLAYTON
            || dml_copula_type(copula) == DML_COPULA_RCLAYTON180) {
        theta = CLAMP(2 * coef / (1 - coef), 0, 50);
    } else { // DML_COPULA_RCLAYTON90 or DML_COPULA_RCLAYTON270.
        theta = CLAMP(2 * coef / (1 + coef), -50, 0);
    }
    params = copula->data;
    params[0] = theta;
}

static void
copula_pdf_clayton(double theta,
                   const gsl_vector *u,
                   const gsl_vector *v,
                   gsl_vector *pdf)
{
    double eps;
    double theta_eps;
    double ui, vi, pdfi;

    theta_eps = 1e-4;
    eps = 1e-10;

    if (theta <= theta_eps) {
        gsl_vector_set_all(pdf, 1);
    } else {
        for (int i = 0; i < u->size; i++) {
            ui = gsl_vector_get(u, i);
            vi = gsl_vector_get(v, i);

            if (ui <= eps || vi <= eps) {
                pdfi = 0;
            } else {
                pdfi = (1 + theta) * pow(ui * vi, - 1 - theta)
                        * pow(pow(ui, -theta) + pow(vi, -theta) - 1, - 2 - 1/theta);
            }
            gsl_vector_set(pdf, i, pdfi);
        }
    }
}

static void
copula_pdf_rclayton(const dml_copula_t *copula,
                    const gsl_vector *u,
                    const gsl_vector *v,
                    gsl_vector *pdf)
{
	double theta;
	gsl_vector *uu, *vv;

	theta = ((double *) copula->data)[0];

	if (dml_copula_type(copula) == DML_COPULA_CLAYTON) {
		copula_pdf_clayton(theta, u, v, pdf);
	} else if (dml_copula_type(copula) == DML_COPULA_RCLAYTON90) {
		vv = gsl_vector_alloc(v->size);
		gsl_vector_set_all(vv, 1);
		gsl_vector_sub(vv, v);
		copula_pdf_clayton(-theta, u, vv, pdf);
		gsl_vector_free(vv);
	} else if (dml_copula_type(copula) == DML_COPULA_RCLAYTON180) {
		uu = gsl_vector_alloc(u->size);
		gsl_vector_set_all(uu, 1);
		gsl_vector_sub(uu, u);
		vv = gsl_vector_alloc(v->size);
		gsl_vector_set_all(vv, 1);
		gsl_vector_sub(vv, v);
		copula_pdf_clayton(theta, uu, vv, pdf);
		gsl_vector_free(uu);
		gsl_vector_free(vv);
	} else { // DML_COPULA_RCLAYTON270.
		uu = gsl_vector_alloc(u->size);
		gsl_vector_set_all(uu, 1);
		gsl_vector_sub(uu, u);
		copula_pdf_clayton(-theta, uu, v, pdf);
		gsl_vector_free(uu);
	}
}

static void
copula_cdf_clayton(double theta,
                   const gsl_vector *u,
                   const gsl_vector *v,
                   gsl_vector *cdf)
{
	double eps;
	double theta_eps;
    double ui, vi, cdfi;

    theta_eps = 1e-4;
    eps = 1e-10;

    if (theta <= theta_eps) {
        gsl_vector_memcpy(cdf, u);
        gsl_vector_mul(cdf, v);
    } else {
    	for (int i = 0; i < u->size; i++) {
    		ui = gsl_vector_get(u, i);
    		vi = gsl_vector_get(v, i);

    		if (ui >= 1 - eps && vi >= 1 - eps) {
    			cdfi = 1;
    		} else if (ui >= 1 - eps) {
    			cdfi = vi;
    		} else if (vi >= 1 - eps) {
    			cdfi = ui;
    		} else if (ui <= eps || vi <= eps) {
    			cdfi = 0;
    		} else {
                ui = CLAMP(ui, eps, 1 - eps);
                vi = CLAMP(vi, eps, 1 - eps);
    			cdfi = pow(pow(ui, -theta) + pow(vi, -theta) - 1, -1 / theta);
    		}
    		gsl_vector_set(cdf, i, cdfi);
    	}
    }
}

static void
copula_cdf_rclayton(const dml_copula_t *copula,
                    const gsl_vector *u,
                    const gsl_vector *v,
                    gsl_vector *cdf)
{
    double theta;
    gsl_vector *tmp, *uu, *vv;

    theta = ((double *) copula->data)[0];

    if (dml_copula_type(copula) == DML_COPULA_CLAYTON) {
        copula_cdf_clayton(theta, u, v, cdf);
    } else if (dml_copula_type(copula) == DML_COPULA_RCLAYTON90) {
        vv = gsl_vector_alloc(v->size);
        gsl_vector_set_all(vv, 1);
        gsl_vector_sub(vv, v);
        tmp = gsl_vector_alloc(cdf->size);
        copula_cdf_clayton(-theta, u, vv, tmp);
        gsl_vector_memcpy(cdf, u);
        gsl_vector_sub(cdf, tmp);
        gsl_vector_free(tmp);
        gsl_vector_free(vv);
    } else if (dml_copula_type(copula) == DML_COPULA_RCLAYTON180) {
        uu = gsl_vector_alloc(u->size);
        gsl_vector_set_all(uu, 1);
        gsl_vector_sub(uu, u);
        vv = gsl_vector_alloc(v->size);
        gsl_vector_set_all(vv, 1);
        gsl_vector_sub(vv, v);
        tmp = gsl_vector_alloc(cdf->size);
        copula_cdf_clayton(theta, uu, vv, tmp);
        gsl_vector_set_all(cdf, -1);
        gsl_vector_add(cdf, u);
        gsl_vector_add(cdf, v);
        gsl_vector_add(cdf, tmp);
        gsl_vector_free(tmp);
        gsl_vector_free(uu);
        gsl_vector_free(vv);
    } else { // DML_COPULA_RCLAYTON270.
        uu = gsl_vector_alloc(u->size);
        gsl_vector_set_all(uu, 1);
        gsl_vector_sub(uu, u);
        tmp = gsl_vector_alloc(cdf->size);
        copula_cdf_clayton(-theta, uu, v, tmp);
        gsl_vector_memcpy(cdf, v);
        gsl_vector_sub(cdf, tmp);
        gsl_vector_free(tmp);
        gsl_vector_free(uu);
    }
}

static void
copula_h_clayton(double theta,
                 const gsl_vector *u,
                 const gsl_vector *v,
                 gsl_vector *h)
{
    double eps;
    double theta_eps;
    double ui, vi, hi;

    theta_eps = 1e-4;
    eps = 1e-10;

    if (theta <= theta_eps) {
        gsl_vector_memcpy(h, u);
    } else {
    	for (int i = 0; i < u->size; i++) {
    		ui = gsl_vector_get(u, i);
    		vi = gsl_vector_get(v, i);

            if (ui <= eps) {
                hi = eps;
            } else if (ui >= 1 - eps || (vi <= eps && theta >= 1)) {
                hi = 1 - eps;
            } else {
                ui = CLAMP(ui, eps, 1 - eps);
                vi = CLAMP(vi, eps, 1 - eps);
                hi = pow(vi, -theta - 1)
						* pow(pow(ui, -theta) + pow(vi, -theta) - 1,
								-1 - 1 / theta);
                hi = CLAMP(hi, eps, 1 - eps);
            }
    		gsl_vector_set(h, i, hi);
    	}
    }
}

static void
copula_h_rclayton(const dml_copula_t *copula,
                  const gsl_vector *u,
                  const gsl_vector *v,
                  gsl_vector *h)
{
    double theta;
    gsl_vector *tmp, *uu, *vv;

    theta = ((double *) copula->data)[0];

    if (dml_copula_type(copula) == DML_COPULA_CLAYTON) {
        copula_h_clayton(theta, u, v, h);
    } else if (dml_copula_type(copula) == DML_COPULA_RCLAYTON90) {
        vv = gsl_vector_alloc(v->size);
        gsl_vector_set_all(vv, 1);
        gsl_vector_sub(vv, v);
        copula_h_clayton(-theta, u, vv, h);
        gsl_vector_free(vv);
    } else if (dml_copula_type(copula) == DML_COPULA_RCLAYTON180) {
        uu = gsl_vector_alloc(u->size);
        gsl_vector_set_all(uu, 1);
        gsl_vector_sub(uu, u);
        vv = gsl_vector_alloc(v->size);
        gsl_vector_set_all(vv, 1);
        gsl_vector_sub(vv, v);
        tmp = gsl_vector_alloc(h->size);
        copula_h_clayton(theta, uu, vv, tmp);
        gsl_vector_set_all(h, 1);
        gsl_vector_sub(h, tmp);
        gsl_vector_free(tmp);
        gsl_vector_free(uu);
        gsl_vector_free(vv);
    } else { // DML_COPULA_RCLAYTON270.
        uu = gsl_vector_alloc(u->size);
        gsl_vector_set_all(uu, 1);
        gsl_vector_sub(uu, u);
        tmp = gsl_vector_alloc(h->size);
        copula_h_clayton(-theta, uu, v, tmp);
        gsl_vector_set_all(h, 1);
        gsl_vector_sub(h, tmp);
        gsl_vector_free(tmp);
        gsl_vector_free(uu);
    }
}

static void
copula_hinv_clayton(double theta,
                    const gsl_vector *u,
                    const gsl_vector *v,
                    gsl_vector *hinv)
{
    double eps;
    double theta_eps;
    double ui, vi, hinvi;

    theta_eps = 1e-4;
    eps = 1e-10;

    if (theta <= theta_eps) {
        gsl_vector_memcpy(hinv, u);
    } else {
		for (int i = 0; i < u->size; i++) {
			ui = gsl_vector_get(u, i);
			vi = gsl_vector_get(v, i);

	        if (ui <= eps) {
	            hinvi = eps;
	        } else if (ui >= 1 - eps) {
	            hinvi = 1 - eps;
	        } else if (vi <= eps) {
	            hinvi = eps;
	        } else {
	            ui = CLAMP(ui, eps, 1 - eps);
	            vi = CLAMP(vi, eps, 1 - eps);
                hinvi = pow(pow(ui * pow(vi, theta + 1), -theta / (theta + 1))
                            + 1 - pow(vi, -theta), -1 / theta);
                hinvi = CLAMP(hinvi, eps, 1 - eps);
	        }
			gsl_vector_set(hinv, i, hinvi);
		}
    }
}

static void
copula_hinv_rclayton(const dml_copula_t *copula,
                     const gsl_vector *u,
                     const gsl_vector *v,
                     gsl_vector *hinv)
{
    double theta;
    gsl_vector *tmp, *uu, *vv;

    theta = ((double *) copula->data)[0];

    if (dml_copula_type(copula) == DML_COPULA_CLAYTON) {
        copula_hinv_clayton(theta, u, v, hinv);
    } else if (dml_copula_type(copula) == DML_COPULA_RCLAYTON90) {
        vv = gsl_vector_alloc(v->size);
        gsl_vector_set_all(vv, 1);
        gsl_vector_sub(vv, v);
        copula_hinv_clayton(-theta, u, vv, hinv);
        gsl_vector_free(vv);
    } else if (dml_copula_type(copula) == DML_COPULA_RCLAYTON180) {
        uu = gsl_vector_alloc(u->size);
        gsl_vector_set_all(uu, 1);
        gsl_vector_sub(uu, u);
        vv = gsl_vector_alloc(v->size);
        gsl_vector_set_all(vv, 1);
        gsl_vector_sub(vv, v);
        tmp = gsl_vector_alloc(hinv->size);
        copula_hinv_clayton(theta, uu, vv, tmp);
        gsl_vector_set_all(hinv, 1);
        gsl_vector_sub(hinv, tmp);
        gsl_vector_free(tmp);
        gsl_vector_free(uu);
        gsl_vector_free(vv);
    } else { // DML_COPULA_RCLAYTON270.
        uu = gsl_vector_alloc(v->size);
        gsl_vector_set_all(uu, 1);
        gsl_vector_sub(uu, u);
        tmp = gsl_vector_alloc(hinv->size);
        copula_hinv_clayton(-theta, uu, v, tmp);
        gsl_vector_set_all(hinv, 1);
        gsl_vector_sub(hinv, tmp);
        gsl_vector_free(tmp);
        gsl_vector_free(uu);
    }
}

static void
copula_aic_rclayton(const dml_copula_t *copula,
                    const gsl_vector *u,
                    const gsl_vector *v,
                    double *aic)
{
    gsl_vector *pdf;
    double loglik;

    pdf = gsl_vector_alloc(u->size);

    loglik = 0;
    dml_copula_pdf(copula, u, v, pdf);
    for (size_t i = 0; i < pdf->size; i++) {
        loglik += log(gsl_vector_get(pdf, i));
    }
    *aic = -2*loglik + 2;

    gsl_vector_free(pdf);
}


static void
copula_free_rclayton(dml_copula_t *copula)
{
    g_free(copula->data);
}

static dml_copula_t *
copula_alloc_rclayton(int degrees, double theta)
{
    dml_copula_t *copula;
    double *params;

    copula = g_malloc(sizeof(dml_copula_t));
    copula->fit = copula_fit_rclayton;
    copula->pdf = copula_pdf_rclayton;
    copula->cdf = copula_cdf_rclayton;
    copula->h = copula_h_rclayton;
    copula->hinv = copula_hinv_rclayton;
    copula->aic = copula_aic_rclayton;
    copula->gof = NULL;
    copula->free = copula_free_rclayton;
    copula->data = g_malloc(sizeof(double));
    params = copula->data;
    switch (degrees) {
    case 0:
    	copula->type = DML_COPULA_CLAYTON;
    	params[0] = CLAMP(theta, 0, 50);
    	break;
    case 90:
    	copula->type = DML_COPULA_RCLAYTON90;
    	params[0] = CLAMP(theta, -50, 0);
    	break;
    case 180:
    	copula->type = DML_COPULA_RCLAYTON180;
    	params[0] = CLAMP(theta, 0, 50);
    	break;
    case 270:
    	copula->type = DML_COPULA_RCLAYTON270;
    	params[0] = CLAMP(theta, -50, 0);
    	break;
    default:
    	break;
    }

    return copula;
}

dml_copula_t *
dml_copula_alloc_clayton(const double theta)
{
	return copula_alloc_rclayton(0, theta);
}

dml_copula_t *
dml_copula_alloc_rclayton90(const double theta)
{
	return copula_alloc_rclayton(90, theta);
}

dml_copula_t *
dml_copula_alloc_rclayton180(const double theta)
{
	return copula_alloc_rclayton(180, theta);
}

dml_copula_t *
dml_copula_alloc_rclayton270(const double theta)
{
    return copula_alloc_rclayton(270, theta);
}
