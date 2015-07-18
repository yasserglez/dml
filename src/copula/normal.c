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

#include <stdbool.h>

#include <glib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>

#include "src/dml.h"

// Evaluation of the bivariate normal distribution functions. Portions
// of code from http://www.math.wsu.edu/faculty/genz/software/mvtdstpack.f.
// The FORTRAN code was translated into C using f2c (version 20090411).

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define max(a,b) ((a) >= (b) ? (a) : (b))

static double
mvphi(double *z)
{
    static double a[44] = { .610143081923200417926465815756,
        -.434841272712577471828182820888, .176351193643605501125840298123,
        -.060710795609249414860051215825, .017712068995694114486147141191,
        -.004321119385567293818599864968, 8.54216676887098678819832055e-4,
        -1.2715509060916274262889394e-4, 1.1248167243671189468847072e-5,
        3.13063885421820972630152e-7, -2.70988068537762022009086e-7,
        3.0737622701407688440959e-8, 2.515620384817622937314e-9,
        -1.02892992132031912759e-9, 2.9944052119949939363e-11,
        2.605178968726693629e-11, -2.634839924171969386e-12,
        -6.43404509890636443e-13, 1.12457401801663447e-13,
        1.7281533389986098e-14, -4.264101694942375e-15, -5.45371977880191e-16,
        1.58697607761671e-16, 2.0899837844334e-17, -5.900526869409e-18,
        -9.41893387554e-19, 2.1497735647e-19, 4.6660985008e-20,
        -7.243011862e-21, -2.387966824e-21, 1.91177535e-22, 1.20482568e-22,
        -6.72377e-25, -5.747997e-24, -4.28493e-25, 2.44856e-25, 4.3793e-26,
        -8.151e-27, -3.089e-27, 9.3e-29, 1.74e-28, 1.6e-29, -8e-30, -2e-30 };

    static double b;
    static int i;
    static double p, t, bm, bp, xa;

    xa = abs(*z) / 1.414213562373095048801688724209;
    if (xa > 100.) {
        p = 0.;
    } else {
        t = (xa * 8 - 30) / (xa * 4 + 15);
        bm = 0.;
        b = 0.;
        for (i = 24; i >= 0; --i) {
            bp = b;
            b = bm;
            bm = t * b - bp + a[i];
        }
        p = exp(-xa * xa) * (bm - bp) / 4;
    }
    if (*z > 0.) {
        p = 1 - p;
    }
    return p;
}

static double
mvbvu(double *sh, double *sk, double *r)
{
    static struct {
        double e_1[3];
        double fill_2[7];
        double e_3[6];
        double fill_4[4];
        double e_5[10];
    } equiv_21 = { { .1713244923791705, .3607615730481384, .4679139345726904 },
        { 0 }, { .04717533638651177, .1069393259953183, .1600783285433464,
            .2031674267230659, .2334925365383547, .2491470458134029 }, { 0 }, {
            .01761400713915212, .04060142980038694, .06267204833410906,
            .08327674157670475, .1019301198172404, .1181945319615184,
            .1316886384491766, .1420961093183821, .1491729864726037,
            .1527533871307259 } };

#define w ((double *)&equiv_21)

    static struct {
        double e_1[3];
        double fill_2[7];
        double e_3[6];
        double fill_4[4];
        double e_5[10];
    } equiv_22 = {
        { -.9324695142031522, -.6612093864662647, -.238619186083197 }, { 0 }, {
            -.9815606342467191, -.904117256370475, -.769902674194305,
            -.5873179542866171, -.3678314989981802, -.1252334085114692 }, { 0 },
        { -.9931285991850949, -.9639719272779138, -.9122344282513259,
            -.8391169718222188, -.7463319064601508, -.636053680726515,
            -.5108670019508271, -.3737060887154196, -.2277858511416451,
            -.07652652113349733 } };

#define x ((double *)&equiv_22)

    int i__1;
    double d__1, d__2, d__3, d__4;
    static double a, b, c, d, h;
    static int i__;
    static double k;
    static int lg;
    static double as;
    static int ng;
    static double bs, hk, hs, sn, rs, xs, bvn, asr;

    if (abs(*r) < .3f) {
        ng = 1;
        lg = 3;
    } else if (abs(*r) < .75f) {
        ng = 2;
        lg = 6;
    } else {
        ng = 3;
        lg = 10;
    }
    h = *sh;
    k = *sk;
    hk = h * k;
    bvn = 0.;
    if (abs(*r) < .925f) {
        hs = (h * h + k * k) / 2;
        asr = asin(*r);
        i__1 = lg;
        for (i__ = 1; i__ <= i__1; ++i__) {
            sn = sin(asr * (x[i__ + ng * 10 - 11] + 1) / 2);
            bvn += w[i__ + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * sn));
            sn = sin(asr * (-x[i__ + ng * 10 - 11] + 1) / 2);
            bvn += w[i__ + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * sn));
        }
        d__1 = -h;
        d__2 = -k;
        bvn = bvn * asr / 12.566370614359172 + mvphi(&d__1) * mvphi(&d__2);
    } else {
        if (*r < 0.) {
            k = -k;
            hk = -hk;
        }
        if (abs(*r) < 1.) {
            as = (1 - *r) * (*r + 1);
            a = sqrt(as);
            d__1 = h - k;
            bs = d__1 * d__1;
            c = (4 - hk) / 8;
            d = (12 - hk) / 16;
            bvn = a * exp(-(bs / as + hk) / 2)
                    * (1 - c * (bs - as) * (1 - d * bs / 5) / 3
                            + c * d * as * as / 5);
            if (hk > -160.) {
                b = sqrt(bs);
                d__1 = -b / a;
                bvn -= exp(-hk / 2) * sqrt(6.283185307179586) * mvphi(&d__1) * b
                        * (1 - c * bs * (1 - d * bs / 5) / 3);
            }
            a /= 2;
            i__1 = lg;
            for (i__ = 1; i__ <= i__1; ++i__) {
                d__1 = a * (x[i__ + ng * 10 - 11] + 1);
                xs = d__1 * d__1;
                rs = sqrt(1 - xs);
                bvn += a * w[i__ + ng * 10 - 11]
                        * (exp(-bs / (xs * 2) - hk / (rs + 1)) / rs
                                - exp(-(bs / xs + hk) / 2)
                                        * (c * xs * (d * xs + 1) + 1));
                d__1 = -x[i__ + ng * 10 - 11] + 1;
                xs = as * (d__1 * d__1) / 4;
                rs = sqrt(1 - xs);
                bvn += a * w[i__ + ng * 10 - 11] * exp(-(bs / xs + hk) / 2)
                        * (exp(-hk * (1 - rs) / ((rs + 1) * 2)) / rs
                                - (c * xs * (d * xs + 1) + 1));
            }
            bvn = -bvn / 6.283185307179586;
        }
        if (*r > 0.) {
            d__1 = -max(h, k);
            bvn += mvphi(&d__1);
        }
        if (*r < 0.) {
            d__3 = -h;
            d__4 = -k;
            d__1 = 0., d__2 = mvphi(&d__3) - mvphi(&d__4);
            bvn = -bvn + max(d__1, d__2);
        }
    }
    return bvn;
}

#undef x
#undef w

static double
mvbvn(double *lower, double *upper, int *infin, double *correl)
{
    double ret_val, d__1, d__2, d__3, d__4;

    --infin;
    --upper;
    --lower;

    if (infin[1] == 2 && infin[2] == 2) {
        ret_val = mvbvu(&lower[1], &lower[2], correl)
                - mvbvu(&upper[1], &lower[2], correl)
                - mvbvu(&lower[1], &upper[2], correl)
                + mvbvu(&upper[1], &upper[2], correl);
    } else if (infin[1] == 2 && infin[2] == 1) {
        ret_val = mvbvu(&lower[1], &lower[2], correl)
                - mvbvu(&upper[1], &lower[2], correl);
    } else if (infin[1] == 1 && infin[2] == 2) {
        ret_val = mvbvu(&lower[1], &lower[2], correl)
                - mvbvu(&lower[1], &upper[2], correl);
    } else if (infin[1] == 2 && infin[2] == 0) {
        d__1 = -upper[1];
        d__2 = -upper[2];
        d__3 = -lower[1];
        d__4 = -upper[2];
        ret_val = mvbvu(&d__1, &d__2, correl) - mvbvu(&d__3, &d__4, correl);
    } else if (infin[1] == 0 && infin[2] == 2) {
        d__1 = -upper[1];
        d__2 = -upper[2];
        d__3 = -upper[1];
        d__4 = -lower[2];
        ret_val = mvbvu(&d__1, &d__2, correl) - mvbvu(&d__3, &d__4, correl);
    } else if (infin[1] == 1 && infin[2] == 0) {
        d__1 = -upper[2];
        d__2 = -(*correl);
        ret_val = mvbvu(&lower[1], &d__1, &d__2);
    } else if (infin[1] == 0 && infin[2] == 1) {
        d__1 = -upper[1];
        d__2 = -(*correl);
        ret_val = mvbvu(&d__1, &lower[2], &d__2);
    } else if (infin[1] == 1 && infin[2] == 1) {
        ret_val = mvbvu(&lower[1], &lower[2], correl);
    } else if (infin[1] == 0 && infin[2] == 0) {
        d__1 = -upper[1];
        d__2 = -upper[2];
        ret_val = mvbvu(&d__1, &d__2, correl);
    } else {
        ret_val = 1.;
    }
    return ret_val;
}

#undef abs
#undef max

static void
copula_fit_normal(dml_copula_t *copula,
                  const gsl_vector *u,
                  const gsl_vector *v,
                  dml_measure_t *measure)
{
    double *params;
    double rho;
    double tau;

    if (measure == NULL) {
        measure = dml_measure_alloc(u, v);
        tau = dml_measure_tau_coef(measure);
        dml_measure_free(measure);
    } else {
        tau = dml_measure_tau_coef(measure);
    }
    rho = sin(0.5 * M_PI * tau);

    params = copula->data;
    params[0] = rho;
}

static void
copula_pdf_normal(const dml_copula_t *copula,
                  const gsl_vector *u,
                  const gsl_vector *v,
                  gsl_vector *pdf)
{
    double rho;
    double b1, b2, ui, vi, pdfi;
    double *params;

    params = copula->data;
    rho = params[0];

    if (rho == -1 || rho == 1) {
        gsl_vector_set_all(pdf, 0);
    } else {
        for (int i = 0; i < u->size; i++) {
            ui = gsl_vector_get(u, i);
            vi = gsl_vector_get(v, i);

            if (ui <= 0 || ui >= 1 || vi <= 0 || vi >= 1) {
                pdfi = 0;
            } else {
                b1 = gsl_cdf_ugaussian_Pinv(ui);
                b2 = gsl_cdf_ugaussian_Pinv(vi);
                pdfi =
                        gsl_ran_bivariate_gaussian_pdf(b1, b2, 1, 1, rho)
                                / (gsl_ran_ugaussian_pdf(b1)
                                        * gsl_ran_ugaussian_pdf(b2));
            }
            gsl_vector_set(pdf, i, pdfi);
        }
    }
}

static void
copula_cdf_normal(const dml_copula_t *copula,
                  const gsl_vector *u,
                  const gsl_vector *v,
                  gsl_vector *cdf)
{
    double rho;
    double ui, vi, cdfi;
    int infin[] = { 0, 0 };
    double lower[] = { 0, 0 };
    double upper[2];
    double *params;

    params = copula->data;
    rho = params[0];

    for (int i = 0; i < u->size; i++) {
        ui = gsl_vector_get(u, i);
        vi = gsl_vector_get(v, i);

        if (ui <= 0 || vi <= 0) {
            cdfi = 0;
        } else if (ui >= 1) {
            cdfi = vi;
        } else if (vi >= 1) {
            cdfi = ui;
        } else {
            upper[0] = gsl_cdf_ugaussian_Pinv(ui);
            upper[1] = gsl_cdf_ugaussian_Pinv(vi);
            cdfi = mvbvn(&lower[0], &upper[0], &infin[0], &rho);
        }
        gsl_vector_set(cdf, i, cdfi);
    }
}

static void
copula_h_normal(const dml_copula_t *copula,
                const gsl_vector *u,
                const gsl_vector *v,
                gsl_vector *h)
{
    double eps;
    double rho;
    double b1, b2, ui, vi, hi;
    double *params;

    params = copula->data;
    rho = params[0];
    eps = 1e-8;

    for (int i = 0; i < u->size; i++) {
        ui = gsl_vector_get(u, i);
        vi = gsl_vector_get(v, i);

        if (ui <= eps || (rho == 1 && ui == vi && ui != 1)) {
            hi = eps;
        } else if (ui >= 1 - eps || (rho == -1 && 1 - (ui + vi) <= eps)) {
            hi = 1 - eps;
        } else {
            vi = CLAMP(vi, eps, 1 - eps);
            b1 = gsl_cdf_ugaussian_Pinv(ui);
            b2 = gsl_cdf_ugaussian_Pinv(vi);
            hi = gsl_cdf_ugaussian_P((b1 - rho * b2) / sqrt(1 - rho * rho));
            hi = CLAMP(hi, eps, 1 - eps);
        }
        gsl_vector_set(h, i, hi);
    }
}

static void
copula_hinv_normal(const dml_copula_t *copula,
                   const gsl_vector *u,
                   const gsl_vector *v,
                   gsl_vector *hinv)
{
    double eps;
    double rho;
    double b1, b2, ui, vi, hinvi;
    double *params;

    params = copula->data;
    rho = params[0];
    eps = 1e-8;

    for (int i = 0; i < u->size; i++) {
        ui = gsl_vector_get(u, i);

        if (ui <= eps) {
            hinvi = eps;
        } else if (ui >= 1 - eps) {
            hinvi = 1 - eps;
        } else {
            vi = gsl_vector_get(v, i);
            vi = CLAMP(vi, eps, 1 - eps);
            b1 = gsl_cdf_ugaussian_Pinv(ui);
            b2 = gsl_cdf_ugaussian_Pinv(vi);
            hinvi = gsl_cdf_ugaussian_P(b1 * sqrt(1 - rho * rho) + rho * b2);
            hinvi = CLAMP(hinvi, eps, 1 - eps);
        }
        gsl_vector_set(hinv, i, hinvi);
    }
}

static void
copula_aic_normal(const dml_copula_t *copula,
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
    *aic = -2 * loglik + 2;

    gsl_vector_free(pdf);
}

static void
copula_free_normal(dml_copula_t *copula)
{
    g_free(copula->data);
}

dml_copula_t *
dml_copula_alloc_normal(const double rho)
{
    dml_copula_t *copula;
    double *params;

    copula = g_malloc(sizeof(dml_copula_t));
    copula->type = DML_COPULA_NORMAL;
    copula->fit = copula_fit_normal;
    copula->pdf = copula_pdf_normal;
    copula->cdf = copula_cdf_normal;
    copula->h = copula_h_normal;
    copula->hinv = copula_hinv_normal;
    copula->aic = copula_aic_normal;
    copula->free = copula_free_normal;
    copula->data = g_malloc(sizeof(double));
    params = copula->data;
    params[0] = rho;

    return copula;
}
