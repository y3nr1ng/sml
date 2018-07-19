#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "pix.h"

/*-------------------------------------------------------------------------
 *
 *  Functions to fit the image pixels to locate spot coordinate.
 *
 *------------------------------------------------------------------------*/

typedef struct {
    double *x_fit;
    double *y_fit;
    int    *intensity;
} fitIntensity_t;

//
//  r0: full width of the spots
//
static void
fit_gaussian2D(int n_data, void *data, double *a, double *f, double *dfda)
{
    fitIntensity_t *d = (fitIntensity_t *)data;
    double  I0, x0, y0, r0, B0;
    double *x, *y, dx, dy, ee;
    int    *II, i;

    I0 = a[0];
    x0 = a[1];
    y0 = a[2];
    r0 = a[3];
    B0 = a[4];
    x  = d->x_fit;
    y  = d->y_fit;
    II = d->intensity;

// Note that the weighting "sig" is set to 1, so it does not appear here.
// #pragma omp parallel for private(i,dx,dy,ee)
    for (i=0; i < n_data; i++) {
	dx   = x[i]-x0;
	dy   = y[i]-y0;
	ee   = exp(-2.0*(dx*dx+dy*dy)/(r0*r0));
	f[i] = I0*ee+B0 - (double)II[i];
	dfda[i*5+0] = ee;
	dfda[i*5+1] = I0*ee*(4.0*dx/(r0*r0));
	dfda[i*5+2] = I0*ee*(4.0*dy/(r0*r0));
	dfda[i*5+3] = I0*ee*(4.0*(dx*dx+dy*dy)/(r0*r0*r0));
	dfda[i*5+4] = 1.0;
    }
}

//
//  wx, wy: full width of the spots
//
static void
fit_gaussian3D(int n_data, void *data, double *a, double *f, double *dfda)
{
    fitIntensity_t *d = (fitIntensity_t *)data;
    double  I0, x0, y0, wx, wy, B0;
    double *x, *y, dx, dy, ee;
    int    *II, i;

    I0 = a[0];
    x0 = a[1];
    y0 = a[2];
    wx = a[3];
    wy = a[4];
    B0 = a[5];
    x  = d->x_fit;
    y  = d->y_fit;
    II = d->intensity;

// Note that the weighting "sig" is set to 1, so it does not appear here.
// #pragma omp parallel for private(i,dx,dy,ee)
    for (i=0; i < n_data; i++) {
	dx   = x[i]-x0;
	dy   = y[i]-y0;
	ee   = exp(-(2.0*dx*dx)/(wx*wx) - (2.0*dy*dy)/(wy*wy));
	f[i] = I0*ee+B0 - (double)II[i];
	dfda[i*6+0] = ee;
	dfda[i*6+1] = I0*ee*(4.0*dx/(wx*wx));
	dfda[i*6+2] = I0*ee*(4.0*dy/(wy*wy));
	dfda[i*6+3] = I0*ee*((4.0*dx*dx)/(wx*wx*wx));
	dfda[i*6+4] = I0*ee*((4.0*dy*dy)/(wy*wy*wy));
	dfda[i*6+5] = 1.0;
    }
}

/*-------------------------------------------------------------------------
 *
 *  Collect the fitting results.
 *
 *------------------------------------------------------------------------*/

static void fit_result(sp_t *sp, double *a, double *da, double chisq)
{
    int  i;

    if ((sp->res = calloc(sizeof(fitres_t), 1)) == NULL)
	pstop("!!! fit_result: not enough memory.\n");
    for (i=0; i < 6; i++) {
	sp->res->a[i]  = a[i];
	sp->res->da[i] = da[i];
	sp->res->chisq = chisq;
    }
}

/*-------------------------------------------------------------------------
 *
 *  Fit the Spots coordinates.
 *
 *------------------------------------------------------------------------*/

int SpotFit(para_t *p, sp_t *sp, framests_t *fsts)
{
    static double *x_fit, *y_fit;
    int    i, k, x, y, x_rng, y_rng, na, mloop, sdim_x, sdim_y, imglen;
    int   *intensity, imax, imin;
    double a[10], da[10], chisq, tol;
    fitIntensity_t fdata;
    void (*fitfunc)(int, void *, double *, double *, double *);

    if (sp == NULL) return -1;
    sp->res = NULL;

// Initial parameters for fitting.
    intensity = sp->img;
    sdim_x = p->x_find_pixels;
    sdim_y = p->y_find_pixels;
    imglen = sdim_x * sdim_y;
    x_rng  = sdim_x / 2;
    y_rng  = sdim_y / 2;
    mloop  = 100;			// maxloop for fitting.
    tol    = 1.E-5;			// stopping criteria for fitting.
    if (p->mode == 0) {
	na      = 5;
	fitfunc = fit_gaussian2D;
    }
    else {
	na      = 6;
	fitfunc = fit_gaussian3D;
    }

// Prepare coord. of pixels of the spot, the coord. are relative to center.
#ifdef USE_OMP
#pragma omp critical
#endif
    if (x_fit == NULL) {
	x_fit = malloc(imglen * sizeof(double));
	y_fit = malloc(imglen * sizeof(double));
	if (! x_fit || ! y_fit)
	    pstop("!!! SpotFit: not enough memory.\n");
	k = 0;
	for (y=-y_rng; y <= y_rng; y++) {
	for (x=-x_rng; x <= x_rng; x++) {
	    x_fit[k] = x;
	    y_fit[k] = y;
	    k ++;
	}}
    }

// Prepare initial parameters for fitting.
    vmaxmin_i(imglen, intensity, &imax, NULL, &imin, NULL);
    if (p->mode == 0) {
	a[0] = (double)(imax-imin);	// initial parameter: intensity.
	a[1] = 0.0;			// initial parameter: x0.
	a[2] = 0.0;			// initial parameter: y0.
	a[3] = 1.0;			// initial parameter: Gaussian width.
	a[4] = (double)imin;		// initial parameter: background.
    }
    else {
	a[0] = (double)(imax-imin);	// initial parameter: intensity.
	a[1] = 0.0;			// initial parameter: x0.
	a[2] = 0.0;			// initial parameter: y0.
	a[3] = 1.0;			// initial parameter: wx.
	a[4] = 1.0;			// initial parameter: wy.
	a[5] = (double)imin;		// initial parameter: background.
    }
    fdata.x_fit     = x_fit;
    fdata.y_fit     = y_fit;
    fdata.intensity = intensity;
    if (nlinfit(imglen, fitfunc, &fdata, na, a, da, &chisq, tol,
		mloop, p->verb) != 0)
	return -1;

// Adjust the spot coordinates and parameters from the fitting results.
    a[1] = (double)(sp->x) + a[1];
    a[2] = (double)(sp->y) + a[2];
    if (p->mode == 0) {
	if (a[0]<0.0 || a[1]<0.0 || a[2]<0.0 || a[4]<0.0) return -1;
	if (a[1]  > p->max_x  || a[2]  > p->max_y)        return -1;
	if (da[1] > p->max_dx || da[2] > p->max_dy || da[3] > p->max_dwx ||
	    a[0]/a[4] < p->min_SN || da[0]/a[0] > p->max_dI_I) return -1;
	a[3] = fabs(a[3]);
    }
    else {
        if (a[0]<0.0 || a[1]<0.0 || a[2]<0.0 || a[5]<0.0) return -1;
	if (a[1]  > p->max_x   || a[2]  > p->max_y)        return -1;
	if (da[1] > p->max_dx  || da[2] > p->max_dy ||
	    da[3] > p->max_dwx || da[4] > p->max_dwy ||
	    a[0]/a[5] < p->min_SN || da[0]/a[0] > p->max_dI_I) return -1;
	a[3] = fabs(a[3]);
	a[4] = fabs(a[4]);
    }

    fit_result(sp, a, da, chisq);
#ifdef USE_OMP
#pragma omp critical
#endif
    if (fsts != NULL) {
	i = sp->fID - p->frameID1;
	fsts[i].n_event ++;
    }
    return 0;
}

