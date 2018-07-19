#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "pix.h"

/*-------------------------------------------------------------------------
 *
 *	Compute:  wx(z) or wy(z),   *dy = errbar^2
 *
 *------------------------------------------------------------------------*/

static int func_ca(void *data, double x, double *y, double *dy)
{
    calb3D_t  *ca = (calb3D_t *)data;
    double     w0, A, B, c, d, dw0, dA, dB, dc, dd;
    double     xi, sq, dxi, dsq, ww;

    A   = ca->A;
    B   = ca->B;
    c   = ca->c;
    d   = ca->d;
    dA  = ca->dA;
    dB  = ca->dB;
    dc  = ca->dc;
    dd  = ca->dd;
    w0  = ca->w0;
    dw0 = ca->dw0;

    xi  = (x-c)/d;
    sq  = 1.0 + xi*xi*(1.0 + A*xi + B*xi*xi);
    if (sq <= 0.0) return -1;
    ww  = w0*w0*sq;
    *y  = sqrt(ww);
    if (dy == NULL) return 0;

    dxi = (x*x+c*c)/(d*d) * (dd*dd)/(d*d) + (dc*dc)/(d*d);
    dsq = xi*xi*((4.0 + xi*xi*(9.0*A*A + xi*xi*16.0*B*B))*dxi
	       + xi*xi*xi*xi*(dA*dA + xi*xi*dB*dB));
    *dy = ww*(dw0*dw0/(w0*w0) + dsq/(4.0*sq*sq));

    return 0;
}

/*-------------------------------------------------------------------------
 *
 *	Compute:  wx(z)/wy(z),   *dy = errbar^2
 *
 *------------------------------------------------------------------------*/
/*
static int func_ca2(void *data, double x, double *y, double *dy)
{
    calb3D_t **cc = (calb3D_t **)data;
    calb3D_t  *cx =  cc[0];
    calb3D_t  *cy =  cc[1];
    double    *pdx, *pdy;
    double     yx, yy, dyx, dyy;
    int        r1, r2;

    if (dy == NULL) {
	pdx = NULL;
	pdy = NULL;
    }
    else {
	pdx = &dyx;
	pdy = &dyy;
    }
    r1 = func_ca(cx, x, &yx, pdx);
    r2 = func_ca(cy, x, &yy, pdy);
    if (r1 != 0 || r2 != 0)
	return -1;

    *y = yx/yy;
    if (dy != NULL)
	*dy = (yx*yx)/(yy*yy)*(dyx/(yx*yx) + dyy/(yy*yy));
    return 0;
}
*/
/*-------------------------------------------------------------------------
 *
 *	Bisection solver
 *
 *------------------------------------------------------------------------*/

static int
bisec(int mode, int (*func)(void *, double, double *, double *),
      void *ca, double x0, double dd, double rhs, double drhs, double *x)
{
    int     n=0, loop;
    double *pdy0, *pdy1, x1, y0, y1, dy0, dy1, tol;

    loop = 100;
    tol  = 1.e-10;
    x1   = 0.0;
    y0   = 0.0;
    y1   = 0.0;
    dy0  = 0.0;
    dy1  = 0.0;
    pdy0 = (mode == 0) ? NULL : &dy0;
    pdy1 = (mode == 0) ? NULL : &dy1;

    if (func(ca, x0, &y0, pdy0) != 0) return -1;
    y0  = y0 - rhs;
    dy0 = sqrt(dy0 + drhs*drhs);

    if (mode == 1)
	y0 = y0 + dy0;
    else if (mode == 2)
	y0 = y0 - dy0;

    while (fabs(y0) > tol && n < loop) {
	dd  = dd / 2.0;
	x1  = x0 + dd;
	if (func(ca, x1, &y1, pdy1) != 0) return -1;
	y1  = y1 - rhs;
	dy1 = sqrt(dy1 + drhs*drhs);

	if (mode == 1)
	    y1 += dy1;
	else if (mode == 2)
	    y1 -= dy1;

	if (y0*y1 > 0.0) {
	    x0  = x1;
	    y0  = y1;
	    dy0 = dy1;
	    x1  = x1 + dd;
	}
	n ++;
    }
    *x = x0;

    return 0;
}

static int
solver(double x0, double x1, double dd,
       int (*func)(void *, double, double *, double *),
       void *ca, double rhs, double drhs, int nx, double *x, double *dx)
{
    int    i, npt, n0, n1, n2, r0, r1;
    double y0, y1, dy0, dy1, xr1[16], xr2[16], dx1[16], dx2[16];

    npt = (int)((x1-x0)/dd)+1;
    r0  = func(ca, x0, &y0, &dy0);
    y0  = y0-rhs;
    dy0 = sqrt(dy0 + drhs*drhs);
    n0  = 0;
    n1  = 0;
    n2  = 0;
    for (i=0; i < nx; i++) {
	x[i]   = 0.0;
	dx[i]  = 0.0;
	xr1[i] = 0.0;
	xr2[i] = 0.0;
	dx1[i] = 0.0;
	dx2[i] = 0.0;
    }
    for (i=0; i < npt; i++) {
	r1  = func(ca, x0+dd, &y1, &dy1);
	y1  = y1-rhs;
	dy1 = sqrt(dy1 + drhs*drhs);

	if (r0 == 0 && r1 == 0) {
	    if (n0 < nx && y0*y1 <= 0.0 &&
	        bisec(0, func, ca, x0, dd, rhs, 0.0,  x+n0) == 0) n0 ++;
	    if (n1 < nx && (y0+dy0)*(y1+dy1) <= 0.0 &&
		bisec(1, func, ca, x0, dd, rhs, drhs, xr1+n1) == 0) n1 ++;
	    if (n2 < nx && (y0-dy0)*(y1-dy1) <= 0.0 &&
		bisec(2, func, ca, x0, dd, rhs, drhs, xr2+n2) == 0) n2 ++;
	}
	if (n0 >= nx && n1 >= nx && n2 >= nx) break;

	x0  = x0 + dd;
	y0  = y1;
	dy0 = dy1;
	r0  = r1;
    }

    if (n0 == n1) {
	for (i=0; i < n0; i++)
	    dx1[i] = fabs(xr1[i]-x[i]);
    }
    if (n0 == n2) {
	for (i=0; i < n0; i++)
	    dx2[i] = fabs(xr2[i]-x[i]);
    }
    for (i=0; i < n0; i++)
	dx[i] = (dx1[i] > dx2[i]) ? dx1[i] : dx2[i];

    return n0;
}

/*-------------------------------------------------------------------------
 *
 *	Determine the best z coordinate.
 *
 *------------------------------------------------------------------------*/

static int match_z(para_t *p, fitres_t *res)
{
    int     i, j;
    double  max_del_z, min_zr, zr, zx, zy, dzx, dzy;

    max_del_z = p->max_del_z;
    min_zr  = 1.0;
    res->z  = 0.0;
    res->dz = 0.0;
    zx  = 0.0;
    zy  = 0.0;
    dzx = 0.0;
    dzy = 0.0;

    for (j=0; j < res->n_zy; j++) {
    for (i=0; i < res->n_zx; i++) {
	zr = fabs(1.0 - res->zx[i] / res->zy[j]);
	if (zr < max_del_z && zr < min_zr) {
	    min_zr = zr;
	    zx  = res->zx[i];
	    zy  = res->zy[j];
	    dzx = res->dzx[i];
	    dzy = res->dzy[j];
	}
    }}
    if (min_zr > max_del_z) return -1;

    res->z  = (zx+zy)/2.0;
    res->dz = sqrt(dzx*dzx + dzy*dzy)/2.0;
    return 0;
}

/*-------------------------------------------------------------------------
 *
 *	Solve z from wx/wy
 *
 *------------------------------------------------------------------------*/

int solve_z_w(para_t *p, fitres_t *res)
{
    int       n, i;
    double    z1, z2, dd, wx, wy, dwx, dwy, z[2], dz[2];
    calb3D_t *ca[2];

    z1    = p->z1;
    z2    = p->z2;
    dd    = p->dz;
    ca[0] = &(p->cax);
    ca[1] = &(p->cay);
    wx    = res->a[3]  * p->nm_px_x;
    wy    = res->a[4]  * p->nm_px_y;
    dwx   = res->da[3] * p->nm_px_x;
    dwy   = res->da[4] * p->nm_px_y;

    n = solver(z1, z2, dd, func_ca, ca[0], wx, dwx, 2, z, dz);
    res->n_zx = n;
    for (i=0; i < n; i++) {
	res->zx[i]  = z[i];
	res->dzx[i] = dz[i];
    }
    n = solver(z1, z2, dd, func_ca, ca[1], wy, dwy, 2, z, dz);
    res->n_zy = n;
    for (i=0; i < n; i++) {
	res->zy[i]  = z[i];
	res->dzy[i] = dz[i];
    }
    return match_z(p, res);

/*
    dwr = wx/wy*sqrt((dwx*dwx)/(wx*wx) + (dwy*dwy)/(wy*wy));
    n   = solver(z1, z2, dd, func_ca2, ca, wx/wy, dwr, 2, z, dz);
    res->n_zr = n;
    for (i=0; i < n; i++) {
	res->zr[i]  = z[i];
	res->dzr[i] = dz[i];
    }
*/
}
