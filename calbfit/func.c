#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calb.h"

/*-------------------------------------------------------------------------
 *  Calibration function:
 *
 *  f(x)   = w0*sqrt(1+((x-c)/d)**2*(1+A*((x-c)/d)+B*((x-c)/d)**2))
 *
 *  df(x)  = (w0**2/(2*f(x)))/d*((x-c)/d)*(2+3*A*((x-c)/d)+4*B*((x-c)/d)**2)
 *
 *  ddf(x) = (-(df(x))**2+(w0/d)**2*(1+3*A*((x-c)/d)+6*B*((x-c)/d)**2))/f(x)
 *
 *
 *  df/dw0 = sqrt(1+((x-c)/d)**2*(1+A*((x-c)/d)+B*((x-c)/d)**2))
 *
 *  df/dc  =-(w0**2/(2*f(x)))/d*((x-c)/d)*(2+3*A*((x-c)/d)+4*B*((x-c)/d)**2)
 *
 *  df/dd  =-(w0**2/(2*f(x)))/d*((x-c)/d)**2*(2+3*A*((x-c)/d)+4*B*((x-c)/d)**2)
 *
 *  df/dA  = (w0**2/(2*f(x))) * ((x-c)/d)**3
 *
 *  df/dB  = (w0**2/(2*f(x))) * ((x-c)/d)**4
 *
 */

/*-------------------------------------------------------------------------
 *      Compute:
 *
 *          f[x] = (func[x]-data[x]) / sqrt(sig[x])
 *
 *       dfda[x] = (df/da)[x] / sqrt(sig[x])
 */
void cbfunc(void *data, double *a, double *f, double *dfda)
{
    calb_t *dd = (calb_t *)data;
    double *x0, *f0, *err;
    double  w0, c, d, A, B;
    double  ff, xx, sq, dsq, fac;
    int     N, i;

    N   = dd->N;
    x0  = dd->x;
    f0  = dd->f;
    err = dd->err;
    w0  = a[0];
    c   = a[1];
    d   = a[2];
    A   = a[3];
    B   = a[4];

    for (i=0; i < N; i++) {
	xx   = (x0[i]-c)/d;
	sq   = sqrt(1.0 + xx*xx*(1.0 + A*xx + B*xx*xx));
	ff   = w0*sq;
	f[i] = (ff - f0[i])/err[i];

	fac  = w0*w0/(2.0*ff);
	dsq  = 2.0 + 3.0*A*xx + 4.0*B*xx*xx;
	dfda[i*5+0] = sq / err[i];
	dfda[i*5+1] = (-fac/d * xx * dsq) / err[i];
	dfda[i*5+2] = (-fac/d * xx*xx * dsq) / err[i];
	dfda[i*5+3] = (fac * xx*xx*xx) / err[i];
	dfda[i*5+4] = (fac * xx*xx*xx*xx) / err[i];
    }
}

/*-------------------------------------------------------------------------
 *      Compute:
 *
 *          f[x] = func[x];
 *
 *         df[x] = (df/dx)[x];
 *
 *        ddf[x] = (d^2 f/dx^2)[x];
 */
void cbddf(void *data, double *f, double *df, double *ddf)
{
    calb_t *dd = (calb_t *)data;
    double *x0, w0, c, d, A, B;
    double  ff, gg, xx, sq, dsq, ddq, fac;
    int     N, i;

    N   = dd->N;
    x0  = dd->x;
    w0  = dd->w0;
    c   = dd->c;
    d   = dd->d;
    A   = dd->A;
    B   = dd->B;

    for (i=0; i < N; i++) {
	xx  = (x0[i]-c)/d;
	sq  = sqrt(1.0 + xx*xx*(1.0 + A*xx + B*xx*xx));
	ff  = w0*sq;
	fac = w0*w0/(2.0*ff);
	dsq = 2.0 + 3.0*A*xx + 4.0*B*xx*xx;
	ddq = (w0/d)*(w0/d)*(1.0 + 3.0*A*xx + 6.0*B*xx*xx);
	gg  = fac/d * xx * dsq;

	if (f   != NULL) f[i]   = ff;
	if (df  != NULL) df[i]  = gg;
	if (ddf != NULL) ddf[i] = (-gg*gg + ddq) / ff;
    }
}
