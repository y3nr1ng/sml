#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "pix.h"

/*-------------------------------------------------------------------------
*
*	Bisection solver
*
*------------------------------------------------------------------------*/

static int bisection(double (*func)(void *, double, int *), void *fd,
                     double x0, double y0, double dx, double tol, double *x) {
    int info=0, n=0;
    double x1, y1;

    while (fabs(y0) > tol && n < 100) {
        dx = dx / 2.0;
        x1 = x0 + dx;
        y1 = func(fd, x1, &info);
        if (info != 0)
            return -1;
        if (y0*y1 > 0.0) {
            x0 = x1;
            y0 = y1;
            x1 = x1 + dx;
        }
        n++;
    }
    *x = x0;
    return 0;
}

static int solver(double x0, double x1, double dx,
                  double (*func)(void *, double, int *),
                  double (*dfunc)(void *, double),
                  void *fd, int nx, double *x, double *xdev) {
    int i, npt, nsol, info0, info1;
    double y0, y1, tol;

    tol  = 1.e-10;
    nsol = 0;

    npt = (int)((x1-x0)/dx)+1;
    y0  = func(fd, x0, &info0);
    for (i=0; i < npt; i++) {
        y1 = func(fd, x0+dx, &info1);
        if (info0 == 0 && info1 == 0 && nsol < nx && y0*y1 <= 0.0) {
            if (bisection(func, fd, x0, y0, dx, tol, x+nsol) == 0) {
                xdev[nsol] = dfunc(fd, x[nsol]);
                nsol++;
            }
        }
        x0 = x0 + dx;
        y0 = y1;
        info0 = info1;
    }
    return nsol;
}

/*-------------------------------------------------------------------------
*
*	Solve z from wx or wy, respectively
*
*------------------------------------------------------------------------*/

static double func_ca(void *data, double x, int *info) {
    calb3D_t      *dd = (calb3D_t *)data;
    double w0, y0, tmp;
    complex double A, B, c, d, zz;

    A   = dd->A[0] + dd->A[1]*I;
    B   = dd->B[0] + dd->B[1]*I;
    c   = dd->c[0] + dd->c[1]*I;
    d   = dd->d[0] + dd->d[1]*I;
    w0  = dd->w0;
    y0  = dd->WW;

    zz  = (x-c)/d;
    tmp = creal(csqrt(1.0 + zz*zz*(1.0 + A*zz + B*zz*zz)));
    *info = 0;
    return (w0*tmp-y0);
}

static double func_dca(void *data, double x) {
    calb3D_t      *dd = (calb3D_t *)data;
    double w0, y0, dy, dz, tmp;
    complex double A, B, c, d, z1, z2, z3, zz;

    w0  = dd->w0;
    y0  = dd->WW;
    dy  = dd->dW;
    A   = dd->A[0] + dd->A[1]*I;
    B   = dd->B[0] + dd->B[1]*I;
    c   = dd->c[0] + dd->c[1]*I;
    d   = dd->d[0] + dd->d[1]*I;

    zz  = (x-c)/d;
    z1  = zz/d;
    z2  = 3.0*A*zz;
    z3  = 4.0*B*zz*zz;
    tmp = sqrt(creal(z1*conj(z1) * (4.0 + z2*conj(z2) + z3*conj(z3))));
    dz  = fabs(2.0*y0*dy) / (w0*w0*tmp);
    return dz;
}

int solve_z_w(para_t *p, calb3D_t *fd, double w, double dw, char *v, char *dv) {
    int n;
    double x1, x2, dx, z[4], dz[4];

    x1 = p->z1;
    x2 = p->z2;
    dx = p->dz;

    fd->WW = w;
    fd->dW = dw;
    n = solver(x1, x2, dx, func_ca, func_dca, fd, 4, z, dz);

    if (n > 0) {
        sprintf(v,  "%13.6E", z[0]);
        sprintf(dv, "%10.3E", dz[0]);
        return 0;
    }else  {
        sprintf(v,  "%13s", "n/a");
        sprintf(dv, "%10s", "n/a");
        return -1;
    }
}

/*-------------------------------------------------------------------------
*
*	Solve z from wx/wy
*
*------------------------------------------------------------------------*/

static double func_ca2(void *data, double x, int *info) {
    calb3D_t     **DD = (calb3D_t **)data;
    calb3D_t      *DX =  DD[0];
    calb3D_t      *DY =  DD[1];
    double y0, tx, ty;
    complex double A, B, c, d, zz;

    A  = DX->A[0] + DX->A[1]*I;
    B  = DX->B[0] + DX->B[1]*I;
    c  = DX->c[0] + DX->c[1]*I;
    d  = DX->d[0] + DX->d[1]*I;
    zz = ((complex double)x-c)/d;
    tx = creal(csqrt(1.0 + zz*zz*(1.0 + A*zz + B*zz*zz)));

    A  = DY->A[0] + DY->A[1]*I;
    B  = DY->B[0] + DY->B[1]*I;
    c  = DY->c[0] + DY->c[1]*I;
    d  = DY->d[0] + DY->d[1]*I;
    zz = ((complex double)x-c)/d;
    ty = creal(csqrt(1.0 + zz*zz*(1.0 + A*zz + B*zz*zz)));

    y0 = DX->WW;
    *info = 0;
    return (tx/ty-y0);
}

static double func_dca2(void *data, double x) {
    calb3D_t     **DD = (calb3D_t **)data;
    calb3D_t      *DX =  DD[0];
    calb3D_t      *DY =  DD[1];
    double y0, dy, err;
    complex double A, B, c, d, zz, z1, z2, z3, d1, d2, y1, y2;

    A  = DX->A[0] + DX->A[1]*I;
    B  = DX->B[0] + DX->B[1]*I;
    c  = DX->c[0] + DX->c[1]*I;
    d  = DX->d[0] + DX->d[1]*I;
    zz = ((complex double)x-c)/d;
    z1 = zz/d;
    z2 = 3.0*A*zz;
    z3 = 4.0*B*zz*zz;
    d1 = z1*conj(z1) * (4.0 + z2*conj(z2) + z3*conj(z3));
    y1 = 1.0 + zz*zz*(1.0 + A*zz + B*zz*zz);

    A  = DY->A[0] + DY->A[1]*I;
    B  = DY->B[0] + DY->B[1]*I;
    c  = DY->c[0] + DY->c[1]*I;
    d  = DY->d[0] + DY->d[1]*I;
    zz = ((complex double)x-c)/d;
    z1 = zz/d;
    z2 = 3.0*A*zz;
    z3 = 4.0*B*zz*zz;
    d2 = z1*conj(z1) * (4.0 + z2*conj(z2) + z3*conj(z3));
    y2 = 1.0 + zz*zz*(1.0 + A*zz + B*zz*zz);

    y0  = DX->WW;
    dy  = DX->dW;
    err = y0 * 0.5 * sqrt(creal(d1/(y1*conj(y1)) + d2/(y2*conj(y2))));
    return (dy/err);
}

int solve_z_wxowy(para_t *p, double *a, double *da, char *v, char *dv) {
    int n;
    double x1, x2, dx, wx, wy, dwx, dwy, z[4], dz[4];
    calb3D_t *fd[2];

    fd[0] = (calb3D_t *)(&p->cax);
    fd[1] = (calb3D_t *)(&p->cay);
    x1    = p->z1;
    x2    = p->z2;
    dx    = p->dz;
    wx    = a[3];
    wy    = a[4];
    dwx   = da[3];
    dwy   = da[4];

    fd[0]->WW = wx/wy;
    fd[0]->dW = wx/wy * sqrt((dwx*dwx)/(wx*wx) + (dwy*dwy)/(wy*wy));
    n = solver(x1, x2, dx, func_ca2, func_dca2, fd, 4, z, dz);

    if (n > 0) {
        sprintf(v,  "%13.6E", z[0]);
        sprintf(dv, "%10.3E", dz[0]);
        return 0;
    }else  {
        sprintf(v,  "%13s", "n/a");
        sprintf(dv, "%10s", "n/a");
        return -1;
    }
}
