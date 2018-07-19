#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "calb.h"

/*-------------------------------------------------------------------------
 *  Read the calibration raw data file.
 */
static int check_nlines(FILE *f)
{
    int  nline=0;
    char buf[BUFLEN];

    fseek(f, 0, SEEK_SET);
    while (fgets(buf, BUFLEN, f) != NULL)
	nline ++;
    fseek(f, 0, SEEK_SET);

    return nline;
}

calb_t *calb_readraw(char *fn)
{
    FILE   *f;
    calb_t *cb;
    char    buf[BUFLEN];
    double  a1, a2, a3;
    int     nline, i, r;

    f = openfile(fn, "rt");
    nline = check_nlines(f);

    cb      = memalloc(sizeof(calb_t), NULL);
    cb->x   = memalloc(sizeof(double)*nline, NULL);
    cb->f   = memalloc(sizeof(double)*nline, NULL);
    cb->err = memalloc(sizeof(double)*nline, NULL);

    for (i=0; i < nline; i++) {
	fgets(buf, BUFLEN, f);
	r = sscanf(buf, "%lf %lf %lf", &a1, &a2, &a3);
	if (r < 2) break;

	cb->x[i]   = a1;
	cb->f[i]   = a2;
	cb->err[i] = (r == 3) ? a3 : 1.0;
    }
    cb->N = i;
    fclose(f);

    return cb;
}

/*-------------------------------------------------------------------------
 *  Destroy the calb_t data structure
 */
void calb_destroy(calb_t *cb)
{
    free(cb->x);
    free(cb->f);
    free(cb->err);
    free(cb);
}

/*-------------------------------------------------------------------------
 *  Initialize the fitting parameters.
 */
void init_para(int verb, calb_t *cb, double *a)
{
    int     i, ii, N;
    double *x, *f, w0, c, d, A, B;
    double  df, dx;

    N  = cb->N;
    x  = cb->x;
    f  = cb->f;
    c  = x[0];
    w0 = f[0];
    ii = 0;
    for (i=1; i < N; i++) {
	if (w0 > f[i]) {
	    c  = x[i];
	    w0 = f[i];
	    ii = i;
	}
    }

    df = ((f[ii-1]-w0) + (f[ii+1]-w0)) / 2.0;
    dx = ((c-x[ii-1])  + (x[ii+1]-c))  / 2.0;
    d  = dx / sqrt(2.0*df/w0);
    A  = 0.5;
    B  = 0.5;

    a[0] = w0;
    a[1] = c;
    a[2] = d;
    a[3] = A;
    a[4] = B;
    if (verb) {
	printf("Init: w0: %13.6E\n", a[0]);
	printf("Init:  A: %13.6E\n", a[3]);
	printf("Init:  B: %13.6E\n", a[4]);
	printf("Init:  c: %13.6E\n", a[1]);
	printf("Init:  d: %13.6E\n", a[2]);
    }
}

/*-------------------------------------------------------------------------
 *  Perform the calibration curve fitting.
 */
void calb_fit(int verb, calb_t *cb, char *name)
{
    int    na, mloop;
    double a[5], da[5], chisq, tol;

    if (verb)
	printf("\nFitting for %s calibration curve ....\n", name);

    na    = 5;
    mloop = 200;
    tol   = 1.0e-6;

    init_para(verb, cb, a);
    nlinfit(cbfunc, cb, cb->N, na, a, da, &chisq, tol, mloop, verb); 
    cb->w0 = a[0];
    cb->c  = a[1];
    cb->d  = a[2];
    cb->A  = a[3];
    cb->B  = a[4];
    cb->dw = da[0];
    cb->dc = da[1];
    cb->dd = da[2];
    cb->dA = da[3];
    cb->dB = da[4];

    if (verb) {
	printf("Res:  w0: %13.6E +- %.4E\n", cb->w0, cb->dw);
	printf("Res:   A: %13.6E +- %.4E\n", cb->A,  cb->dA);
	printf("Res:   B: %13.6E +- %.4E\n", cb->B,  cb->dB);
	printf("Res:   c: %13.6E +- %.4E\n", cb->c,  cb->dc);
	printf("Res:   d: %13.6E +- %.4E\n", cb->d,  cb->dd);
	printf("Res:  chisq/dof: %E\n", chisq);
    }
}

/*-------------------------------------------------------------------------
 *  Output the fitting results
 */
void calb_output(char *outfn, calb_t *cbX, calb_t *cbY)
{
    FILE *f;

    f = openfile(outfn, "wt");
    fprintf(f, "w0x = %17.10E +- %10.4E\n", cbX->w0, cbX->dw);
    fprintf(f, "WxA = %17.10E +- %10.4E\n", cbX->A,  cbX->dA);
    fprintf(f, "WxB = %17.10E +- %10.4E\n", cbX->B,  cbX->dB);
    fprintf(f, "Wxc = %17.10E +- %10.4E\n", cbX->c,  cbX->dc);
    fprintf(f, "Wxd = %17.10E +- %10.4E\n", cbX->d,  cbX->dd);
    fprintf(f, "w0y = %17.10E +- %10.4E\n", cbY->w0, cbY->dw);
    fprintf(f, "WyA = %17.10E +- %10.4E\n", cbY->A,  cbY->dA);
    fprintf(f, "WyB = %17.10E +- %10.4E\n", cbY->B,  cbY->dB);
    fprintf(f, "Wyc = %17.10E +- %10.4E\n", cbY->c,  cbY->dc);
    fprintf(f, "Wyd = %17.10E +- %10.4E\n", cbY->d,  cbY->dd);
    fclose(f);
}
