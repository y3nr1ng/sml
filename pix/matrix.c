#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "pix.h"

void *vdup(void *v, int vlen)
{
    void *r;

    if (v == NULL || vlen == 0)
	return NULL;
    if ((r = malloc(vlen)) == NULL)
	pstop("!!! vdup: not enough memory.\n");
    memcpy(r, v, vlen);

    return r;
}

void vadd_s(int vlen, int op, short *v1, short *v2, short *r)
{
    int  i;

    if (op == 1) {
#pragma omp parallel for private(i)
	for (i=0; i < vlen; i++)
	    r[i] = (short)(v1[i] + v2[i]);
    }
    else if (op == 0) {
#pragma omp parallel for private(i)
	for (i=0; i < vlen; i++)
	    r[i] = (short)((v1[i] > v2[i]) ? (v1[i] - v2[i]) : 0);
    }
    else {
#pragma omp parallel for private(i)
	for (i=0; i < vlen; i++)
	    r[i] = (short)(v1[i] - v2[i]);
    }
}

int vmax_s(int vlen, short *v, int *i_max)
{
    int   i, im;
    short max;

    im  = 0;
    max = v[0];
    for (i=1; i < vlen; i++) {
	if (max < v[i]) {
	    im  = i;
	    max = v[i];
	}
    }
    if (i_max != NULL)
	*i_max = im;

    return (int)max;
}

void vmaxmin_i(int vlen, int *v, int *vmax, int *i_max, int *vmin,
	       int *i_min)
{
    int i, im1, im2;
    int max, min;

    im1 = 0;
    im2 = 0;
    max = v[0];
    min = v[0];
    for (i=1; i < vlen; i++) {
	if (max < v[i]) {
	    im1 = i;
	    max = v[i];
	}
	if (min > v[i]) {
	    im2 = i;
	    min = v[i];
	}
    }
    if (i_max != NULL) *i_max = im1;
    if (i_min != NULL) *i_min = im2;
    *vmax = max;
    *vmin = min;
}

void mx_sub_s(int dim_x, int dim_y, int *sdim, short *mx, short *mr)
{
    int    x1, x2, y1, y2, sdim_x;
    int    i, j, xx, yy;

    x1 = sdim[0];
    x2 = sdim[1];
    y1 = sdim[2];
    y2 = sdim[3];
    sdim_x = x2-x1+1;

#pragma omp parallel for private(i,j,xx,yy)
    for (j=y1; j <= y2; j++) {
    for (i=x1; i <= x2; i++) {
	yy = (i-x1) + (j-y1)*sdim_x;
	xx = i + j*dim_x;
	mr[yy] = mx[xx];
    }}
}

// void mx_subT_s(int dim_x, int dim_y, int *sdim, short *mx, short *mr)
// {
//     int    x1, x2, y1, y2, sdim_y;
//     int    i, j, xx, yy;
// 
//     x1 = sdim[0];
//     x2 = sdim[1];
//     y1 = sdim[2];
//     y2 = sdim[3];
//     sdim_y = y2-y1+1;
// 
// #pragma omp parallel for private(i,j,xx,yy)
//     for (j=y1; j <= y2; j++) {
//     for (i=x1; i <= x2; i++) {
// 	yy = (i-x1)*sdim_y + (j-y1);
// 	if (i >= 0 && j >= 0 && i < dim_x && j < dim_y) { 
// 	    xx = i + j*dim_x;
// 	    mr[yy] = mx[xx];
// 	}
// 	else
// 	    mr[yy] = 0;
//     }}
// }

void mx_rsub_s(int dim_x, int dim_y, int *sdim, short *mx, short *mr)
{
    int  x1, x2, y1, y2;
    int  i, j, xx;

    x1 = sdim[0];
    x2 = sdim[1];
    y1 = sdim[2];
    y2 = sdim[3];

#pragma omp parallel for private(i,j,xx)
    for (j=y1; j <= y2; j++) {
    for (i=x1; i <= x2; i++) {
	xx = i + j*dim_x;
	mr[xx] = mx[xx];
    }}
}

int mx_inv(int N, double *A, int NRHS, double *B)
{
/*
    double *work;
    int    *ipiv, lwork, NB=16, info;
    char   *uplo="U";

    lwork = N*NB;
    work  = malloc(lwork*sizeof(double));
    ipiv  = malloc(N*sizeof(int));
    if (work==NULL || ipiv==NULL)
	pstop("!!! mx_inv: not enough memory.\n");
    dsysv_(uplo, &N, &NRHS, A, &N, ipiv, B, &N, work, &lwork, &info, 1);
    free(work);
    free(ipiv);
*/

    int  info=0;
    dsysv(N, A, NRHS, B, B);

    return info;
}
