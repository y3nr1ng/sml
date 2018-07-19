#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "aimg.h"

/*-------------------------------------------------------------------------
 *
 *	Generate parameters of artificial spots.
 *
 *------------------------------------------------------------------------*/

static void gen_spot(para_t *p)
{
    int     i;
    double  r[4], rng;

    for (i=0; i < p->np; i++) {
	r[0] = (double)rand() / (double)RAND_MAX;
	r[1] = (double)rand() / (double)RAND_MAX;
	r[2] = (double)rand() / (double)RAND_MAX;
	r[3] = (double)rand() / (double)RAND_MAX;
	p->sp[i].x  = r[0] * (double)p->iW;
	p->sp[i].y  = r[1] * (double)p->iH;
	rng = p->I2 - p->I1;
	p->sp[i].II = r[2] * rng + p->I1;
	rng = p->w2 - p->w1;
	p->sp[i].w  = r[3] * rng + p->w1;
    }
}

/*-------------------------------------------------------------------------
 *
 *	Generate spot pixels on the output image.
 *
 *------------------------------------------------------------------------*/

static double sp_gaussian(spot_t *sp, double x, double y)
{
    double  xx, yy, w;

    xx = x - sp->x;
    yy = y - sp->y;
    w  = sp->w;
    return (sp->II * exp(-(xx*xx+yy*yy) / (2.0*w*w)) * MAX_PIXEL);
}

static void gen_spot_mesh(para_t *p, spot_t *sp)
{
    int     lmesh, nmesh, x0, y0;
    int     i, j, k, ii, jj, ix, iy, xx;
    double *spot, II, x, y, tic;

    lmesh = p->spixel * p->smesh;
    nmesh = p->smesh  * p->smesh;
    tic   = 1.0 / (double)p->smesh;
    x0    = (int)round(sp->x);
    y0    = (int)round(sp->y);
    if ((spot = calloc(lmesh*lmesh, sizeof(double))) == NULL)
	pstop("!!! gen_spot_mesh: not enough memory.\n");

    for (j=0; j < lmesh; j++) {
    for (i=0; i < lmesh; i++) {
	k = i + j*lmesh;
	x = (double)(i-lmesh/2) * tic + sp->x;
	y = (double)(j-lmesh/2) * tic + sp->y;
	spot[k] = sp_gaussian(sp, x, y);
    }}
    for (jj=0; jj < p->spixel; jj++) {
    for (ii=0; ii < p->spixel; ii++) {
	ix = x0 + ii - p->spixel/2;
	iy = y0 + jj - p->spixel/2;
	if (ix >= 0 && ix < p->iW && iy >= 0 && iy < p->iH) {
	    II = 0.0;
	    for (j=0; j < p->smesh; j++) {
	    for (i=0; i < p->smesh; i++) {
	        k   = (i+ii*p->smesh) + (j+jj*p->smesh)*lmesh;
	        II += spot[k];
	    }}
	    II /= (double)nmesh;

	    xx = ix + iy * p->iW;
	    p->img[xx] += ((unsigned char)II);
	}
    }}
    free(spot);
}

void gen_spot_pixels(para_t *p)
{
    int     i;
    spot_t *sp = p->sp;

    for (i=0; i < p->np; i++) {
	gen_spot_mesh(p, sp);
	sp ++;
    }
}

/*-------------------------------------------------------------------------
 *
 *	Generate the background noise on the image.
 *
 *------------------------------------------------------------------------*/

void gen_noise(para_t *p)
{
    int     i, imglen;
    double  r, rr, rng, ra;

    imglen = p->iW * p->iH;
    rng    = p->N2 - p->N1;
    ra     = 0.0;
    for (i=0; i < imglen; i++) {
	r  = ((double)rand() / (double)RAND_MAX * rng + p->N1);
	rr = (double)(p->img[i]) + r * MAX_PIXEL;
	p->img[i] = (unsigned char)((rr > MAX_PIXEL) ? MAX_PIXEL : rr);
	ra += r;
    }
    p->noise = ra / (double)imglen;
}

/*-------------------------------------------------------------------------
 *
 *	Generate the artificial image with spots.
 *
 *------------------------------------------------------------------------*/

void aimg(para_t *p)
{
    int  imglen;

    imglen = p->iW * p->iH;
    if ((p->img = calloc(imglen, sizeof(char))) == NULL)
	pstop("!!! aimg: not enough memory for image.\n");
    if ((p->sp = calloc(p->np, sizeof(spot_t))) == NULL)
	pstop("!!! aimg: not enough memory for spot list.\n");

    gen_spot(p);
    gen_spot_pixels(p);
    gen_noise(p);
}
