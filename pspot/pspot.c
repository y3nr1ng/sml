/*-------------------------------------------------------------------------
 *
 *	Plot and generate JPEG image for spots
 *
 *------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "pspot.h"

/*-------------------------------------------------------------------------
 *
 *	Read the frame statistics data file.
 *
 *------------------------------------------------------------------------*/

void readfsts(para_t *p)
{
    FILE   *f;
    char    s[1024];
    int     fID, n_evt;
    double  density;

    if ((f = fopen(p->inpfnF, "rt")) == NULL)
	pstop("!!! readfsts: cannot open file: %s\n", p->inpfnF);
    fgets(s, 1024, f);
    fgets(s, 1024, f);
    while (fgets(s, 1024, f) != NULL) {
	sscanf(s, "%d %d %lf", &fID, &n_evt, &density);
	if (density > p->m_density)
	    p->Fstart = fID+1;
    }
    fclose(f);
}

/*-------------------------------------------------------------------------
 *
 *	Read the spot source data file.
 *
 *------------------------------------------------------------------------*/
/*
static void push_spot(double x, double y, para_t *p)
{
    if (p->np <= p->mp) {
	p->mp += 128;
	p->x = realloc(p->x, p->mp*sizeof(double));
	p->y = realloc(p->y, p->mp*sizeof(double));
	if (p->x == NULL || p->y == NULL)
	    pstop("push_spot: not enough memory.\n");
    }
    p->x[p->np] = x;
    p->y[p->np] = y;
    p->np ++;
}
*/

void readspot(para_t *p)
{
    FILE   *f;
    char    s[1024];
    int     n, i, sID, fID, xp, yp, cnt;
    double  II, dII, x, y, dx, dy;

    printf("pspot: starting frame: %d\n", p->Fstart);
    fflush(stdout);
    if ((f = fopen(p->inpfn, "rt")) == NULL)
	pstop("!!! readspot: cannot open file: %s\n", p->inpfn);
    fgets(s, 1024, f);
    fgets(s, 1024, f);
    n = 0;
    while (fgets(s, 1024, f) != NULL) n++;
    fseek(f, 0, SEEK_SET);

    p->x = malloc(n*sizeof(double));
    p->y = malloc(n*sizeof(double));
    for (i=0; i < n; i++) {
	if (fgets(s, 1024, f) == NULL)
	    pstop("!!! readspot: reading data error: %s\n", p->inpfn);
	sscanf(s, "%d %d %d %d %d %lf %lf %lf %lf %lf %lf",
		   &sID, &fID, &xp, &yp, &cnt, &II, &dII, &x, &dx, &y, &dy);
	if (fID >= p->Fstart) {
    	    p->x[p->np] = x;
    	    p->y[p->np] = y;
	    p->np ++;
	}
    }
    fclose(f);
    printf("pspot: total spots read: %d\n", p->np);

    if (p->imgscale == 0) {
	p->imgWidth  = (int)p->iW;
	p->imgHeight = (int)p->iH;
    }
}

/*-------------------------------------------------------------------------
 *
 *	Read the pixel data file.
 *
 *------------------------------------------------------------------------*/

static void push_pixel(int x, int y, int pixel, para_t *p)
{
    if (p->np >= p->mp) {
	p->mp += 1024;
	p->x   = realloc(p->x, p->mp * sizeof(double));
	p->y   = realloc(p->y, p->mp * sizeof(double));
	p->pix = realloc(p->pix, p->mp * sizeof(double));
	if (p->x == NULL || p->y == NULL || p->pix == NULL)
	    pstop("push_pixel: not enough memory.\n");
    }
    p->x[p->np]   = x;
    p->y[p->np]   = y;
    p->pix[p->np] = (double)pixel;
    p->np ++;    
}

void readpixel(para_t *p)
{
    FILE   *f;
    char    s[1024];
    int     x, y, pixel, iW, iH, pmax, pmin;

    if ((f = fopen(p->inpfn, "rt")) == NULL)
	pstop("!!! readpixel: cannot open file: %s\n", p->inpfn);
    fgets(s, 1024, f);
    fgets(s, 1024, f);

    iW   = 0;
    iH   = 0;
    pmax = 0;
    pmin = INT_MAX;
    while (fgets(s, 1024, f) != NULL) {
	sscanf(s, "%d %d %d", &x, &y, &pixel);
	push_pixel(x, y, pixel, p);
	iW   = (iW   < x)     ? x : iW;
	iH   = (iH   < y)     ? y : iH;
	pmax = (pmax < pixel) ? pixel : pmax;
	pmin = (pmin > pixel) ? pixel : pmin;
    }
    p->iW   = (double)(iW+1);
    p->iH   = (double)(iH+1);
    p->pmax = (double)pmax;
    p->pmin = (double)pmin;
    if (p->imgscale == 0) {
	p->imgWidth  = (int)p->iW;
	p->imgHeight = (int)p->iH;
    }
}

/*-------------------------------------------------------------------------
 *
 *	Fill the spots in a image buffer.
 *
 *------------------------------------------------------------------------*/

void spot_image(para_t *p)
{
    unsigned char *img;
    double x, y, iW, iH;
    int    i, px, py, pp, npx;

    npx = p->imgWidth * p->imgHeight;
    if ((img = malloc(npx*3)) == NULL)
	pstop("!!! spot_image: not enough memory.\n");
    memset(img, 255, npx*3);
    iW = p->iW;
    iH = p->iH;
    for (i=0; i < p->np; i++) {
        x  = p->x[i];
        y  = p->y[i];
	if (x >= iW || x < 0.0 || y >= iH || y < 0.0) continue;

        px = (int)(x/iW * p->imgWidth);
        py = p->imgHeight - (int)(y/iH * p->imgHeight);
	if (px < p->imgWidth  && px >= 0 &&
	    py < p->imgHeight && py >= 0) {
	    pp = px*3 + py*3*p->imgWidth;
	    img[pp  ] = (unsigned char)0;
	    img[pp+1] = (unsigned char)0;
	}
    }
    if (p->gimg) free(p->gimg);
    p->gimg = img;
}

/*-------------------------------------------------------------------------
 *
 *	Count the event in a given mesh.
 *
 *------------------------------------------------------------------------*/

void mesh_events(para_t *p)
{
    int     nx, ny, np, ix, iy, ii, i, imax, *img;
    double *x, *y, *pix, mesh_size;

    mesh_size = p->mesh_size;
    nx  = (int)floor(p->iW / mesh_size) + 1;
    ny  = (int)floor(p->iH / mesh_size) + 1;
    np  = nx*ny;
    img = calloc(np, sizeof(int));
    pix = malloc(np * sizeof(double));
    x   = malloc(np * sizeof(double));
    y   = malloc(np * sizeof(double));
    if (! img || ! pix || ! x || ! y)
	pstop("!!! mesh_events: not enough memory.\n");

    for (i=0; i < p->np; i++) {
	ix = (int)floor(p->x[i] / mesh_size);
	iy = (int)floor(p->y[i] / mesh_size);
	if (ix >= 0 && ix < nx && iy >= 0 && iy < ny) {
	    ii = ix + iy*nx;
	    img[ii]++;
	}
    }
    imax = 0;
    for (iy=0; iy < ny; iy++) {
    for (ix=0; ix < nx; ix++) {
	ii = ix + iy*nx;
	x[ii] = (ix+0.5)*mesh_size;
	y[ii] = (iy+0.5)*mesh_size;
	imax  = (imax < img[ii]) ? img[ii] : imax;
	pix[ii] = (double)img[ii];
    }}

    if (p->x)   free(p->x);
    if (p->y)   free(p->y);
    if (p->pix) free(p->pix);
    free(img);
    p->pix  = pix;
    p->x    = x;
    p->y    = y;
    p->np   = np;
    p->mp   = np;
    p->pmax = (double)imax;
    p->pmin = 0.0;
    if (p->imgscale == 0) {
	p->imgWidth  = nx;
	p->imgHeight = ny;
    }
}

void mesh_output(para_t *p)
{
    FILE *f;
    int   i;

    if ((f = fopen(p->outmesh, "wt")) == NULL) {
	printf("!!! mesh_output: cannot open file: %s\n", p->outmesh);
	exit(1);
    }
    fprintf(f, "%12s  %12s  %8s\n", "x", "y", "n_events");
    fprintf(f, "------------------------------------\n");
    for (i=0; i < p->np; i++) {
	if (p->pix[i] == 0.0) continue;
	fprintf(f, "%12.4E  %12.4E  %8d\n", p->x[i], p->y[i], (int)(p->pix[i]));
    }
    fclose(f);
}

/*-------------------------------------------------------------------------
 *
 *	Fill the pixels in a gray scale image buffer.
 *
 *------------------------------------------------------------------------*/

void gray_image(para_t *p)
{
    int    i, j, k, px, py, pp, npx;
    double x, y, iW, iH;
    double scale, *img2=NULL, pix;

    npx = p->imgWidth * p->imgHeight;
    if (p->imgscale == 1) {
	if ((img2 = calloc(npx, sizeof(double))) == NULL)
	    pstop("!!! gray_image: not enough memory.\n");
    }
    if ((p->gimg = malloc(npx)) == NULL)
	pstop("!!! gray_image: not enough memory.\n");

    iW = p->iW;
    iH = p->iH;
    if (p->imgscale == 1) {
	scale = 0.0;
	for (i=0; i < p->np; i++) {
	    x = p->x[i];
	    y = p->y[i];
	    if (x >= iW || x < 0.0 || y >= iH || y < 0.0) continue;

            px = (int)(x/iW * p->imgWidth);
            py = (int)(p->imgHeight - y/iH * p->imgHeight);
	    pp = px + py*p->imgWidth;
	    img2[pp] += p->pix[i];
	    scale = (img2[pp] > scale) ? img2[pp] : scale;
	}
	for (i=0; i < npx; i++) {
	    pix = round(pow(img2[i]/scale, p->gamma) * 255.0);
	    p->gimg[i] = (unsigned char)pix;
	}
    }
    else {
	scale = p->pmax;
	for (j=0; j < p->imgHeight; j++) {
	for (i=0; i < p->imgWidth; i++) {
	    k = i + j * p->imgWidth;
	    x = p->x[k];
	    y = p->y[k];
	    if (x >= iW || x < 0.0 || y >= iH || y < 0.0) continue;

	    pix = round(pow(p->pix[k]/scale, p->gamma) * 255.0);
	    k = i + (p->imgHeight - j - 1) * p->imgWidth;
	    p->gimg[k] = (unsigned char)pix;
	}}
    }
    if (img2) free(img2);
}
