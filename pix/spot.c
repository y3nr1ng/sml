#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_OMP
#include <omp.h>
#endif
#include "pix.h"

/*-------------------------------------------------------------------------
 *
 *	Used to sort spots in descent of frame IDs.
 *
 *------------------------------------------------------------------------*/

int spot_cmp(const void *a, const void *b)
{
    sp_t *aa = *((sp_t **)a);
    sp_t *bb = *((sp_t **)b);

    if (aa != NULL && bb != NULL) {
	if (aa->fID > bb->fID)
	    return -1;
	else if (aa->fID < bb->fID)
	    return  1;
	else
	    return  0;
    }
    else {
	if (aa == NULL && bb == NULL)
	    return  0;
	else if (aa == NULL)
	    return  1;
	else
	    return -1;
    }
}

/*-------------------------------------------------------------------------
 *
 *	Spot finding: from a set of successive frames.
 *
 *------------------------------------------------------------------------*/

void spot_dframe(para_t *p)
{
    frameIO_t  *fio;
    frameloc_t *fm0=NULL, *fm1=NULL;
    int         fID, fID0, fIDN, oneprec, n;

    fID0 = p->frameID1;
    fIDN = p->frameID2;
    fio  = frameIO_Open(p);
    oneprec = (fIDN-fID0+1) / 100;

    fm1 = frameIO(p, fio, fIDN);		// load the *next* frame
    n   = 0;
    for (fID=fIDN-1; fID >= fID0; fID--) {
	fm0 = frameIO(p, fio, fID);		// load the *this* frame
	if (p->alg == 0)
	    frameSpots(p, fm0, fm1);
	else
	    frameSpot2(p, fm0, fm1);
	frameDelete(fm1);
	frame_sum(p, fm0);
	fm1 = fm0;
	n ++;
	if (oneprec > 0 && n % oneprec == 0) {
	    fprintf(stderr, ".");
	    fflush(stderr);
	}
    }
    frameIO_Close(p, fio);
    frameDelete(fm1);
}

/*-------------------------------------------------------------------------
 *
 *	Spot finding: from a single frame.
 *
 *------------------------------------------------------------------------*/

void spot_sframe(para_t *p)
{
    frameIO_t  *fio;
    frameloc_t *fm=NULL;

    fio = frameIO_Open(p);
    fm  = frameIO(p, fio, p->frameID2);
    if (p->alg == 0)
	frameSpots(p, fm, NULL);
    else
	frameSpot2(p, fm, NULL);
    frame_sum(p, fm);
    frameDelete(fm);
}

/*-------------------------------------------------------------------------
 *
 *	Spot handling: fitting for all the candidate spots.
 *
 *------------------------------------------------------------------------*/

void spot_fitting(para_t *p)
{
    FILE    *f0, *f1;
    double  *x_fit, *y_fit, dt;
    sp_t   **sp;
    int      i, x, y, r, n, oneperc; 
    int      sdim_x, sdim_y, x_rng, y_rng, imglen;

//  Prepare to show the progress of running.
    printf("Found candidate spots: %d\n", p->n_sp1);
    fflush(stdout);
    oneperc = p->n_sp1 / 100;

//  Prepare coordinate of pixels of the spot, which is relative to its center.
    sdim_x = p->x_find_pixels;
    sdim_y = p->y_find_pixels;
    imglen = sdim_x * sdim_y;
    x_rng  = sdim_x / 2;
    y_rng  = sdim_y / 2;
    x_fit  = malloc(imglen * sizeof(double));
    y_fit  = malloc(imglen * sizeof(double));
    if (! x_fit || ! y_fit)
        pstop("!!! spot_fitting: not enough memory.\n");
    i = 0;
    for (y=-y_rng; y <= y_rng; y++) {
    for (x=-x_rng; x <= x_rng; x++) {
        x_fit[i] = x;
        y_fit[i] = y;
        i ++;
    }}

//  Fitting for the normal spots.
    sp = p->sp1;
    n  = 0;
#pragma omp parallel for private(i,r) reduction(+:n)
    for (i=0; i < p->n_sp1; i++) {
	r = SpotFit(p, x_fit, y_fit, sp[i]);
	if (r == 0) n ++;
	if (oneperc > 0 && i % oneperc == 0) {
	    fprintf(stderr, ".");
	    fflush(stderr);
	}
    }
    printf("\nTotal valid particles: %d\n", n);

//  Fitting for the high-intensity spots.
    sp = p->sp2;
    n  = 0;
#pragma omp parallel for private(i,r) reduction(+:n)
    for (i=0; i < p->n_sp2; i++) {
	r = SpotFit(p, x_fit, y_fit, sp[i]);
	if (r == 0) n++;
    }
    printf("Total high-intensity particles: %d\n", n);

//  Output the results of normal spots, and update the statistics.
    f0 = out_fit_init(p, p->outfn, 0);
    f1 = out_fit_init(p, p->outfn, 1);
    sp = p->sp1;
    for (i=0, n=0; i < p->n_sp1; i++) {
	if (sp[i] && sp[i]->res) {
	    out_fit(p, f0, n, sp[i]);
	    out_fit(p, f1, n, sp[i]);
	    x = sp[i]->fID - p->frameID1;
	    p->fsts[x].n_event ++;
	    n ++;
	}
    }
    dt = get_realtime();
    fprintf(f1, "ExecTime: %E sec\n", dt);
    printf("ExecTime: %E sec\n", dt);
    out_fit_close(f0);
    out_fit_close(f1);

//  Output the results of high-intensity spots.
    f0 = out_fit_init(p, p->outfnH, 0);
    f1 = out_fit_init(p, p->outfnH, 1);
    sp = p->sp2;
    for (i=0, n=0; i < p->n_sp2; i++) {
	if (sp[i] && sp[i]->res) {
	    out_fit(p, f0, n, sp[i]);
	    out_fit(p, f1, n, sp[i]);
	    n ++;
	}
    }
    out_fit_close(f0);
    out_fit_close(f1);

//  Output the statistics and the sum of pixels.
    out_framests(p);
    out_framesum(p);
}

/*-------------------------------------------------------------------------
 *
 *	Spot handling: output the images of all candidate spots.
 *
 *------------------------------------------------------------------------*/

void spot_output_img(para_t *p)
{
    FILE  *f;
    float *data;
    int    i, j, n, imglen, *img;
    const char *fn;

    fn     = "spot.img";
    imglen = p->x_find_pixels * p->y_find_pixels;
    if ((f = fopen(fn, "wb")) == NULL)
       	pstop("!!! spot_output_img: cannot open file: %s\n", fn);
    if ((data = malloc(imglen*sizeof(float))) == NULL)
        pstop("!!! spot_output_img: not enough memory.\n");

    n = 0;
    for (i=0; i < p->n_sp1; i++) {
	if (p->sp1[i] == NULL || p->sp1[i]->img == NULL) continue;
	img = p->sp1[i]->img;

	for (j=0; j < imglen; j++)
	    data[j] = (float)img[j];
	fwrite(data, sizeof(float), imglen, f);
	n ++;
    }
    printf("Number of found spots: %d\n", n);
    fclose(f);
    free(data);
}

