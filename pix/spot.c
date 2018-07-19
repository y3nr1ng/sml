#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_OMP
#include <omp.h>
#endif
#include "pix.h"

/*---------------------------------------------------------------------------
 *
 *      Distribute frames to each threads or parallel processes
 *
 *--------------------------------------------------------------------------*/

static void
frameDistribution(int myid, int nprc, int tot_fID0, int tot_fIDN,
		  int *fID0, int *fIDN)
{
    int  nframe, nhandle, ID0, IDN;

    nframe  =  tot_fIDN - tot_fID0 + 1;
    nhandle =  nframe / nprc;
    nhandle = (nframe % nprc == 0) ? nhandle : nhandle+1;
    ID0     =  myid * nhandle + tot_fID0;
    IDN     = (ID0+nhandle <= tot_fIDN) ? ID0+nhandle : tot_fIDN;
    *fID0   =  ID0;
    *fIDN   =  IDN;
}

/*-------------------------------------------------------------------------
 *
 *	Initialization of spot finding for each thread.
 *
 *------------------------------------------------------------------------*/

static parath_t *frame_spot_thinit(para_t *p, int myid, int nprc)
{
    parath_t *pp;

    if ((pp = calloc(sizeof(parath_t), 1)) == NULL)
	pstop("!!! frame_spot_thinit: not enough memory.\n");
    pp->myid  = myid;
    pp->nprc  = nprc;
    pp->nfsep = p->nfsep;
    pp->x_find_pixels = p->x_find_pixels;
    pp->y_find_pixels = p->y_find_pixels;

    return pp;
}

/*-------------------------------------------------------------------------
 *
 *	Merge the found candidate spots for each thread.
 *
 *------------------------------------------------------------------------*/

static void merge_spots(para_t *p, parath_t *pp, int mode)
{
    int    nsp, nspp;
    sp_t **sp, **spp;

    if (mode == 0) {
	nsp  = p->n_sp1;
	nspp = pp->n_sp1;
	sp   = p->sp1;
	spp  = pp->sp1;
    }
    else {
	nsp  = p->n_sp2;
	nspp = pp->n_sp2;
	sp   = p->sp2;
	spp  = pp->sp2;
    }

    if (nsp == 0) {
	nsp = nspp;
	sp  = spp;
    }
    else {
	if ((sp = realloc(sp, sizeof(sp_t *)*(nsp+nspp))) == NULL)
	    pstop("!!! merge_spots: not enough memory for sp.\n");
	memcpy(sp+nsp, spp, sizeof(sp_t *)*nspp);
	nsp += nspp;
	free(spp);
    }

    if (mode == 0) {
	p->n_sp1 = nsp;
	p->sp1   = sp;
    }
    else {
	p->n_sp2 = nsp;
	p->sp2   = sp;
    }
}

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

void dframe_spot(para_t *p)
{
    frameIO_t  *fio;
    frameloc_t *fm0=NULL, *fm1=NULL;
    int         myid, nprc, fID, fID0, fIDN, oneprec, n;
    parath_t   *pp;

#ifdef USE_OMP
    myid = omp_get_thread_num();
    nprc = omp_get_num_threads();
#else
    myid = p->myid;
    nprc = p->nprc;
#endif
    frameDistribution(myid, nprc, p->frameID1, p->frameID2, &fID0, &fIDN);
    pp  = frame_spot_thinit(p, myid, nprc);
    fio = frameIO_Open(p);
    oneprec = (fIDN-fID0+1) / (100 / nprc);
 
    fm1 = frameIO(p, fio, fIDN);		// load the *next* frame
    n   = 0;
    for (fID=fIDN-1; fID >= fID0; fID--) {
	fm0 = frameIO(p, fio, fID);		// load the *this* frame
	if (p->alg == 0)
	    frameSpots(p, pp, fm0, fm1);
	else
	    frameSpot2(p, pp, fm0, fm1);
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
#ifdef USE_OMP
#pragma omp barrier
#endif

#ifdef USE_OMP
#pragma omp critical
#endif
    {
	merge_spots(p, pp, 0);
	merge_spots(p, pp, 1);
    }
}

/*-------------------------------------------------------------------------
 *
 *	Spot finding: from a single frame.
 *
 *------------------------------------------------------------------------*/

void sframe_spot(para_t *p)
{
    parath_t   *pp;
    frameIO_t  *fio;
    frameloc_t *fm=NULL;

    pp  = frame_spot_thinit(p, p->myid, p->nprc);
    fio = frameIO_Open(p);
    fm  = frameIO(p, fio, p->frameID2);
    if (p->alg == 0)
	frameSpots(p, pp, fm, NULL);
    else
	frameSpot2(p, pp, fm, NULL);
    frame_sum(p, fm);
    frameDelete(fm);

    p->n_sp1 = pp->n_sp1;
    p->n_sp2 = pp->n_sp2;
    p->sp1   = pp->sp1;
    p->sp2   = pp->sp2;
}

/*-------------------------------------------------------------------------
 *
 *	Spot handling: fitting for all the candidate spots.
 *
 *------------------------------------------------------------------------*/

void spot_fitting(para_t *p)
{
    FILE    *f;
    double  *x_fit, *y_fit, dt;
    sp_t   **sp;
    int      i, x, y, r, n, oneperc; 
    int      sdim_x, sdim_y, x_rng, y_rng, imglen;

//  Prepare to show the progress of running.
    printf("%d: found candidate spots: %d\n", p->myid, p->n_sp1);
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
        pstop("!!! SpotFit: not enough memory.\n");
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
#ifdef USE_OMP
#pragma omp parallel for private(i,r) reduction(+:n)
#endif
    for (i=0; i < p->n_sp1; i++) {
	r = SpotFit(p, x_fit, y_fit, sp[i]);
	if (r == 0) n ++;
	if (oneperc > 0 && i % oneperc == 0) {
	    fprintf(stderr, ".");
	    fflush(stderr);
	}
    }
    printf("\n%d: total valid particles: %d\n", p->myid, n);

//  Fitting for the high-intensity spots.
    sp = p->sp2;
    n  = 0;
#ifdef USE_OMP
#pragma omp parallel for private(i,r) reduction(+:n)
#endif
    for (i=0; i < p->n_sp2; i++) {
	r = SpotFit(p, x_fit, y_fit, sp[i]);
	if (r == 0) n++;
    }
    printf("%d: total high-intensity particles: %d\n", p->myid, n);

//  Output the results of normal spots, and update the statistics.
    f  = out_fit_init(p, p->outfn);
    sp = p->sp1;
    for (i=0, n=0; i < p->n_sp1; i++) {
	if (sp[i] && sp[i]->res) {
	    out_fit(p, f, n, sp[i]);
	    x = sp[i]->fID - p->frameID1;
	    p->fsts[x].n_event ++;
	    n ++;
	}
    }
    if (p->myid == 0) {
	dt = get_realtime();
	fprintf(f, "ExecTime: %E sec\n", dt);
	printf("ExecTime: %E sec\n", dt);
    }
    out_fit_close(f);

//  Output the results of high-intensity spots.
    f  = out_fit_init(p, p->outfnH);
    sp = p->sp2;
    for (i=0, n=0; i < p->n_sp2; i++) {
	if (sp[i] && sp[i]->res) {
	    out_fit(p, f, n, sp[i]);
	    n ++;
	}
    }
    out_fit_close(f);

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

