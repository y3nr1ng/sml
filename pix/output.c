#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pix.h"

/*-------------------------------------------------------------------------
 *
 *  Output the list of found spots (candidates)
 *
 *------------------------------------------------------------------------*/

void out_spotlist(char *outfn, int nsp, int sqsz, sp_t **sp)
{
    FILE *f;
    char  buf[1024];
    int   i, ii, sID, fID0;

    if (outfn == NULL) return;

    sprintf(buf, "%s.txt",  outfn);
    if ((f = fopen(buf, "wt")) == NULL)
        pstop("!!! inputs: cannot open output file (pixel): %s\n", buf);

    fprintf(f, "%-8s %-6s %4s %4s %8s %6s\n",
	       "frame", "spot", "x", "y", "II", "cnt");
    fprintf(f, "-----------------------------------------\n");
    fID0 = -1;
    sID  =  0;
    for (i=0; i < nsp; i++) {
	if (sp[i] == NULL || sp[i]->img == NULL) continue;

	if (fID0 != sp[i]->fID) {
	    fID0 = sp[i]->fID;
	    sID  = 0;
	}
	ii = sqsz*sqsz/2;
	fprintf(f, "%08d %06d %04d %04d %8d %6d\n",
		sp[i]->fID, sID, sp[i]->x, sp[i]->y,
		sp[i]->img[ii]/sp[i]->cnt, sp[i]->cnt);
	sID ++;
    }
    fclose(f);
}

/*-------------------------------------------------------------------------
 *
 *  Output the fitting results.
 *
 *------------------------------------------------------------------------*/

static void out_z(int n_z, double *z, double *dz, char *z1, char *z2,
		  char *dz1, char *dz2)
{
    if (n_z > 0) {
	sprintf(z1,  "%13.6E", z[0]);
	sprintf(dz1, "%10.3E", dz[0]);
    }
    else {
	sprintf(z1,  "%13.6E", 0.0);
	sprintf(dz1, "%10.3E", 0.0);
    }
    if (z2 != NULL && dz2 != NULL) {
	if (n_z > 1 && z2 != NULL && dz2 != NULL) {
	    sprintf(z2,  "%13.6E", z[1]);
	    sprintf(dz2, "%10.3E", dz[1]);
	}
	else {
	    sprintf(z2,  "%13.6E", 0.0);
	    sprintf(dz2, "%10.3E", 0.0);
	}
    }
}

FILE *out_fit_init(para_t *p, char *outfn, int head)
{
    FILE *f;
    char  buf[1024];

    if (outfn == NULL) return NULL;

    if (head == 0)
	sprintf(buf, "%s0.txt",  outfn);
    else
	sprintf(buf, "%s.txt",  outfn);
    if ((f = fopen(buf, "wt")) == NULL)
        pstop("!!! inputs: cannot open output file (pixel): %s\n", buf);
    if (head == 0) return f;

    fprintf(f, "%-8s %-6s  %4s %4s %4s  %13s %10s %13s %10s %13s %10s",
	       "spot", "frame", "x(p)", "y(p)", "cnt",
	       "Intensity", "dI", "x", "dx", "y", "dy");
    if (p->mode == 0)
	fprintf(f, " %13s %10s %13s %10s %13s %11s\n",
		   "w", "dw", "Background", "dB", "S/N", "chisq/ndof");
    else {
	fprintf(f, " %13s %10s %13s %10s %13s %10s %13s %11s", "wx", "dwx",
		   "wy", "dwy", "Background", "dB", "S/N", "chisq/ndof");
	fprintf(f, " %13s %10s %13s %10s %13s %10s %13s %10s %13s %10s\n",
		   "z", "dz", "z1(wx)", "dz1(wx)", "z2(wx)", "dz2(wx)",
		   "z1(wy)", "dz1(wy)", "z2(wy)", "dz2(wy)"); 
    }

    fprintf(f, "%s%s%s%s",
	       "----------------------------------------------------",
	       "----------------------------------------------------",
	       "----------------------------------------------------",
	       "---------------------------");
    if (p->mode == 0)
	fprintf(f, "\n");
    else {
	fprintf(f, "%s%s%s\n",
		   "---------------------------------------------------",
	    	   "---------------------------------------------------",
		   "-------------------------------------------------");
    }
    return f;
}

void out_fit(para_t *p, FILE *f, int Sid, sp_t *sp)
{
    fitres_t *res;
    double   *a, *da, chisq;
    double    an[10], dan[10];

    if (f == NULL) return;

    res   = sp->res;
    a     = res->a;
    da    = res->da;
    chisq = res->chisq;
    if (p->mode == 0) {				// 2D fit
	an[0]  = a[0] * p->i_photon;		// I0:  intensity
	an[1]  = a[1] * p->nm_px_x;		// x0:  x position
	an[2]  = a[2] * p->nm_px_y;		// y0:  y position
	an[3]  = a[3] * p->nm_px_x;		// w0:  width
	an[4]  = a[4] * p->i_photon;		// B0:  background
	an[5]  = 0.0;
	dan[0] = da[0] * p->i_photon;
	dan[1] = da[1] * p->nm_px_x;
	dan[2] = da[2] * p->nm_px_y;
	dan[3] = da[3] * p->nm_px_x;
	dan[4] = da[4] * p->i_photon;
	dan[5] = 0.0;
    }
    else {					// 3D fit
//	if (res->n_zx == 0 && res->n_zy == 0 && res->n_zr == 0) return;
	if (res->z == 0.0 && res->dz == 0.0) return;

	an[0]  = a[0] * p->i_photon;		// I0:  intesity
	an[1]  = a[1] * p->nm_px_x;		// x0:  x position
	an[2]  = a[2] * p->nm_px_y;		// y0:  y position
	an[3]  = a[3] * p->nm_px_x;		// wx:  x width
	an[4]  = a[4] * p->nm_px_y;		// wy:  y width
	an[5]  = a[5] * p->i_photon;		// B0:  background
	dan[0] = da[0] * p->i_photon;
	dan[1] = da[1] * p->nm_px_x;
	dan[2] = da[2] * p->nm_px_y;
	dan[3] = da[3] * p->nm_px_x;
	dan[4] = da[4] * p->nm_px_y;
	dan[5] = da[5] * p->i_photon;
    }

    fprintf(f, "%08d %06d  %4d %4d %4d", Sid, sp->fID, sp->x, sp->y, sp->cnt);
    fprintf(f, "  %13.6E %10.3E %13.6E %10.3E %13.6E %10.3E",
		an[0], dan[0], an[1], dan[1], an[2], dan[2]);

    if (p->mode == 0) {
	fprintf(f, " %13.6E %10.3E %13.6E %10.3E %13.6E  %10.3E\n",
		    an[3], dan[3], an[4], dan[4], an[0]/an[4], chisq);
    }
    else if (p->mode == 1) {
	char  buf[1024];
	char *zx1, *zx2, *zy1, *zy2, *z;
	char *dx1, *dx2, *dy1, *dy2, *dz;

	zx1 = buf;
	zx2 = zx1 + 30;
	zy1 = zx2 + 30;
	zy2 = zy1 + 30;
	dx1 = zy2 + 30;
	dx2 = dx1 + 30;
	dy1 = dx2 + 30;
	dy2 = dy1 + 30;
	z   = dy2 + 30;
	dz  = z   + 30;
	out_z(res->n_zx, res->zx, res->dzx, zx1, zx2, dx1, dx2);
	out_z(res->n_zy, res->zy, res->dzy, zy1, zy2, dy1, dy2);
	out_z(1, &(res->z), &(res->dz), z, NULL, dz, NULL);

	fprintf(f, " %13.6E %10.3E %13.6E %10.3E %13.6E %10.3E %13.6E  %10.3E",
	       an[3], dan[3], an[4], dan[4], an[5], dan[5], an[0]/an[5], chisq);
	fprintf(f, " %s %s %s %s %s %s %s %s %s %s\n",
	       z, dz, zx1, dx1, zx2, dx2, zy1, dy1, zy2, dy2);
    }
}

void out_fit_close(FILE *f)
{
    if (f == NULL) return;
    fclose(f);
}

/*-------------------------------------------------------------------------
 *
 *  Output the frame statistics.
 *
 *------------------------------------------------------------------------*/

void out_framests(para_t *p)
{
    FILE  *f;
    char   fn[1024];
    int    i, ntotal, n_evt, n_pxl;
    double density;

    if (p->outfnF == NULL) return;

    sprintf(fn, "%s.txt",  p->outfnF);
    if ((f = fopen(fn, "wt")) == NULL)
	pstop("!!! out_framests: cannot open file: %s\n", fn);
    fprintf(f, "%-8s  %8s  %8s\n", "frame", "n_event", "density");
    fprintf(f, "----------------------------\n");

    ntotal = p->frameID2 - p->frameID1;
    n_pxl  = (p->frame_x2 - p->frame_x1 + 1) * (p->frame_y2 - p->frame_y1 + 1);
    for (i=0; i < ntotal; i++) {
	n_evt   = p->fsts[i].n_event;
	density = ((double)(n_evt * p->x_find_pixels * p->y_find_pixels))
		/ ((double)n_pxl);
	fprintf(f, "%08d  %8d  %.6f\n", i+p->frameID1, n_evt, density);
    }
    fclose(f);
}

/*-------------------------------------------------------------------------
 *
 *      Output the pixel values of the summed frames.
 *
 *------------------------------------------------------------------------*/

void out_framesum(para_t *p)
{
    FILE *f;
    int   i, j, k, dim_x, dim_y;
    char  fn[1024];

    if (p->outfnS == NULL) return;

    dim_x = p->frame_x2 - p->frame_x1 + 1;
    dim_y = p->frame_y2 - p->frame_y1 + 1;
    sprintf(fn, "%s.txt",  p->outfnS);
    if ((f = fopen(fn, "wt")) == NULL)
        pstop("!!! out_imgtxt: cannot open file: %s\n", fn);
    fprintf(f, "%4s  %4s  %8s\n", "x", "y", "pixel");
    fprintf(f, "--------------------\n");

    k = 0;
    for (j=0; j < dim_y; j++) {
    for (i=0; i < dim_x; i++) {
        fprintf(f, "%04d  %04d  %8d\n", i, j, p->psum[k++]);
    }}
    fclose(f);
}
