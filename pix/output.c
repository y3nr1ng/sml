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

FILE *out_fit_init(para_t *p, char *outfn)
{
    FILE *f;
    char  buf[1024];

    if (outfn == NULL) return NULL;

    if (p->nprc > 1)
	sprintf(buf, "%s.%04d", outfn, p->myid);
    else
	sprintf(buf, "%s.txt",  outfn);
    if ((f = fopen(buf, "wt")) == NULL)
        pstop("!!! inputs: cannot open output file (pixel): %s\n", buf);

    fprintf(f, "%-8s %-6s  %4s %4s %4s  %13s %10s %13s %10s %13s %10s",
	       "spot", "frame", "x(p)", "y(p)", "cnt",
	       "Intensity", "dI", "x", "dx", "y", "dy");
    if (p->mode == 0)
	fprintf(f, " %13s %10s %13s %10s %13s %11s\n",
		   "w", "dw", "Background", "dB", "S/N", "chisq/ndof");
    else {
	fprintf(f, " %13s %10s %13s %10s %13s %10s %13s %11s %13s %10s",
		   "wx", "dwx", "wy", "dwy", "Background", "dB", "S/N",
		   "chisq/ndof", "z(wx/wy)", "dz(wx/wy)");
	if (p->verb)
	    fprintf(f, " %13s %10s %13s %10s",
		       "z(wx)", "dz(wx)", "z(wy)", "dz(wy)"); 
	fprintf(f, "\n");
    }

    fprintf(f, "%s%s%s%s",
	       "----------------------------------------------------",
	       "----------------------------------------------------",
	       "----------------------------------------------------",
	       "---------------------------");
    if (p->mode == 0)
	fprintf(f, "\n");
    else {
	fprintf(f, "---------------------------------------------------");
	if (p->verb)
	    fprintf(f, "--------------------------------------------------");
	fprintf(f, "\n");
    }
    return f;
}

void out_fit(para_t *p, FILE *f, int Sid, sp_t *sp)
{
    int     r;
    double *a, *da, chisq;
    double  an[10], dan[10];
    char    r6[30], r7[30], r8[30];
    char    e6[30], e7[30], e8[30];

    if (f == NULL) return;

    a     = sp->res->a;
    da    = sp->res->da;
    chisq = sp->res->chisq;
    if (p->mode == 0) {				// 2D fit
	an[0]  = a[0] * p->i_photon;		// I0:  intensity
	an[1]  = a[1] * p->nm_px_x;		// x0:  x position
	an[2]  = a[2] * p->nm_px_y;		// y0:  y position
	an[3]  = a[3] * p->nm_px_x;		// w0:  width
	an[4]  = a[4] * p->i_photon;		// B0:  background
	dan[0] = da[0] * p->i_photon;
	dan[1] = da[1] * p->nm_px_x;
	dan[2] = da[2] * p->nm_px_y;
	dan[3] = da[3] * p->nm_px_x;
	dan[4] = da[4] * p->i_photon;
    }
    else {					// 3D fit
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

	r = solve_z_wxowy(p, an, dan, r6, e6);	// z(wx/wy)
	if (p->verb == 0 && r != 0) return;
    }

    fprintf(f, "%08d %06d  %4d %4d %4d", Sid, sp->fID, sp->x, sp->y, sp->cnt);
    fprintf(f, "  %13.6E %10.3E %13.6E %10.3E %13.6E %10.3E",
		an[0], dan[0], an[1], dan[1], an[2], dan[2]);

    if (p->mode == 0) {
	fprintf(f, " %13.6E %10.3E %13.6E %10.3E %13.6E  %10.3E\n",
		    an[3], dan[3], an[4], dan[4], an[0]/an[4], chisq);
    }
    else if (p->mode == 1) {
	fprintf(f, " %13.6E %10.3E %13.6E %10.3E %13.6E %10.3E %13.6E",
	            an[3], dan[3], an[4], dan[4], an[5], dan[5], an[0]/an[5]);
	fprintf(f, "  %10.3E %s %s", chisq, r6, e6);

	if (p->verb >= 1) {
	    solve_z_w(p, &(p->cax), an[3], dan[3], r7, e7);
	    solve_z_w(p, &(p->cay), an[4], dan[4], r8, e8);
	    fprintf(f, " %s %s %s %s", r7, e7, r8, e8);
	}
	fprintf(f, "\n");
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

    if (p->nprc > 1)
	sprintf(fn, "%s.%04d", p->outfnF, p->myid);
    else
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
    if (p->nprc > 1)
	sprintf(fn, "%s.%04d", p->outfnS, p->myid);
    else
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
