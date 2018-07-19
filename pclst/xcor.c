/*--------------------------------------------------------------------------
 *
 *	Use OpenGL to display the image on the screen
 *
 *-------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pclst.h"

/*--------------------------------------------------------------------------
 *
 *	Record the rectangles.
 *
 *-------------------------------------------------------------------------*/

rec_t *push_rec(para_t *p, double x1, double x2, double y1, double y2)
{
    double  xx, tol=1.E-10;
    rec_t  *rec;

    if (p->n_rec >= p->m_rec) {
	p->m_rec += 10;
	if ((p->rec = realloc(p->rec, p->m_rec*sizeof(rec_t))) == NULL)
	    pstop("!!! push_rec: not enough memory.\n");
    }
    if (fabs(x1-x2) < tol && fabs(y1-y2) < tol)
	return NULL;
    if (x1 > x2) {
	xx = x1;
	x1 = x2;
	x2 = xx;
    }
    if (y1 > y2) {
	xx = y1;
	y1 = y2;
	y2 = xx;
    }
    p->rec[p->n_rec].x1 = x1;
    p->rec[p->n_rec].x2 = x2;
    p->rec[p->n_rec].y1 = y1;
    p->rec[p->n_rec].y2 = y2;
    p->rec[p->n_rec].hist = NULL;
    rec = p->rec + p->n_rec;
    p->n_rec ++;

    return rec;
}

/*--------------------------------------------------------------------------
 *
 *	Find all the spots in each rectangle.
 *
 *-------------------------------------------------------------------------*/

void xcor_spots(para_t *p, clst_t *cl, rec_t *r)
{
    long    i, j, k, sid1, sid2;
    double  maxr, dr, rr, x1, x2, y1, y2;

    if (r == NULL) return;
    x1   = r->x1;
    x2   = r->x2;
    y1   = r->y1;
    y2   = r->y2;
    maxr = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
    dr   = maxr / ((double)(p->n_intvl));
    if ((r->hist = calloc(p->n_intvl, sizeof(int))) == NULL)
	pstop("!!! xcor_spots: not enough memory.\n");
    r->nsp = cl->np;
    r->sp  = cl->sID;

    for (j=0; j < cl->np; j++) {
	sid2 = cl->sID[j];
	x2   = p->sp[sid2].x;
	y2   = p->sp[sid2].y;
	for (i=0; i < cl->np; i++) {
	    sid1 = cl->sID[i];
	    x1   = p->sp[sid1].x;
	    y1   = p->sp[sid1].y;
	    rr   = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
	    k    = (int)(rr / dr);
	    if (k < p->n_intvl) r->hist[k] ++;
	}
    }
}

/*--------------------------------------------------------------------------
 *
 *	Data I/O for spot correlation
 *
 *-------------------------------------------------------------------------*/

void xcor_output(para_t *p)
{
    char  *fn;
    FILE  *f;
    long   i, j;
    double x1, x2, y1, y2, maxr, dr, rr;

    fn = "xcor.dat";
    if ((f = fopen(fn, "wt")) == NULL)
	return;
    fprintf(f, "N_REC:  %ld\n\n", p->n_rec);
    for (i=0; i < p->n_rec; i++) {
	x1   = p->rec[i].x1;
	x2   = p->rec[i].x2;
	y1   = p->rec[i].y1;
	y2   = p->rec[i].y2;
	maxr = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
	dr   = maxr / ((double)(p->n_intvl));

	fprintf(f, "REC_ID: %ld\n", i);
	fprintf(f, "REC_C1: %.16E %.16E\n", x1, y1);
	fprintf(f, "REC_C2: %.16E %.16E\n", x2, y2);
	fprintf(f, "REC_NP: %ld\n", p->rec[i].nsp);
	fprintf(f, "REC_NH: %d\n", p->n_intvl);
	fprintf(f, "REC_F0: %d\n", p->Fstart);
	fprintf(f, "REC_DR: %.16E\n", dr);
	fprintf(f, "REC_SP:");
	for (j=0; j < p->rec[i].nsp; j++)
	    fprintf(f, " %08ld", p->rec[i].sp[j]);
	fprintf(f, "\n");

	for (j=0; j < p->n_intvl; j++) {
	    rr = j * dr;
	    fprintf(f, "%.16E  %d\n", rr, p->rec[i].hist[j]);
	}
	fprintf(f, "\n");
    }
    fclose(f);
}

