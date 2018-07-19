#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pclst.h"

/*-------------------------------------------------------------------------
 *
 *	Output the cluster result.
 *
 *------------------------------------------------------------------------*/

void output(para_t *p)
{
    long    i, j, k=0;
    rec_t  *rec;
    clst_t *cl = p->cl;

    for (i=0; i < p->ncl; i++, cl++) {
	if (cl->np == 0) continue;

	printf("cluster %ld: %ld\n", k, cl->np);
	printf("border: %14.6E %14.6E, %14.6E %14.6E\n",
		cl->x1[0], cl->x1[1], cl->x2[0], cl->x2[1]);
	for (j=0; j < p->cl[i].np; j++)
	    printf("%06ld ", cl->sID[j]);
	printf("\n\n");

	if (cl->np >= p->nsp_draw && p->n_rec < p->nrec_draw) {
	    rec = push_rec(p, cl->x1[0], cl->x2[0], cl->x1[1], cl->x2[1]);
	    xcor_spots(p, cl, rec);
	}
	k ++;
    }
}

/*-------------------------------------------------------------------------
 *
 *	Cluster management.
 *
 *------------------------------------------------------------------------*/

static void push_spot_cluster(clst_t *cl, long sID)
{
    if (cl->np >= cl->mp) {
	cl->mp += 16;
	if ((cl->sID = realloc(cl->sID, cl->mp*sizeof(long))) == NULL)
	    pstop("push_spot_cluster: not enough memory.\n");
    }
    cl->sID[cl->np] = sID;
    cl->np ++;
}

static long new_cluster(para_t *p, long sID1, long sID2)
{
    long cID;

    if (p->ncl >= p->mcl) {
	p->mcl += 16;
	if ((p->cl = realloc(p->cl, p->mcl*sizeof(clst_t))) == NULL)
	    pstop("new_cluster: not enough memory.\n");
    }
    cID = p->ncl;
    p->cl[cID].np  = 0;
    p->cl[cID].mp  = 0;
    p->cl[cID].sID = NULL;
    push_spot_cluster(p->cl+cID, sID1);
    push_spot_cluster(p->cl+cID, sID2);
    p->sp[sID1].cID = cID;
    p->sp[sID2].cID = cID;
    p->ncl ++;

    return cID;
}

static void cluster_add_spot(para_t *p, long cID, long sID)
{
    push_spot_cluster(p->cl+cID, sID);
    p->sp[sID].cID = cID;
}

/*-------------------------------------------------------------------------
 *
 *	Identify the boundary of clusters.
 *
 *------------------------------------------------------------------------*/

static void clst_boundary(clst_t *cl, spot_t *sp)
{
    long   i, sID;
    double x, y, x1[2], x2[2];

    if (cl->np == 0) return;
    x1[0] = 1.E+100;
    x1[1] = 1.E+100;
    x2[0] = 0.0;
    x2[1] = 0.0;
    for (i=0; i < cl->np; i++) {
	sID = cl->sID[i];
	x   = sp[sID].x;
	y   = sp[sID].y;
	x1[0] = (x < x1[0]) ? x : x1[0];
	x1[1] = (y < x1[1]) ? y : x1[1];
	x2[0] = (x > x2[0]) ? x : x2[0];
	x2[1] = (y > x2[1]) ? y : x2[1];
    }
    cl->x1[0] = x1[0];
    cl->x1[1] = x1[1];
    cl->x2[0] = x2[0];
    cl->x2[1] = x2[1];
}

/*-------------------------------------------------------------------------
 *
 *	Identify clusters among the spots.
 *
 *	First stage: identify clusters based on half range.
 *
 *------------------------------------------------------------------------*/

static void cluster_identify_stage1(para_t *p)
{
    long    i, j, cID, np;
    double  x1, x2, y1, y2, xc, yc, dx, dy, dd, spd_h;

    spd_h = p->spd / 4.0;
    for (i=0; i < p->np; i++) {
	if (p->sp[i].cID != -1) continue;

	x1  = p->sp[i].x;
	y1  = p->sp[i].y;
        xc  = x1;
        yc  = y1;
	cID = -1;
	for (j=i+1; j < p->np; j++) {
	    if (p->sp[j].cID != -1) continue;

	    x2 = p->sp[j].x;
	    y2 = p->sp[j].y;
	    dx = x2-x1;
	    dy = y2-y1;
	    dd = dx*dx + dy*dy;
	    if (dd > spd_h) continue;

	    if (cID == -1)
		cID = new_cluster(p, i, j);
	    else
		cluster_add_spot(p, cID, j);
	    xc += x2;
	    yc += y2;
	}
	if (cID != -1) {
	    np = p->cl[cID].np;
	    p->cl[cID].xc[0] = xc / ((double)np);
	    p->cl[cID].xc[1] = yc / ((double)np);
	}
    }
}

/*-------------------------------------------------------------------------
 *
 *	Clean out the small clusters (only keep the largest ncls1 clusters)
 *
 *------------------------------------------------------------------------*/

static int cmp_cl(const void *a, const void *b)
{
    clst_t  *cl1, *cl2;

    cl1 = (clst_t *)a;
    cl2 = (clst_t *)b;
    if (cl1->np < cl2->np)
	return 1;
    else if (cl1->np > cl2->np)
	return -1;
    else
	return 0;
}

static void cluster_clean(para_t *p)
{
    long    i, j, sID, ncl;
    clst_t *cl;

    qsort(p->cl, p->ncl, sizeof(clst_t), cmp_cl);

    ncl = p->ncl;
    for (i=p->ncls1; i < ncl; i++) {
	cl = p->cl+i;
	for (j=0; j < cl->np; j++) {
	    sID = cl->sID[j];
	    p->sp[sID].cID = -1;
	}
	free(cl->sID);
	cl->np  = 0;
	cl->mp  = 0;
	cl->sID = NULL;
	p->ncl --;
    }
}

/*-------------------------------------------------------------------------
 *
 *	Identify clusters among the spots.
 *
 *	Second stage: identify far-away spots into existing clusters.
 *
 *------------------------------------------------------------------------*/

static void cluster_identify_stage2(para_t *p)
{
    long    i, j, cmin;
    double  x1, x2, y1, y2, dx, dy, dd, dmin=0.0;

    for (i=0; i < p->np; i++) {
	if (p->sp[i].cID != -1) continue;

	x1 = p->sp[i].x;
	y1 = p->sp[i].y;
	cmin = -1;
	for (j=0; j < p->ncl; j++) {
	    x2 = p->cl[j].xc[0];
	    y2 = p->cl[j].xc[1];
	    dx = x2-x1;
	    dy = y2-y1;
	    dd = dx*dx + dy*dy;
	    if (dd > p->spd) continue;

	    if (cmin == -1 || dd < dmin) {
		dmin = dd;
		cmin = j;
	    }
	}
	if (cmin != -1)
	    cluster_add_spot(p, cmin, i);
    }
}

/*-------------------------------------------------------------------------
 *
 *	Identify clusters among the spots: main routine.
 *
 *------------------------------------------------------------------------*/

void cluster_identify(para_t *p)
{
    long  i;

    printf("Identifying clusters for %ld particles ....\n", p->np);
    cluster_identify_stage1(p);
    cluster_clean(p);
    cluster_identify_stage2(p);
    for (i=0; i < p->ncl; i++)
	clst_boundary(p->cl+i, p->sp);
    qsort(p->cl, p->ncl, sizeof(clst_t), cmp_cl);
}
