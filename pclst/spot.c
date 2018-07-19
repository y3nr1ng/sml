#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "pclst.h"

/*-------------------------------------------------------------------------
 *
 *      Read the frame statistics data file.
 *
 *------------------------------------------------------------------------*/

void readfsts(para_t *p)
{
    FILE   *f;
    char    s[1024];
    int     fID, n_evt;
    double  density;

    if ((f = fopen(p->Fstsfn, "rt")) == NULL)
        pstop("!!! readfsts: cannot open file: %s\n", p->Fstsfn);
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
 *	Read the spots data file.
 *
 *------------------------------------------------------------------------*/

static void push_spot(int sID, int fID, double x, double y, para_t *p)
{
    if (p->bx1 > x || p->bx2 < x || p->by1 > y || p->by2 < y) {
	if (p->verb > 0)
	    printf("spot %d in frame %d is out of the image.\n", sID, fID);
	return;
    }

    if (p->np >= p->mp) {
        p->mp += 1024;
        if ((p->sp = realloc(p->sp, p->mp*sizeof(spot_t))) == NULL)
            pstop("push_spot: not enough memory.\n");
    }
    p->sp[p->np].cID = -1;
    p->sp[p->np].x = x;
    p->sp[p->np].y = y;
    p->np ++;
}

void readspot(para_t *p)
{
    FILE   *f;
    char    s[1024];
    long    sID, fID;
    int     xp, yp, cnt;
    double  II, dII, x, y, dx, dy;

    printf("starting frame: %d\n", p->Fstart);
    if ((f = fopen(p->spotfn, "rt")) == NULL)
	pstop("!!! readspot: cannot open file: %s\n", p->spotfn);
    fgets(s, 1024, f);
    fgets(s, 1024, f);
    while (fgets(s, 1024, f) != NULL) {
	sscanf(s, "%ld %ld %d %d %d %lf %lf %lf %lf %lf %lf",
		   &sID, &fID, &xp, &yp, &cnt, &II, &dII, &x, &dx, &y, &dy);
	if (fID >= p->Fstart)
	    push_spot(sID, fID, x, y, p);
    }
    fclose(f);
}

