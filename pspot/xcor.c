/*--------------------------------------------------------------------------
 *
 *	Use OpenGL to display the image on the screen
 *
 *-------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pspot.h"

/*--------------------------------------------------------------------------
 *
 *	Record the rectangles.
 *
 *-------------------------------------------------------------------------*/

void push_rec(para_t *p, double x1, double x2, double y1, double y2)
{
    double xx;

    if (p->n_rec >= p->m_rec) {
	p->m_rec += 10;
	if ((p->rec = realloc(p->rec, p->m_rec*sizeof(rec_t))) == NULL)
	    pstop("!!! push_rec: not enough memory.\n");
    }
    if (x1 == x2 && y1 == y2)
	return;
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
    p->rec[p->n_rec].update  = 1;
    p->rec[p->n_rec].n_intvl = p->n_intvl;
    p->rec[p->n_rec].Fstart  = p->Fstart;
    p->rec[p->n_rec].n_sp    = 0;
    p->rec[p->n_rec].m_sp    = 0;
    p->rec[p->n_rec].sp      = NULL;
    p->rec[p->n_rec].hist    = NULL;
    p->n_rec ++;
}

void del_rec(para_t *p, double x, double y)
{
    int i, j;

    for (i=0; i < p->n_rec; i++) {
	if (x >= p->rec[i].x1 && x <= p->rec[i].x2 &&
	    y >= p->rec[i].y1 && y <= p->rec[i].y2) break;
    }
    if (i >= p->n_rec) return;

    for (j=i+1; j < p->n_rec; j++)
	memcpy(p->rec+j-1, p->rec+j, sizeof(rec_t));
    p->n_rec --;
}

rec_t *search_rec(para_t *p, double x, double y)
{
    int i;

    for (i=0; i < p->n_rec; i++) {
	if (x >= p->rec[i].x1 && x <= p->rec[i].x2 &&
	    y >= p->rec[i].y1 && y <= p->rec[i].y2) break;
    }
    return (i >= p->n_rec) ? NULL : p->rec+i;
}

/*--------------------------------------------------------------------------
 *
 *	Find all the spots in each rectangle.
 *
 *-------------------------------------------------------------------------*/

static void find_spots(para_t *p, rec_t *r)
{
    int    i;
    double x1, x2, y1, y2;

    x1 = r->x1;
    y1 = r->y1;
    x2 = r->x2;
    y2 = r->y2;
    for (i=0; i < p->np; i++) {
	if (p->x[i] >= x1 && p->x[i] <= x2 &&
	    p->y[i] >= y1 && p->y[i] <= y2) {
	    if (r->n_sp >= r->m_sp) {
		r->m_sp += 256;
		if (! (r->sp = realloc(r->sp, r->m_sp * sizeof(int))))
		    pstop("!!! find_spots: not enough memory.\n");
	    }
	    r->sp[r->n_sp] = i;
	    r->n_sp ++;
	}
    }
}

static void xcor_spots(para_t *p, rec_t *r)
{
    int    i, j, k, jdx, idx;
    double maxr, dr, rr, x1, x2, y1, y2;

    x1   = r->x1;
    x2   = r->x2;
    y1   = r->y1;
    y2   = r->y2;
    maxr = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
    dr   = maxr / ((double)(r->n_intvl));
    if ((r->hist = calloc(r->n_intvl, sizeof(int))) == NULL)
	pstop("!!! xcor_spots: not enough memory.\n");
    for (j=0; j < r->n_sp; j++) {
	jdx = r->sp[j];
	x2  = p->x[jdx];
	y2  = p->y[jdx];
	for (i=0; i < r->n_sp; i++) {
	    idx = r->sp[i];
	    x1  = p->x[idx];
	    y1  = p->y[idx];
	    rr = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
	    k  = (int)(rr / dr);
	    if (k < r->n_intvl) r->hist[k] ++;
	}
    }
}

void recspot(para_t *p)
{
    int  i;

    for (i=0; i < p->n_rec; i++) {
	if (p->rec[i].update == 0) continue;

	printf("xcor: (x1,y1)=(%.4e,%.4e), (x2,y2)=(%.4e,%.4e)\n",
		p->rec[i].x1, p->rec[i].y1, p->rec[i].x2, p->rec[i].y2);
	if (p->rec[i].hist != NULL) free(p->rec[i].hist);
	if (p->rec[i].sp   != NULL) free(p->rec[i].sp);
	p->rec[i].hist = NULL;
	p->rec[i].sp   = NULL;
	p->rec[i].n_sp = 0;
	p->rec[i].m_sp = 0;

	find_spots(p, p->rec+i);
	xcor_spots(p, p->rec+i);
	p->rec[i].update = 0;
    }
}

/*--------------------------------------------------------------------------
 *
 *	Data I/O for spot correlation
 *
 *-------------------------------------------------------------------------*/

static char *xcor_token(FILE *f)
{
    static char buf[1024];
    char       *s, *s1;

    if (fgets(buf, 1023, f) == NULL)
	pstop("!!! xcor_token: reading error.\n");
    s = strchr(buf, ':');
    s = (s == NULL) ? buf : s+1;
    while ((*s == ' ' || *s == '\t') && *s != '\0')
	s ++;
    if ((s1 = strrchr(buf, '\n')) != NULL)
	*s1 = '\0';
    return (s[0] == '\0') ? NULL : s;
}

static void xcor_skip_newline(FILE *f)
{
    int  c;

    while ((c = fgetc(f)) != '\0' && c != '\n') ;
}

void xcor_input(para_t *p)
{
    char  *fn, *s;
    FILE  *f;
    int    i, j;
    double rr;

    fn = "xcor.dat";
    if ((f = fopen(fn, "rt")) == NULL) return;

    s = xcor_token(f);
    sscanf(s, "%d", &(p->n_rec));
    p->m_rec = p->n_rec;
    if ((p->rec = malloc(p->m_rec * sizeof(rec_t))) == NULL)
	pstop("!!! xcor_input: not enough memory.\n");
    xcor_token(f);
    
    for (i=0; i < p->n_rec; i++) {
	s = xcor_token(f);		// REC_ID:
	s = xcor_token(f);		// REC_C1:
	sscanf(s, "%lf %lf", &(p->rec[i].x1), &(p->rec[i].y1));
	s = xcor_token(f);		// REC_C2:
	sscanf(s, "%lf %lf", &(p->rec[i].x2), &(p->rec[i].y2));
	s = xcor_token(f);		// REC_NP:
	sscanf(s, "%d", &(p->rec[i].n_sp));
	s = xcor_token(f);		// REC_NH:
	sscanf(s, "%d", &(p->rec[i].n_intvl));
	s = xcor_token(f);		// REC_F0:
	sscanf(s, "%d", &(p->rec[i].Fstart));
	s = xcor_token(f);		// REC_DR:

	fseek(f, 7, SEEK_CUR);		// REC_SP:
	if ((p->rec[i].sp = malloc(p->rec[i].n_sp * sizeof(int))) == NULL)
	    pstop("!!! xcor_input: not enough memory for sp: %d\n", i);
	for (j=0; j < p->rec[i].n_sp; j++)
	    fscanf(f, " %d", p->rec[i].sp+j);
	xcor_skip_newline(f);
	    
	p->rec[i].m_sp   = p->rec[i].n_sp;
	p->rec[i].update = 0;
	if (p->rec[i].n_intvl != p->n_intvl) {
	    p->rec[i].update  = 1;
	    p->rec[i].n_intvl = p->n_intvl;
	    p->rec[i].hist    = NULL;
	}
	if (p->rec[i].Fstart != p->Fstart) {
	    p->rec[i].update = 1;
	    p->rec[i].Fstart = p->Fstart;
	    p->rec[i].hist   = NULL;
	}

	if ((p->rec[i].hist = malloc(p->rec[i].n_intvl * sizeof(int))) == NULL)
	    pstop("!!! xcor_input: not enough memory.\n");
	for (j=0; j < p->rec[i].n_intvl; j++) {
	    if ((s = xcor_token(f)) == NULL)
		pstop("!!! xcor_input: REC_ID=%d reading error.\n", i);
	    sscanf(s, "%lf %d", &rr, p->rec[i].hist+j);
	}
	xcor_token(f);
	printf("read xcor: REC_ID=%d\n", i);
    }
}

void xcor_output(para_t *p)
{
    char  *fn;
    FILE  *f;
    int    i, j;
    double x1, x2, y1, y2, maxr, dr, rr;

    fn = "xcor.dat";
    if ((f = fopen(fn, "wt")) == NULL)
	return;
    fprintf(f, "N_REC:  %d\n\n", p->n_rec);
    for (i=0; i < p->n_rec; i++) {
	x1   = p->rec[i].x1;
	x2   = p->rec[i].x2;
	y1   = p->rec[i].y1;
	y2   = p->rec[i].y2;
	maxr = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
	dr   = maxr / ((double)(p->rec[i].n_intvl));

	fprintf(f, "REC_ID: %d\n", i);
	fprintf(f, "REC_C1: %.16E %.16E\n", x1, y1);
	fprintf(f, "REC_C2: %.16E %.16E\n", x2, y2);
	fprintf(f, "REC_NP: %d\n", p->rec[i].n_sp);
	fprintf(f, "REC_NH: %d\n", p->rec[i].n_intvl);
	fprintf(f, "REC_F0: %d\n", p->rec[i].Fstart);
	fprintf(f, "REC_DR: %.16E\n", dr);
	fprintf(f, "REC_SP:");
	for (j=0; j < p->rec[i].n_sp; j++)
	    fprintf(f, " %08d", p->rec[i].sp[j]);
	fprintf(f, "\n");
	for (j=0; j < p->rec[i].n_intvl; j++) {
	    rr = j * dr;
	    fprintf(f, "%.16E  %d\n", rr, p->rec[i].hist[j]);
	}
	fprintf(f, "\n");
    }
    fclose(f);
}

