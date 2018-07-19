#ifndef _PCLST_H
#define _PCLST_H

typedef struct {
    long     cID;		// cluster ID.
    double   x, y;		// spot coordinate
} spot_t;

typedef struct {
    long     np, mp;		// # of spots, buffer size
    long    *sID;		// list of spots
    double   x1[2], x2[2];	// boundary of this cluster.
    double   xc[2];		// weight center of this cluster.
} clst_t;

typedef struct {
    double   x1, x2, y1, y2;
    int     *hist;
    long    *sp, nsp;
} rec_t;

typedef struct {
    char    *spotfn;		// input spot filename
    char    *Fstsfn;		// input frame statistics filename
    double   bx1, by1;		// left-bottom corner of the image
    double   bx2, by2;		// right-top   corner of the image
    double   m_density;		// max. event density for starting frame
    int      Fstart;		// starting frame
    double   spd;		// max. distances of spots for a cluster
    int      nsp_draw;		// min. # of spots of clusters to draw
    int      nrec_draw;		// max. # of rectangles to draw
    int      ncls1;		// # of clusters to keep in the 1st stage.
    int      n_intvl;		// # of intervals for histogram
    int      verb;		// verbose message output.

    long     np, mp;		// # of spots, spot buffer size
    spot_t  *sp;		// spot list
    long     ncl, mcl;		// # of clusters, size of buffer
    clst_t  *cl;		// clusters
    long     n_rec, m_rec;	// # of rectagles, size of buffer
    rec_t   *rec;		// retangles
} para_t;

void pstop(char *fmt, ...);
void readfsts(para_t *p);
void readspot(para_t *p);
void cluster_identify(para_t *p);
void output(para_t *p);
rec_t *push_rec(para_t *p, double x1, double x2, double y1, double y2);
void xcor_spots(para_t *p, clst_t *cl, rec_t *r);
void xcor_output(para_t *p);

#endif
