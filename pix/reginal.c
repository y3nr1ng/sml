#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pix.h"

#define UNSET 255

typedef struct {
    int    dim_x, dim_y;	// dimension of the image.
    int    n_nei;		// # of neighbors (4 or 8).
    int    small;		// indicate that the cluster is small or not.
    int   *pos, n_pos;		// pixels belong to this cluster.
    int   *n1, *n2;		// # of forward / backword neighbors.
    int   *nei1, *nei2;		// indice of forward / backward neighbors.
    short *chk, *res;
    char   imgtype;		// i:int, s:short
    void  *image;
} regmm_t;
 
/*
 *  Get the neighbor pixels.
 *      n1, nei1[]: neighbor pixels in the forward directions.
 *      n2, nei2[]: neighbor pixels in the backward directions.
 *
 *  n_nei==4:      d          n_nei==8:     f g h
 *               c x a                      e x a
 *                 b                        d c b
 */
static void get_neighbor(regmm_t *rg)
{
    int  nn1, nn2, x, y, dim_x, dim_y;
    int *n1, *n2, *nei1, *nei2;

    dim_x = rg->dim_x;
    dim_y = rg->dim_y;
    n1    = rg->n1;
    n2    = rg->n2;
    nei1  = rg->nei1;
    nei2  = rg->nei2;

    for (y=0; y < dim_y; y++) {
    for (x=0; x < dim_x; x++) {
	nn1 = 0;
	nn2 = 0;

	if (x+1 < dim_x)			// forward neighbors.
	    nei1[nn1++] = (x+1) + y*dim_x;
	if (y+1 < dim_y)
	    nei1[nn1++] = x + (y+1)*dim_x;
	if (rg->n_nei == 8) {
	    if (x+1 < dim_x && y+1 < dim_y)
		nei1[nn1++] = (x+1) + (y+1)*dim_x;
	    if (x+1 < dim_x && y-1 >= 0)
		nei1[nn1++] = (x+1) + (y-1)*dim_x;
	}

	if (x-1 >= 0)				// backward neighbors.
	    nei2[nn2++] = (x-1) + y*dim_x;
	if (y-1 >= 0)
	    nei2[nn2++] = x + (y-1)*dim_x;
	if (rg->n_nei == 8) {
	    if (x-1 >= 0 && y+1 < dim_y)
		nei2[nn2++] = (x-1) + (y+1)*dim_x;
	    if (x-1 >= 0 && y-1 >= 0)
		nei2[nn2++] = (x-1) + (y-1)*dim_x;
	}

	*n1 = nn1;
	*n2 = nn2;
	n1 ++;
	n2 ++;
	nei1 += 4;
	nei2 += 4;
    }}
}

/*
 *  Check the pixel xx with its neighbors to determine whether it is
 *  a regional maximum.
 */
static void
check_neighbor(regmm_t *rg, int xx)
{
    int    i, yy, *pos;
    int   *imgi, c, d;
    short *imgs, *res, *chk;
    char   type;

    imgi = (int *)(rg->image);
    imgs = (short *)(rg->image);
    type = rg->imgtype;
    res  = rg->res;
    chk  = rg->chk;
    pos  = rg->pos;
    if (res[xx] != UNSET) return;
    c = (type == 'i') ? imgi[xx] : (int)(imgs[xx]);

    for (i=0; i < rg->n1[xx] && res[xx] == UNSET; i++) {
	yy = rg->nei1[i+xx*4];
	d  = (type == 'i') ? imgi[yy] : (int)(imgs[yy]);
	if (c < d) {
	    rg->small ++;
	    res[xx] = 0;
	}
	else if (c == d && chk[yy] == UNSET) {
	    pos[rg->n_pos ++] = yy;
	    check_neighbor(rg, yy);
	    if (res[xx] == UNSET && res[yy] != UNSET)
		res[xx] = res[yy];
	}
	chk[yy] = 1;
    }
    if (res[xx] == UNSET) {
	for (i=0; i < rg->n2[xx]; i++) {
	    yy = rg->nei2[i+xx*4];
	    d  = (type == 'i') ? imgi[yy] : (int)(imgs[yy]);
	    if (c < d) {
		rg->small ++;
	        res[xx] = 0;
	        break;
	    }
	    else if (c == d && res[yy] != UNSET) {
		res[xx] = res[yy];
		break;
	    }
	}
    }
}

/*-------------------------------------------------------------------------
 *
 *  Get the regional maximum (API)
 *
 *------------------------------------------------------------------------*/

void regional_max(int dim_x, int dim_y, int n_neighbor, char imgtype,
		  void *image, short *res)
{
    int      xx, yy, x, y, i;
    int     *pos, *nn;
    short   *chk;
    regmm_t  rg;

    chk = malloc(dim_x*dim_y*sizeof(short));		// record checked pixels
    pos = malloc(dim_x*dim_y*sizeof(int));		// pixels in cluster
    nn  = malloc(dim_x*dim_y*10*sizeof(int));
    if (chk==NULL || pos==NULL || nn==NULL)
	pstop("!!! regional_max: not enough memory.\n");
    for (xx=0; xx < dim_x*dim_y; xx++)
	res[xx] = UNSET;
    rg.dim_x   = dim_x;
    rg.dim_y   = dim_y;
    rg.image   = image;
    rg.imgtype = imgtype;
    rg.pos     = pos;
    rg.chk     = chk;
    rg.res     = res;
    rg.n_nei   = n_neighbor;
    rg.n1      = nn;
    rg.n2      = nn+dim_x*dim_y;
    rg.nei1    = nn+dim_x*dim_y*2;
    rg.nei2    = nn+dim_x*dim_y*6;
    get_neighbor(&rg);
    
    for (y=0; y < dim_y; y++) {
    for (x=0; x < dim_x; x++) {
	memcpy(chk, res, dim_x*dim_y*sizeof(short));
	rg.n_pos = 0;
	rg.small = 0;

	xx = x+y*dim_x;
	check_neighbor(&rg, xx);	// check neighbor pixels.
	if (res[xx] == UNSET)
	    res[xx] = (short)((rg.small > 0) ? 0 : 1);
	for (i=0; i < rg.n_pos; i++) {
	    yy = pos[i];		// set all pixels in the cluster.
	    if (res[yy] == UNSET)
		res[yy] = res[xx];
	}
    }}
    free(chk);
    free(pos);
    free(nn);
}
