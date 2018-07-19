#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tiffio.h>
#include "pix.h"

/*---------------------------------------------------------------------------
 *
 *	Distribute frames to each threads or parallel processes
 *
 *--------------------------------------------------------------------------*/

void frameDistribution(int myid, int nprc, int tot_fID0, int tot_fIDN,
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

/*---------------------------------------------------------------------------
 *
 *	Frame I/O initialization.
 *
 *--------------------------------------------------------------------------*/

static void FIO_matlab(frameIO_t *fio, char *imgfn)
{
    int        ndig;
    char      *s1, *s2, fn[1024], fmt[1024];

    fio->imgf = NULL;
    strncpy(fn, imgfn, 1023);
    if ((s1 = strstr(fn, "XX")) != NULL) {
	s2 = s1+2;
	while (*s2 == 'X') s2++;
	ndig = (int)(s2 - s1);
	*s1 = '\0';
	sprintf(fmt, "%s%%0%dd%s", fn, ndig, s2);
	if ((fio->basefn = strdup(fmt)) == NULL)
	    pstop("!!! FIO_matlab: not enough memory.\n");
    }
    else {
	fio->basefn = imgfn;
    }
}

static void FIO_tiff(frameIO_t *fio, char *imgfn)
{
    TIFF *tif;

    fio->basefn = imgfn;
    TIFFSetWarningHandler(NULL);
    if ((tif = TIFFOpen(fio->basefn, "r")) == NULL)
        pstop("!!! frameIO: cannot open tiff file: %s\n", fio->basefn);
    fio->imgf = (void *)tif;
}

static void tiff_Nframes(frameIO_t *fio, para_t *p)
{
    TIFF *tif;
    int   idx=1024, sdx=1024, r=0, fID1, fID2;

    fID1 = p->frameID1;
    fID2 = p->frameID2;
    if (fID1 < 0)
	fID1 = 0;
    if (fID2 < 0) {
	tif = (TIFF *)(fio->imgf);
	while (TIFFSetDirectory(tif, (tdir_t)idx) == 1) idx += sdx;
	while (sdx > 1) {
	    sdx /= 2;
	    idx  = (r == 0) ? idx-sdx : idx+sdx;
	    r = TIFFSetDirectory(tif, (tdir_t)idx);
	}
	if (r == 0)
	    idx = idx-sdx;
	else if (TIFFSetDirectory(tif, (tdir_t)(idx+sdx)) == 1)
	    idx = idx+sdx;
	fID2 = idx;
	printf("%d: Total available frames: %d\n", p->myid, fID2-fID1+1);
	fflush(stdout);
	TIFFClose(tif);
	fio->imgf = NULL;
    }
    frameDistribution(p->myid, p->nprc, fID1, fID2, &(p->frameID1),
		      &(p->frameID2));
}

static void FIO_raw(frameIO_t *fio, char *imgfn)
{
    FILE *f;

    if ((f = fopen(imgfn, "rt")) == NULL)
	pstop("!!! FIO_raw: cannot open file: %s.\n", imgfn);
    fio->basefn = imgfn;
    fio->imgf   = (void *)f;
}

void frameIO_init(para_t *p)
{
    frameIO_t *fio;
    int  ntotal;

    if ((fio = malloc(sizeof(frameIO_t))) == NULL)
	pstop("!!! frameIO_init: not enough memory.\n");
    fio->type = p->imgfmt;

    switch (p->imgfmt) {
    case 0:
//	FIO_matlab(fio, p->imgfn);
	break;

    case 1:
	FIO_tiff(fio, p->imgfn);
	tiff_Nframes(fio, p);
	break;

    case 2:
//	FIO_raw(fio, p->imgfn);
	p->frameID1 = 0;
	p->frameID2 = 0;
	break;

    default:
	pstop("!!! frameIO: unknown image file format: %d\n", p->imgfmt);
    }
    ntotal = p->frameID2 - p->frameID1 + 1;
    if ((p->fsts = calloc(ntotal, sizeof(framests_t))) == NULL)
	pstop("!!! frameIO_init: not enough memory for frame sts.\n");

    ntotal = (p->frame_x2 - p->frame_x1 + 1)*(p->frame_y2 - p->frame_y1 + 1);
    if ((p->psum = calloc(ntotal, sizeof(int))) == NULL)
	pstop("!!! frameIO_init: not enough memory for frame sum.\n");
}

frameIO_t *frameIO_ReOpen(para_t *p)
{
    frameIO_t *fio;

    if ((fio = malloc(sizeof(frameIO_t))) == NULL)
	pstop("!!! frameIO_ReOpen: not enough memory.\n");
    fio->type = p->imgfmt;

    switch (p->imgfmt) {
    case 0:
	FIO_matlab(fio, p->imgfn);
	break;

    case 1:
	FIO_tiff(fio, p->imgfn);
	break;

    case 2:
	FIO_raw(fio, p->imgfn);
	break;
    }
    return fio;
}

void frameIO_Close(para_t *p, frameIO_t *fio)
{
    switch (p->imgfmt) {
    case 1:
	TIFFClose((TIFF *)(fio->imgf));
	break;
    }
    free(fio);
}

/*----------------------------------------------------------------------------
 *
 *	Frame reading.
 *
 *---------------------------------------------------------------------------*/

static frameloc_t *FR_matlab(para_t *p, frameIO_t *fio, int idx)
{
    frameloc_t *fm=NULL;
    matfile_t  *matf=NULL;
    char        fn[1024];

    sprintf(fn, fio->basefn, idx);
    matf = matRead(fn);
    fm   = frameCreate(p, idx, 0, &(matf->mx), 's');
    matFree(matf);

    return fm;
}

static frameloc_t *FR_tiff(para_t *p, frameIO_t *fio, int idx)
{
    frameloc_t *fm=NULL;
    TIFF       *tif=NULL;
    tsize_t     s_strip;
    tstrip_t    n_strip;
    uint32     *bc, w, h, i;
    char       *data, *dp;
    matmx_t     mx;

    tif = (TIFF *)fio->imgf;
    if (TIFFSetDirectory(tif, (tdir_t)idx) != 1)
        pstop("TIFF: cannot set directory: %d\n", idx);
    s_strip = TIFFStripSize(tif);
    n_strip = TIFFNumberOfStrips(tif);
    data    = _TIFFmalloc(s_strip*n_strip);
    dp      = data;
    if (data == NULL)
	pstop("!!! TIFF: not enough memory.\n");
    if (TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w) != 1)
        pstop("!!! TIFF: cannot get tag: TIFFTAG_IMAGEWIDTH\n");
    if (TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h) != 1)
        pstop("!!! TIFF: cannot get tag: TIFFTAG_IMAGELENGTH\n");
    if (TIFFGetField(tif, TIFFTAG_STRIPBYTECOUNTS, &bc) != 1)
        pstop("!!! TIFF: cannot get tag: TIFFTAG_STRIPBYTECOUNTS\n");
    for (i=0; i < n_strip; i++) {
        if (TIFFReadRawStrip(tif, i, dp, bc[i]) != bc[i])
	    pstop("!!! TIFF: cannot read completed strip.\n");
        dp += s_strip;
    }
    mx.dim_x = w;
    mx.dim_y = h;
    mx.data  = (void *)data;
    fm = frameCreate(p, idx, 1, &mx, 's');
    _TIFFfree(data);
    return fm;
}

static frameloc_t *FR_raw(para_t *p, frameIO_t *fio, int idx)
{
    frameloc_t *fm = NULL;
    matmx_t     mx;
    FILE       *f  = (FILE *)(fio->imgf);
    short      *data;
    char        buf[4096];
    int         xmax=0, ymax=0, i, x, y, II;

    if (f == NULL) return NULL;

    fgets(buf, 4096, f);
    fgets(buf, 4096, f);
    while (fgets(buf, 4096, f) != NULL) {
	sscanf(buf, "%d %d %d", &x, &y, &II);
	xmax = (xmax < x) ? x : xmax;
	ymax = (ymax < y) ? y : ymax;
    }
    xmax ++;
    ymax ++;
    if ((data = malloc(xmax*ymax*sizeof(short))) == NULL)
	pstop("!!! FR_raw: not enough memory.\n");

    fseek(f, 0, SEEK_SET);
    fgets(buf, 4096, f);
    fgets(buf, 4096, f);
    for (i=0; i < xmax*ymax; i++) {
	fgets(buf, 4096, f);
	sscanf(buf, "%d %d %d", &x, &y, &II);
	data[y+x*ymax] = (short)II;
    }
    fclose(f);
    fio->imgf = NULL;

    mx.dim_x = xmax;
    mx.dim_y = ymax;
    mx.data  = (void *)data;
    fm = frameCreate(p, idx, 1, &mx, 's');

    return fm;
}

frameloc_t *frameIO(para_t *p, frameIO_t *fio, int idx)
{
    frameloc_t *fm=NULL;

    switch (p->imgfmt) {
    case 0:
	fm = FR_matlab(p, fio, idx);
	break;

    case 1:
	fm = FR_tiff(p, fio, idx);
	break;

    case 2:
	fm = FR_raw(p, fio, idx);
	break;
    }
    return fm;
}
