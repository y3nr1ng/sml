#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tiffio.h>
#include "pix.h"

/*---------------------------------------------------------------------------
*
*	Prepare filename format for separated image frame files.
*
*--------------------------------------------------------------------------*/

static void FIO_separate(frameIO_t *fio, char *imgfn) {
    int ndig;
    char      *s1, *s2, *fnp, fnbase[1024], fn[1024];

    fio->imgf = NULL;
    strncpy(fnbase, imgfn, 1023);
    if ((s1 = strstr(fnbase, "XX")) != NULL) {
        s2 = s1+2;
        while (*s2 == 'X') s2++;
        ndig = (int)(s2 - s1);
        *s1 = '\0';
        sprintf(fn, "%s%%0%dd%s", fnbase, ndig, s2);
        fnp = fn;
    }else
        fnp = imgfn;

    if ((fio->basefn = strdup(fnp)) == NULL)
        pstop("!!! FIO_separate: not enough memory.\n");
}

/*---------------------------------------------------------------------------
*
*	Frame I/O from TIFF files.
*
*--------------------------------------------------------------------------*/

static void FIO_tiff(frameIO_t *fio, char *imgfn) {
    TIFF *tif;

    TIFFSetWarningHandler(NULL);
    if ((tif = TIFFOpen(fio->basefn, "r")) == NULL)
        pstop("!!! frameIO: cannot open tiff file: %s\n", fio->basefn);
    fio->imgf = (void *)tif;
    if ((fio->basefn = strdup(imgfn)) == NULL)
        pstop("!!! FIO_tiff: not enough memory.\n");
}

static void tiff_info(para_t *p, frameIO_t *fio) {
    TIFF *tif;
    int idx=1024, sdx=1024, r=0, w, h;

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

    if (TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w) != 1)
        pstop("!!! TIFF: cannot get tag: TIFFTAG_IMAGEWIDTH\n");
    if (TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h) != 1)
        pstop("!!! TIFF: cannot get tag: TIFFTAG_IMAGELENGTH\n");

    p->tot_frames = idx;
    p->img_W = w;
    p->img_H = h;
}

static frameloc_t *FR_tiff(para_t *p, frameIO_t *fio, int idx) {
    tsize_t     s_strip;
    tstrip_t    n_strip;
    frameloc_t *fm=NULL;
    TIFF       *tif=NULL;
    uint32      w, h, i, psize;
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
        pstop("!!! TIFF %d: not enough memory.\n", idx);
    if (TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w) != 1)
        pstop("!!! TIFF %d: cannot get tag: TIFFTAG_IMAGEWIDTH\n", idx);
    if (TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h) != 1)
        pstop("!!! TIFF %d: cannot get tag: TIFFTAG_IMAGELENGTH\n", idx);
    psize = (size_t)(s_strip*n_strip) / (size_t)(w*h);
    if (psize != 2)
        pstop("!!! TIFF: %d: pixel size should be 2 bytes.\n", idx);

    for (i=0; i < h; i++) {
        TIFFReadScanline(tif, dp, i, 1);
        dp += (w*psize);
    }
    mx.dim_x = w;
    mx.dim_y = h;
    mx.data  = (void *)data;
    fm = frameCreate(p, idx, &mx, 's');
    _TIFFfree(data);

    return fm;
}

/*---------------------------------------------------------------------------
*
*	Frame I/O from RAW data files.
*
*--------------------------------------------------------------------------*/

static void raw_info(para_t *p, frameIO_t *fio) {
    FILE  *f;
    char buf[4096];
    int xmax=0, ymax=0, x, y, II;

    sprintf(buf, fio->basefn, p->frameID1);
    if ((f = fopen(buf, "rt")) == NULL)
        pstop("!!! raw_info: cannot open file: %s.\n", buf);
    fgets(buf, 4096, f);
    fgets(buf, 4096, f);
    while (fgets(buf, 4096, f) != NULL) {
        sscanf(buf, "%d %d %d", &x, &y, &II);
        xmax = (xmax < x) ? x : xmax;
        ymax = (ymax < y) ? y : ymax;
    }
    xmax++;
    ymax++;
    fclose(f);

    p->tot_frames = 1;
    p->frameID1   = 0;
    p->frameID2   = 0;
    p->img_W      = xmax;
    p->img_H      = ymax;
}

static frameloc_t *FR_raw(para_t *p, frameIO_t *fio, int idx) {
    frameloc_t *fm = NULL;
    matmx_t mx;
    FILE       *f;
    short      *data;
    char buf[4096];
    int w, h, i, x, y, II;

    sprintf(buf, fio->basefn, idx);
    if ((f = fopen(buf, "rt")) == NULL)
        pstop("!!! FR_raw: cannot open file: %s.\n", buf);

    w = p->img_W;
    h = p->img_H;
    if ((data = malloc(w*h*sizeof(short))) == NULL)
        pstop("!!! FR_raw: not enough memory.\n");

    fseek(f, 0, SEEK_SET);
    fgets(buf, 4096, f);
    fgets(buf, 4096, f);
    for (i=0; i < w*h; i++) {
        fgets(buf, 4096, f);
        sscanf(buf, "%d %d %d", &x, &y, &II);
        data[y+x*h] = (short)II;
    }
    fclose(f);

    mx.dim_x = w;
    mx.dim_y = h;
    mx.data  = (void *)data;
    fm = frameCreate(p, idx, &mx, 's');

    return fm;
}


/*----------------------------------------------------------------------------
*
*	Frame I/O driver routings
*
*---------------------------------------------------------------------------*/

frameIO_t *frameIO_Open(para_t *p) {
    frameIO_t *fio;

    if ((fio = malloc(sizeof(frameIO_t))) == NULL)
        pstop("!!! frameIO_ReOpen: not enough memory.\n");
    if ((fio->basefn = strdup(p->imgfn)) == NULL)
        pstop("!!! frameIO_ReOpen: not enough memory.\n");
    fio->type = p->imgfmt;

    switch (p->imgfmt) {
    case 1:
        FIO_tiff(fio, p->imgfn);
        break;

    case 2:
        FIO_separate(fio, p->imgfn);
        break;

    default:
        pstop("!!! frameIO_Open: unknown image file format: %d\n", p->imgfmt);
    }
    return fio;
}

void frameIO_Close(para_t *p, frameIO_t *fio) {
    if (fio->imgf) {
        switch (p->imgfmt) {
        case 1:
            TIFFClose((TIFF *)(fio->imgf));
            break;
        }
    }
    fio->imgf = NULL;
    free(fio->basefn);
    free(fio);
}

void frameIO_init(para_t *p) {
    frameIO_t *fio;
    int ntotal;
//
//  Open the image frames data file to obtain the frame information.
//
    fio = frameIO_Open(p);
    switch (p->imgfmt) {
    case 1:
        tiff_info(p, fio);
        break;

    case 2:
        raw_info(p, fio);
        break;
    }
    frameIO_Close(p, fio);
//
//  Initialize the frame parameters.
//
    if (p->frameID1 < 0) p->frameID1=0;
    if (p->frameID2 < 0) p->frameID2=p->tot_frames-1;
    if (p->frameID1 >= p->tot_frames || p->frameID2 >= p->tot_frames)
        pstop("!!! frameIO: frameID out of range, max=%d\n", p->tot_frames-1);
    if (p->frame_x1 < 0 || p->frame_x2 <= p->frame_x1) {
        p->frame_x1 = 0;
        p->frame_x2 = p->img_W-1;
    }
    if (p->frame_y1 < 0 || p->frame_y2 <= p->frame_y1) {
        p->frame_y1 = 0;
        p->frame_y2 = p->img_H-1;
    }
    p->max_x   = p->frame_x2 - p->frame_x1;
    p->max_y   = p->frame_y2 - p->frame_y1;
    p->max_dx  = p->max_dx  / p->nm_px_x;
    p->max_dy  = p->max_dy  / p->nm_px_y;
    p->max_dwx = p->max_dwx / p->nm_px_x;
    p->max_dwy = p->max_dwy / p->nm_px_y;

    ntotal = p->frameID2 - p->frameID1 + 1;
    if ((p->fsts = calloc(ntotal, sizeof(framests_t))) == NULL)
        pstop("!!! frameIO_init: not enough memory for frame sts.\n");

    ntotal = (p->frame_x2 - p->frame_x1 + 1)*(p->frame_y2 - p->frame_y1 + 1);
    if ((p->psum = calloc(ntotal, sizeof(int))) == NULL)
        pstop("!!! frameIO_init: not enough memory for frame sum.\n");

    printf("Total frames:  %d\n", p->tot_frames);
    printf("Image size:    %d x %d\n", p->img_W, p->img_H);
    printf("Image ROI:  x=(%d,%d), y=(%d,%d)\n",
           p->frame_x1, p->frame_x2, p->frame_y1, p->frame_y2);
    fflush(stdout);
}

frameloc_t *frameIO(para_t *p, frameIO_t *fio, int idx) {
    frameloc_t *fm=NULL;

    switch (p->imgfmt) {
    case 1:
        fm = FR_tiff(p, fio, idx);
        break;

    case 2:
        fm = FR_raw(p, fio, idx);
        break;
    }
    return fm;
}
