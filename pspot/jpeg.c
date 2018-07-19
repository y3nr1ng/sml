#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <jpeglib.h>
#include "pspot.h"

/*-------------------------------------------------------------------------
 *
 *	Output the JPEG image file from a given image buffer.
 *
 *      color_mode=0:  RGB format
 *      color_mode=1:  Gray scale format
 *
 *------------------------------------------------------------------------*/

void output_JPEG(para_t *p, int color_mode, char *fn)
{
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr       jerr;
    FILE    *f;
    char     fn2[1024], *s;
    JSAMPROW row_pointer[1];
    int      iW, iH, quality, row_stride;

    iW = p->imgWidth;
    iH = p->imgHeight;
    quality = 100;
    strncpy(fn2, fn, 1023);
    if ((s = strrchr(fn2, '.')) != NULL)
	*s = '\0';
    strncat(fn2, ".jpg", 1023);

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);

    if ((f = fopen(fn2, "wb")) == NULL)
	pstop("output_JPEG: cannot open file: %s\n", fn2);
    jpeg_stdio_dest(&cinfo, f);
    cinfo.image_width  = iW;
    cinfo.image_height = iH;
    if (color_mode == 0) {
	cinfo.input_components = 3;
	cinfo.in_color_space   = JCS_RGB;
    }
    else {
	cinfo.input_components = 1;
	cinfo.in_color_space   = JCS_GRAYSCALE;
    }
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE);

    jpeg_start_compress(&cinfo, TRUE);
    row_stride = iW * cinfo.input_components;
    while (cinfo.next_scanline < cinfo.image_height) {
	row_pointer[0] = &(p->gimg[cinfo.next_scanline*row_stride]);
	jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }
    jpeg_finish_compress(&cinfo);
    fclose(f);

    jpeg_destroy_compress(&cinfo);
}
