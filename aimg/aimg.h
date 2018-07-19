#ifndef _AIMG_H
#define _AIMG_H

#define   MAX_PIXEL  255.0

typedef struct {
    double  x, y, w, II;	// parameters for an artificial spot
} spot_t;

typedef struct {
    char          *outf;	// output image filename
    unsigned int   seed;	// random number seed
    int            iW, iH;	// image size (pixel)
    double         I1, I2;	// signal intensity range
    double         N1, N2;	// noise intensity range
    double         w1, w2;	// spot width range (pixel)
    int            np;		// number of spots
    int            spixel;	// size of spot pixel square
    int            smesh;	// mesh size of each pixel

    double         noise;	// averaged intensity of noise
    spot_t        *sp;		// list of artificial spots
    unsigned char *img;		// the generated image.
} para_t;

void pstop(char *fmt, ...);
void aimg(para_t *p);
void output_JPEG(para_t *p);
void output_spot(para_t *p);
void output_raw(para_t *p);

#endif
