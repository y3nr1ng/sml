#ifndef _PSPOT_H
#define _PSPOT_H

typedef struct {
    double  x1, x2, y1, y2;
    int     update, n_intvl;
    int     n_sp, m_sp, *sp;
    int    *hist, Fstart;
} rec_t;

typedef struct {
    int     datafmt;		// data format: 0:spot, 1:spotmesh, 2:pixel
    char   *inpfn;		// input data filename
    char   *inpfnF;		// input frame sts filename
    char   *outmesh;		// output spot meshed image filename
    double  m_density;		// max. event density for starting frame
    int     imgscale;		// image size: 0:full, 1:scaled
    int     imgWidth;		// image size (pixel)
    int     imgHeight;		// image size (pixel)
    double  iW, iH;             // image range for spot (nm)
    int     n_intvl;		// n_intvls for histogram for xcor
    double  mesh_size;		// grid size of the mesh for spot (nm)
    int     Fstart;		// starting frame
    int     view;		// view the image or not
    char   *fontfn;		// font filename.

    int     mp;                 // buffer size of number of spots
    int     np;                 // number of spots
    double *x, *y;              // coordinate of spots
    double *pix;		// pixel array for gray image
    double  pmax, pmin;		// max and min value of pixels
    double  gamma;		// gamma correction for gray scale image.
    unsigned char *gimg;	// the generated gray scale image.

    int     n_rec, m_rec;	// count of rectanglars for xcor
    rec_t  *rec;
} para_t;

void readfsts(para_t *p);
void readspot(para_t *p);
void readpixel(para_t *p);
void spot_image(para_t *p);
void gray_image(para_t *p);
void mesh_events(para_t *p);
void mesh_output(para_t *p);
void output_JPEG(para_t *p, int color_mode, char *fn);
void display_image(int argc, char **argv, para_t *p);
void pstop(char *fmt, ...);

void push_rec(para_t *p, double x1, double x2, double y1, double y2);
void del_rec(para_t *p, double x, double y);
rec_t *search_rec(para_t *p, double x, double y);
void recspot(para_t *p);
void xcor_input(para_t *p);
void xcor_output(para_t *p);

#endif
