#ifndef _PIX_H
#define _PIX_H

/*-------------------------------------------------------------------------
 *
 *  Image frame I/O.
 *
 *------------------------------------------------------------------------*/

typedef struct {                // structure for the MAT-Matrix variable.
    unsigned int m_class;
    int          type;
    int          dim_x, dim_y;
    char        *name;
    void        *data;
} matmx_t;

typedef struct {                // main structure for MAT-file.
    char    *fn;
    char     desc[128];
    int      version;
    int      endian;
    matmx_t  mx;
} matfile_t;

typedef struct {
    int      type;
    char    *basefn;
    void    *imgf;
} frameIO_t;

/*-------------------------------------------------------------------------
 *
 *  Calibration curve for solving z-coordinate of 3D algorithm.
 *
 *------------------------------------------------------------------------*/

typedef struct {		// Coefficient for the calibration curve
    double   w0, dw0;		//   to find the z position of spots.
    double   A[2], dA;
    double   B[2], dB;
    double   c[2], dc;
    double   d[2], dd;

    double   WW, dW;		// obtained from spot fitting.
} calb3D_t;

/*-------------------------------------------------------------------------
 *
 *  Image frame structure.
 *
 *------------------------------------------------------------------------*/

typedef struct {		// main structure for frame localization.
    int      ID;		//   ID of this frame.
    int      dim_x, dim_y;	//   dimension of the ROI image.
    char     type;		//   the data type of the frame.
    void    *frame;		//   the image of this frame.
} frameloc_t;

/*-------------------------------------------------------------------------
 *
 *  The candidate spot.
 *
 *------------------------------------------------------------------------*/

typedef struct {
    double    a[6], da[6];	// fitting parameters.
    double    chisq;		// fitting chisq/ndof.
} fitres_t;

typedef struct {
    int       x, y;		// spot coordinate.
    int       cnt;		// number of accumulates
    int       fID;		// the last frame index
    int      *img;		// pointer to frame image.
    fitres_t *res;		// the fitting results.
} sp_t;

/*-------------------------------------------------------------------------
 *
 *  For frame statistics.
 *
 *------------------------------------------------------------------------*/

typedef struct {
    int   n_event;		// number of events (spots).
} framests_t;

/*-------------------------------------------------------------------------
 *
 *  Code parameters.
 *
 *------------------------------------------------------------------------*/

typedef struct {
    int         imgfmt;		// image file format.
    char       *imgfn;		// image filename.
    char       *outfn;		// output filename (for normal spots).
    char       *outfnH;		// output filename (for high intensity spots).
    char       *outfnp;		// output filename (for spot pixel list).
    char       *outfnF;		// output filename (for frame event statistics).
    char       *outfnS;		// output filename (frame summation).
    int         frameID1;	// starting frame ID.
    int         frameID2;	// stopping frame ID.
    int         frame_x1;	// starting x for a frame.
    int         frame_x2;	// stopping x for a frame.
    int         frame_y1;	// starting y for a frame.
    int         frame_y2;	// stopping y for a frame.
    int         x_find_pixels;	// x range for finding pixels.
    int         y_find_pixels;	// y range for finding pixels.
    double      threshold1;	// lower bound threshold for pixels.
    double      threshold2;	// upper bound threshold for pixels.
    int         nfsep;		// max separation of frames of each particle
    double      nm_px_x;	// nm per pixel (x-direction).
    double      nm_px_y;	// nm per pixel (y-direction).
    double      i_photon;	// factor to convert intensity to photon
    int         mode;		// 0:2D, 1:3D
    int         rmode;		// 0:spot fitting, 1:fit & output spot
    int         alg;		// 0:reginal max; 1:max intensity
    double      max_dx;		// maximum allowed dx in fitting
    double      max_dy;		// maximum allowed dy in fitting
    double      max_dwx;	// maximum allowed dw(x) in fitting
    double      max_dwy;	// maximum allowed dw(y) in fitting
    double      min_SN;		// minimum allowed S/N ratio
    double      max_dI_I;	// maximum allowed dII/II in fitting
    double      z1, z2, dz;	// parameters for bisection solver
    calb3D_t    cax;		// parameters for calibration curve for wx
    calb3D_t    cay;		// parameters for calibration curve for wy
    int         verb;		// verbose message

    int         tot_frames;	// total number of frames.
    int         img_W, img_H;   // image width, image height

    int         m_sp1;		// current buffer size of sp (normal spots)
    int         n_sp1;		// current number of sp (for normal spots)
    sp_t      **sp1;		// candidate spots (for normal spots)
    int         m_sp2;		// current buffer size of sp (high intensity)
    int         n_sp2;		// current number of sp (for high intensity)
    sp_t      **sp2;		// candidate spots (for high intensity)

    framests_t *fsts;		// for frame statistics.
    int        *psum;		// sum over all pixels.
    double      max_x;		// width of the image (nm).
    double      max_y;		// width of the image (nm).
} para_t;


/*-------------------------------------------------------------------------
 *
 *  Subroutines.
 *
 *------------------------------------------------------------------------*/

void spot_dframe(para_t *p);
void spot_sframe(para_t *p);
void spot_fitting(para_t *p);
void spot_output_img(para_t *p);
int  spot_cmp(const void *a, const void *b);

frameIO_t  *frameIO_Open(para_t *p);
void        frameIO_Close(para_t *p, frameIO_t *fio);
void        frameIO_init(para_t *p);
frameloc_t *frameIO(para_t *p, frameIO_t *fio, int idx);
matfile_t  *matRead(char *fn);
void        matFree(matfile_t *m);

frameloc_t *frameCreate(para_t *p, int idx, matmx_t *mx, char type);
void  frameDelete(frameloc_t *fm);
void  frameSpots(para_t *p, frameloc_t *fm0, frameloc_t *fm1);
void  frameSpot2(para_t *p, frameloc_t *fm0, frameloc_t *fm1);
int   SpotFit(para_t *p, double *x_fit, double *y_fit, sp_t *sp);
int   solve_z_w(para_t *p, calb3D_t *fd, double w, double dw,
		char *v, char *dv);
int   solve_z_wxowy(para_t *p, double *a, double *da, char *v, char *dv);
void  frame_sum(para_t *p, frameloc_t *fm);

void  out_spotlist(char *outfn, int nsp, int sqsz, sp_t **sp);
FILE *out_fit_init(para_t *p, char *outfn);
void  out_fit(para_t *p, FILE *f, int Sid, sp_t *sp);
void  out_fit_close(FILE *f);
void  out_framests(para_t *p);
void  out_framesum(para_t *p);

void  regional_max(int dim_x, int dim_y, int n_neighbor, char imgtype,
		   void *image, short *res);
int   nlinfit(int n_data, void (*func)(int, void*, double*, double*, double*),
              void *data, int na, double *a, double *da, double *chisq,
	      double tol, int maxloop, int verb);
void  dsysv_(char *, int *, int *, double *, int *, int *,
	      double *, int *, double *, int *, int *, int);

void *vdup(void *v, int vlen);
int   vmax_s(int vlen, short *v, int *i_max);
void  vmaxmin_i(int vlen, int *v, int *vmax, int *i_max, int *vmin, int *i_min);
void  vadd_s(int vlen, int op, short *v1, short *v2, short *r);
void  mx_sub_s(int dim_x, int dim_y, int *sdim, short *mx, short *mr);
void  mx_subT_s(int dim_x, int dim_y, int *sdim, short *mx, short *mr);
void  mx_rsub_s(int dim_x, int dim_y, int *sdim, short *mx, short *mr);
int   mx_inv(int N, double *A, int NRHS, double *B);

int    skipline(FILE *f);
int    inp_cutcmt(char *buf);
void   inp_getINT(char *buf, int *iarg, int n, int nline);
void   inp_getDBL(char *buf, double *darg, int n, int nline);
void   inp_getSTR(char *buf, char **pstr, int nline);
double get_realtime(void);
void   pstop(char *fmt, ...);

#ifdef DEBUG
void  spotdump(int sxdim, int sydim, int *sp);
void  get_md5sum(const void *buf, size_t n);
#endif

#endif
