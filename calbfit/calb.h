#ifndef _CALBFIT_H
#define _CALBFIT_H

#define BUFLEN 1024

// Data structure for calibration curve.
typedef struct {
    int     N;
    double *x, *f, *err;
    double  w0, c, d, A, B;
    double  dw, dc, dd, dA, dB;
} calb_t;


// Tool functions.
void   *memalloc(size_t n_bytes, void *v0);
FILE   *openfile(const char *fn, const char *mode);
void    dsysv(int N, double *A, int Nrhs, double *B, double *X);

// Fitting formulas.
void    cbfunc(void *data, double *a, double *f, double *dfda);
void    cbddf(void *data, double *f, double *df, double *ddf);

// Fitting functions.
int     nlinfit(void (*func)(void*, double*, double*, double*),
		void *data, int n_data, int na, double *a, double *da,
		double *chisq, double tol, int maxloop, int verb);
double  get_chisq(void (*func)(void*, double*, double*, double*),
		  void *data, int n_data, int na, double *a);

// calb_t functions.
calb_t *calb_readraw(char *fn);
void    calb_destroy(calb_t *cb);
void    calb_fit(int verb, calb_t *cb, char *name); 
void    calb_output(char *outfn, calb_t *cbX, calb_t *cbY);

#endif
