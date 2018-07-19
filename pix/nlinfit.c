#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pix.h"

/*-------------------------------------------------------------------------
 *
 *	Non-linear non-correlated fit.
 *
 *	Note that the weighting sig is absorbed in fmod and dfda, so
 *	that they can be handled by the user supplied function:
 *
 *		fmod[i] -> fmod[i] / sqrt(sig[i]);
 *		dfda[i] -> dfda[i] / sqrt(sig[i]);
 *
 *------------------------------------------------------------------------*/

static double
mrqcof(int n_data, int na, double *fmod, double *dfda, double *alpha,
       double *beta)
{	
    int    i, i1, j1;
    double c;

    for (j1=0; j1 < na; j1++) {
    for (i1=0; i1 < na; i1++) {
	c = 0.0;
	for (i=0; i < n_data; i++)
	    c += ( dfda[i1+i*na] * dfda[j1+i*na] );
	alpha[i1+j1*na] = c;
    }}
    for (i1=0; i1 < na; i1++) {
	c = 0.0;
	for (i=0; i < n_data; i++)
	    c += ( dfda[i1+i*na] * fmod[i] );
	beta[i1] = -c;
    }

    c = 0.0;				// compute the chisq;
    for (i=0; i < n_data; i++)
	c += ( fmod[i] * fmod[i] );
    return c;
}

#define maxlam  1.E+50
#define minlam  1.E-50
#define maxchi  1.E+12

int nlinfit(int n_data, void (*func)(int, void*, double*, double*, double*),
	    void *data, int na, double *a, double *da, double *chisq,
	    double tol, int maxloop, int verb)
{
    int     i, nloop, n_dof;
    double  alambda, _chisq, ochisq, dchisq=0.0;
    double *buf, *covar, *alpha, *fmod, *dfda, *atry, *beta, *da2;

    alambda = 0.001;
    n_dof   = n_data - na;
    buf     = malloc((2*na*na + n_data*(na+1) + 3*na)*sizeof(double));
    covar   = buf;			// size: na x na
    alpha   = covar + na*na;		// size: na x na
    fmod    = alpha + na*na;		// size: n_data
    dfda    = fmod  + n_data;		// size: na x n_data
    atry    = dfda  + n_data*na;	// size: na
    beta    = atry  + na;		// size: na
    da2     = beta  + na;		// size: na
    if (buf == NULL)
	pstop("!!! nlinfit: not enough memory.\n");

// Initial alpha, beta, and chisq.
    func(n_data, data, a, fmod, dfda);
    _chisq = mrqcof(n_data, na, fmod, dfda, alpha, beta);
    ochisq = _chisq;
    memcpy(covar, alpha, na*na*sizeof(double));
    memcpy(da2, beta, na*sizeof(double));

    for (nloop=0; nloop < maxloop; nloop++) {
        for (i=0; i < na; i++)
	    covar[i+i*na] = alpha[i+i*na] * (1.0+alambda);

// |da2> = (alpha*(1+alambda))^{-1}|da2>
	if (mx_inv(na, covar, 1, da2) != 0) {
	    if (verb >= 2)
		printf("!!! Lapack dsysv failed.\n");
	    free(buf);
	    return 1;
	}
        for (i=0; i < na; i++)
	    atry[i] = a[i] + da2[i];

	func(n_data, data, atry, fmod, dfda);
        _chisq = mrqcof(n_data, na, fmod, dfda, covar, da2);
	if (_chisq < 0.0 || _chisq > maxchi) {
	    if (verb >= 2)
		printf("!!! nlinfit: abnormal chisq: %E\n", _chisq);
	    free(buf);
	    return 1;
	}
	else if (_chisq < ochisq) {	// accept new alpha, beta, a, chisq.
	    if (alambda > minlam)
		alambda *= 0.1;
	    memcpy(alpha, covar, na*na*sizeof(double));
	    memcpy(beta, da2, na*sizeof(double));
	    memcpy(a, atry, na*sizeof(double));
	    dchisq = (ochisq - _chisq)/_chisq;
	    ochisq = _chisq;
	    if (dchisq <= tol) break;
	}
	else {				// reject the trial.
	    if (alambda < maxlam)
		alambda *= 10.0;
	    memcpy(covar, alpha, na*na*sizeof(double));
	    memcpy(da2, beta, na*sizeof(double));
	    _chisq = ochisq;
	}
    }
    if (nloop >= maxloop && dchisq > tol) {
	if (verb >= 2)
	    printf("!!! nlinfit: cannot converge in maxloop=%d, err=%E\n",
		    maxloop, dchisq);
	free(buf);
	return 1;
    }

// Compute:  covar = (alpha*(1+alambda))^{-1}, in order to get |da>.
    if (da != NULL) {
	memset(covar, 0, na*na*sizeof(double));
	for (i=0; i < na; i++) {
	    covar[i+i*na] = 1.0;
	    alpha[i+i*na] = alpha[i+i*na] * (1.0+alambda);
	}
	if (mx_inv(na, alpha, na, covar) != 0) {
	    if (verb >= 2)
		printf("Lapack dsysv failed.\n");
	    free(buf);
	    return 1;
	}
	*chisq = _chisq / ((double)n_dof);
	for (i=0; i < na; i++) {
	    if (covar[i+i*na] < 0.0) {
		if (verb >= 2)
		    printf("Negative covar matrix at %d\n", i);
		return 1;
	    }
	    da[i] = sqrt(covar[i+i*na] * (*chisq));
	}
    }
    free(buf);

    return 0;
}
