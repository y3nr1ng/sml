#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "calb.h"

#define EPS  1.E-16

/*-------------------------------------------------------------------------
 *
 *  LU decomposition:
 *
 *  A = L1 (L1^{-1} A0)
 *    = L1 A1 P1
 *    = L1 L2 (L2^{-1} A1) P1
 *    = L1 L2 A2 P2 P1
 *    = ....
 *    = (L1 L2 ... LN) AN (PN ... P2 P1)
 *
 *  where
 *
 *    L = L1 L2 ... LN
 *    U = AN = LN^{-1} ... L2^{-1} L1^{-1} A0
 *    P = PN ... P2 P1
 *
 *        | 1                   |              | 1                    |
 *    Li: |      1              |     Li^{-1}: |       1              |
 *        |   l_{i+1,i}   1     |              |   -l_{i+1,i}   1     |
 *        |   l_{i+2,i}       1 |              |   -l_{i+2,i}       1 |
 *
 *        l_{i,j} = a^{(j-1)}_{i,j} / a^{(j-1)}_{j,j}
 *
 *
 *        | 1       u_{i-2,i}   |              | 1       -u_{i-2,i}   |
 *    Ui: |      1  u_{i-1,i}   |     Ui^{-1}: |       1 -u_{i-1,i}   |
 *        |           1         |              |            1         |
 *        |                   1 |              |                    1 |
 *
 *        | U_{1,1}                |
 *    Ud: |      U_{2,2}           |
 *        |           U_{3,3}      |
 *        |                U_{N,N} |
 *
 *        u_{i,j} = U_{i,j} / U_{i,i}
 *
 *    Pi: the permutation matrix for partial pivoting of A, Pi = Pi^{-1}
 *
 *                                  | 0  1  0  0 |
 *        for permutation 1 and 2:  | 1  0  0  0 |
 *                                  | 0  0  1  0 |
 *                                  | 0  0  0  1 |
 *
 *  A = L U P
 *  L = L1 L2 ... LN          L^{-1} = LN^{-1} ... L2^{-1} L1^{-1}
 *  U = Ud UN ... U2 U1       U^{-1} = U1^{-1} U2^{-1} ... UN^{-1} Ud^{-1}
 *  P = PN ... P2 P1          P^{-1} = P1 P2 ... PN
 *
 *------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
 *  Partial pivoting for matrix A.
 */
static void get_ipiv(int N, int irow, double *A, int *ipiv)
{
    int    i, j1, j2;
    double c=0.0;

    for (i=irow; i < N; i++)
	c = c + fabs(A[irow+i*N]);
    c /= ((double)(N-irow));

    ipiv[irow] = irow;
    for (i=irow; i < N; i++) {
	if (fabs(A[irow+i*N]/c) > EPS) {
	    ipiv[irow] = i;
	    break;
	}
    }

    if (ipiv[irow] != irow) {
	j1 = irow;
	j2 = ipiv[irow];
	for (i=0; i < N; i++) {
	    c = A[i+j1*N];
	    A[i+j1*N] = A[i+j2*N];
	    A[i+j2*N] = c;
	}
    }
}

/*-------------------------------------------------------------------------
 *  Interchange the elements of B to restore the correct ordering.
 */
static void ipiv_rhs(int N, int *ipiv, int Nrhs, double *X)
{
    int     i, j, k;
    double  c;

    for (i=N-1; i >= 0; i--) {
	j = ipiv[i];
	if (i == j) continue;

	for (k=0; k < Nrhs; k++) {
	    c = X[i+k*N];
	    X[i+k*N] = X[j+k*N];
	    X[j+k*N] = c;
	}
    }
}

/*-------------------------------------------------------------------------
 *  Perform the LU decomposition with partial pivoting for matrix A.
 */
static void mxLU(int N, double *A, int *ipiv, double *LU)
{
    int     i, j, k;
    double *L, lij;

    L = (double *)memalloc(sizeof(double)*N*N, NULL);
    memcpy(LU, A, sizeof(double)*N*N);

    for (j=0; j < N; j++) {
	get_ipiv(N, j, LU, ipiv);

	for (i=j+1; i < N; i++) {
	    lij = LU[i+j*N] / LU[j+j*N];
	    L[i+j*N] = lij;
	    for (k=0; k < N; k++)
		LU[i+k*N] -= ( lij * LU[j+k*N] );
	}
    }
    for (j=0; j < N; j++) {
	for (i=j+1; i < N; i++)
	    LU[i+j*N] = L[i+j*N];
    }
    free(L);
}

/*-------------------------------------------------------------------------
 *  Perform the L^{-1} times a vector.
 */
static void mxLinv_mul(int N, double *LU, double *x)
{
    int  i, j;

    for (j=0; j < N; j++) {
	for (i=j+1; i < N; i++)
	    x[i] -= (LU[i+j*N] * x[j]);
    }
}

/*-------------------------------------------------------------------------
 *  Perform the U^{-1} times a vector.
 */
static void mxUinv_mul(int N, double *LU, double *x)
{
    int  i, j;

    for (i=0; i < N; i++)
	x[i] /= LU[i+i*N];

    for (j=N-1; j >= 0; j--) {
	for (i=0; i < j; i++)
	    x[i] -= (LU[i+j*N]/LU[i+i*N] * x[j]);
    }
}

/*-------------------------------------------------------------------------
 *  Solve the linear system:  A X = B, with LU decomposition.
 */
void dsysv(int N, double *A, int Nrhs, double *B, double *X)
{
    int    *ipiv, i;
    double *LU;

    ipiv = (int    *)memalloc(sizeof(int)*N, NULL);
    LU   = (double *)memalloc(sizeof(double)*N*N, NULL);
    if (X != B)
	memcpy(X, B, sizeof(double)*N*Nrhs);
    mxLU(N, A, ipiv, LU);

    for (i=0; i < Nrhs; i++) {
	mxLinv_mul(N, LU, X+i*N);
	mxUinv_mul(N, LU, X+i*N);
    }
    ipiv_rhs(N, ipiv, Nrhs, X);
}
