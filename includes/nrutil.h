/*
**
** 	Utility Routines
** 	from Numerical Recipies in C, W. Press, et al, Cambridge 1988.
**
**	Second edition
*/

#ifndef _NRUTIL_

#define _NRUTIL_ 1

#include "VAX_Alloc.h"

#ifdef __cplusplus

extern "C" {

#endif /* _cplusplus */

#ifndef __GNUC__

// #pragma unused( farg1, farg2, darg1, darg2, larg1, larg2, iarg1, iarg2)
static float  farg1 , farg2 ;
static double darg1 , darg2 ;
static long   larg1 , larg2 ;
static int    iarg1 , iarg2 ;

#else

static float  farg1 __attribute__ ((unused)), farg2 __attribute__ ((unused));
static double darg1 __attribute__ ((unused)), darg2 __attribute__ ((unused));
static long   larg1 __attribute__ ((unused)), larg2 __attribute__ ((unused));
static int    iarg1 __attribute__ ((unused)), iarg2 __attribute__ ((unused));

#endif

#define SQR(a)  ((farg1=(a)) == 0.0 ? 0.0 : farg1*farg1)
#define DSQR(a) ((darg1=(a)) == 0.0 ? 0.0 : darg1*darg1)

#define DMAX(a,b) (darg1=(a),darg2=(b),(darg1) > (darg2) ? (darg1) : (darg2))
#define DMIN(a,b) (darg1=(a),darg2=(b),(darg1) < (darg2) ? (darg1) : (darg2))

#define FMAX(a,b) (farg1=(a),farg2=(b),(farg1) > (farg2) ? (farg1) : (farg2))
#define FMIN(a,b) (farg1=(a),farg2=(b),(farg1) < (farg2) ? (farg1) : (farg2))

#define LMAX(a,b) (larg1=(a),larg2=(b),(larg1) > (larg2) ? (larg1) : (larg2))
#define LMIN(a,b) (larg1=(a),larg2=(b),(larg1) < (larg2) ? (larg1) : (larg2))

#define IMAX(a,b) (iarg1=(a),iarg2=(b),(iarg1) > (iarg2) ? (iarg1) : (iarg2))
#define IMIN(a,b) (iarg1=(a),iarg2=(b),(iarg1) < (iarg2) ? (iarg1) : (iarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/*
**	P R O T O T Y P E S
*/

void nrerror(const char *error_text);
void nrinfo(const char *error_text);

float *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
double *dvector0(long nl, long nh);

float **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix0(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);


#ifdef __cplusplus

}

#endif /* _cplusplus */


#endif
