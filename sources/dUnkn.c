/*
**	dUnkn.c
**
**	Computes a random vector representing the unknowns of
**	a least-squares data fitting problem.
**
**	October 18, 1993
**
**	M. E. Mauel -- Dept. of Applied Physics, Columbia University
**
**
*/

#include <math.h>
#include "nrutil.h"
#include "dUnkn.h"

long          GLOBAL_RAN3_SEED;	/* initilize to some negative number */

/*
**	r a n 1
**
**	The recommended random number generator in Numerical Recipies
**
*/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double        ran1(long *idum)
{
	int           j;
	long          k;
	static long   iy = 0;
	static long   iv[NTAB];
	double        temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1)
			*idum = 1;
		else
			*idum = -(*idum);
		for (j = NTAB + 7; j >= 0; j--) {
			k = (*idum) / IQ;
			*idum = IA * (*idum - k * IQ) - IR * k;
			if (*idum < 0)
				*idum += IM;
			if (j < NTAB)
				iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ;
	*idum = IA * (*idum - k * IQ) - IR * k;
	if (*idum < 0)
		*idum += IM;
	j = iy / NDIV;
	iy = iv[j];
	iv[j] = *idum;
	if ((temp = AM * iy) > RNMX)
		return RNMX;
	else
		return temp;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

/*
**	r a n 3
**
**	From Numerical Recipies in C, p. 283
**	A random number from 0.0 to 1.0.
**
**	Set idum to any negative number to initialize or reinitilize the sequence.
**
*/

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double        ran3(long *idum)
{
	static int    inext, inextp;
	static long   ma[56];
	static int    iff = 0;
	long          mj, mk;
	int           i, ii, k;

	if (*idum < 0 || iff == 0) {
		iff = 1;
		mj = MSEED - (*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;
		for (i = 1; i <= 54; i++) {
			ii = (21 * i) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if (mk < MZ)
				mk += MBIG;
			mj = ma[ii];
		}
		for (k = 1; k <= 4; k++)
			for (i = 1; i <= 55; i++) {
				ma[i] -= ma[1 + (i + 30) % 55];
				if (ma[i] < MZ)
					ma[i] += MBIG;
			}
		inext = 0;
		inextp = 31;
		*idum = 1;
	}
	if (++inext == 56)
		inext = 1;
	if (++inextp == 56)
		inextp = 1;
	mj = ma[inext] - ma[inextp];
	if (mj < MZ)
		mj += MBIG;
	ma[inext] = mj;
	return (double) mj *FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

/*
** g a s d e v
**
** Normal Gaussian deviate with zero mean and unit standard deviation
**
*/

double        gasdev(long *idum)
{
	static int    iset = 0;
	static double gset;
	double        fac, rsq, v1, v2;

	if (iset == 0) {
		do {
			v1 = 2.0 * ran1(idum) - 1.0;
			v2 = 2.0 * ran1(idum) - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0 * log(rsq) / rsq);
		gset = v1 * fac;
		iset = 1;
		return v2 * fac;
	} else {
		iset = 0;
		return gset;
	}
}

/*
**	d U n k n
**
**	del is the bounding value of ChiSqr.
**	The length of an axis of the ellipsoid is sqrt(del)/w
**
*/

#define TOL		1.0e-6

void          dUnkn(double **V, double *w, double del, double *unkn0, double *unkn1, int nu)
{
	int           i, j;
	double       *rana, *du, *V2;
	double        sdel;
	double        wmin;

	wmin = w[1];
	for (i = 1; i <= nu; i++)
		if (wmin < w[i])
			wmin = w[i];
	wmin = TOL * nu * wmin;		/* let wmin be equal to "machine precision" times wmax */

	sdel = sqrt(del / nu);

	rana = dvector(1, nu);
	du = dvector0(1, nu);
	V2 = dvector0(1, nu);

	for (i = 1; i <= nu; i++)
		rana[i] = gasdev(&GLOBAL_RAN3_SEED);

	for (i = 1; i <= nu; i++)
		if (w[i] < wmin)
			rana[i] = 0.0;		/* ignore singular values near zero */
		else
			rana[i] = rana[i] * sdel / w[i];

	for (i = 1; i <= nu; i++)
		for (j = 1; j <= nu; j++)
			V2[i] += V[j][i] * V[j][i];	/* V2 contains the magnitude of each column V */

	for (i = 1; i <= nu; i++)
		for (j = 1; j <= nu; j++)
			du[i] += rana[j] * V[i][j] / V2[j];

	for (i = 1; i <= nu; i++)
		unkn1[i] = unkn0[i] + du[i];

	free_dvector(V2, 1, nu);
	free_dvector(du, 1, nu);
	free_dvector(rana, 1, nu);
}

#undef TOL
