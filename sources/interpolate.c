/*
** TokaMac v2.0
**
** interpolate.c
**
**
** File:		interpolate.c
** Date:		March 12, 1993
**
** Revisions:
**
**		August 3, 1993		Added interpolate_int
**		August 8, 1993		Fixed small typo
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <math.h>
#include "psigrid.h"
#include "interpolate.h"

#define PT_FLOOR(x, xmin, dx)	(int) floor( ((x)-(xmin))/(dx) )

double limit(double x, double xmin, double xmax) {
    if (x < xmin) {
        return xmin;
    } else if (x > xmax) {
        return xmax;
    } else {
        return x;
    }
}

/***************************************************************
 * interpolate
 *
 * This function interpolates/extrapolates an array A in order to
 * return its value at the requested point, (x, z).
 *
 * We use a simple bilinear interpolation.
 *
 ***************************************************************/
double        interpolate(PSIGRID * pg, double **A, double x, double z)
{
	int           ix, iz;
	double        xm, zm;
	double        t, u;

    xm = limit(x, pg->Xmin, pg->Xmax);
    zm = limit(z, pg->Zmin, pg->Zmax);

	ix = PT_FLOOR(xm, pg->Xmin, pg->dx);
	iz = PT_FLOOR(zm, pg->Zmin, pg->dz);
	t = (xm - pg->X[ix]) / pg->dx;
	u = (zm - pg->Z[iz]) / pg->dz;

	return ((1.0 - t) * (1.0 - u) * A[ix][iz] + t * (1.0 - u) * A[ix + 1][iz] +
			t * u * A[ix + 1][iz + 1] + (1.0 - t) * u * A[ix][iz + 1]);
}

/***************************************************************
 * interpolate_int
 *
 * This function interpolates/extrapolates an array A in order to
 * return its value at the requested point, (x, z).
 *
 * We use a simple bilinear interpolation.
 *
 ***************************************************************/
double        interpolate_int(PSIGRID * pg, int **A, double x, double z)
{
	int           ix, iz;
	double        xm, zm;
	double        t, u;

    xm = limit(x, pg->Xmin, pg->Xmax);
    zm = limit(z, pg->Zmin, pg->Zmax);

	ix = PT_FLOOR(xm, pg->Xmin, pg->dx);
	iz = PT_FLOOR(zm, pg->Zmin, pg->dz);
	t = (xm - pg->X[ix]) / pg->dx;
	u = (zm - pg->Z[iz]) / pg->dz;

	return ((1.0 - t) * (1.0 - u) * A[ix][iz] + t * (1.0 - u) * A[ix + 1][iz] +
			t * u * A[ix + 1][iz + 1] + (1.0 - t) * u * A[ix][iz + 1]);
}
