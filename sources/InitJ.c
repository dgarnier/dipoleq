/*
** TokaMac v2.0
**
** InitJ.c
**
** This file contains routines which initializes the plasma
** current using info within the Plasma data structure.
**
**
** File:		InitJ.c
** Date:		March 20, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdio.h>
#include <math.h>
#include "plasma.h"
#include "psigrid.h"
#include "InitJ.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

#define SQUARE(x)	((x) * (x))

extern FILE  *LogFile;

/*
** InitJ
**
** This method computes the initial current density that can be used
** when beginning a reconstruction without knowing a better solution.
**
** We use the simple formula...J(x,z) = (n/x)*(1 - r^2/d^2)
** where
**  n = normalization
**  r = minor radial coordinate
**  d = minor radius
**
** **Remember: the plasma current will be **POSITIVE** & we're defining
** it exactly as in Eq. 7 of Johnson, et al.
**
** In other words,
**  Æ*Psi = 2piX J(x,z)
** and
**  J(x,z) = -2pi*(X ¶p/¶Psi + ((R0*B0)^2/2X)*¶g^2/¶Psi)
**
**
*/

void          InitJ(PSIGRID * pg, PLASMA * pl)
{
	int           ix, iz, nmax;
	double        r2, d2;
	double        sum = 0.0;
	double      **J;

	J = pg->Current;
	d2 = SQUARE(pl->a0);
	nmax = pg->Nsize;

	printf("INFO:	Initializing current to %g (A).\n", pl->Ip0);
	fprintf(LogFile, "INFO:	Initializing current to %g (A).\n", pl->Ip0);

	for (ix = 0; ix <= nmax; ix++)
		for (iz = 0; iz <= nmax; iz++) {
			r2 = SQUARE(pg->X[ix] - pl->R0) + SQUARE(pg->Z[iz] - pl->Z0);
			J[ix][iz] = ((ix != 0) && (ix != nmax) && (iz != 0) && (iz != nmax) && (r2 < d2)) ?
				pg->X[ix] * (1.0 - r2 / d2) : 0.0;
			sum += J[ix][iz];
		}

	sum *= pg->dx * pg->dz;
	sum = pl->Ip0 / sum;
	sum *= MU0;

	/* Re-Scale the current density...*/

	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			J[ix][iz] *= sum;
}
