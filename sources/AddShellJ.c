/*
** TokaMac v2.0
**
** AddShellJ.c
**
** This file contains routines which add any shell currents
** located within the computational domain.
**
** NOTE:  All shells should be within the computational domain.
**        But, everything should work when they are outside the domain.
**
**
** File:		AddShellJ.c
** Date:		August 3, 1993
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <math.h>
#include "psigrid.h"
#include "shell.h"
#include "AddShellJ.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

#define PT_FLOOR(x, xmin, xmax, nsize)	floor((nsize)*((x)-(xmin))/((xmax)-(xmin)))

/*
**	AddSubShellJ
**
**
*/

void          AddSubShellJ(PSIGRID * pg, SHELL * c)
{
	int           isc, nmax;
	int           ixc, izc;		/* the closest grid points for the subshell location */
	SUBSHELL     *sc;
	double        xc, zc;
	double        Jc;			/* the current density of the subshell */
	double        dxn, dzn;
	double      **J;

	J = pg->Current;
	nmax = pg->Nsize;

	for (isc = 0; isc < c->NumSubShells; isc++) {
		sc = c->SubShells[isc];
		xc = sc->X;
		zc = sc->Z;
		Jc = sc->Current / (pg->dx * pg->dz);
		ixc = (int) PT_FLOOR(xc, pg->Xmin, pg->Xmax, nmax);
		izc = (int) PT_FLOOR(zc, pg->Zmin, pg->Zmax, nmax);
		if ((ixc > 0) && (ixc < nmax) && (izc > 0) && (izc < nmax)) {
			dxn = (xc - pg->X[ixc]) / pg->dx;
			dzn = (zc - pg->Z[izc]) / pg->dz;
			J[ixc][izc] += Jc * (1.0 - dxn) * (1.0 - dzn);
			J[ixc + 1][izc] += Jc * dxn * (1.0 - dzn);
			J[ixc][izc + 1] += Jc * (1.0 - dxn) * dzn;
			J[ixc + 1][izc + 1] += Jc * dxn * dzn;
			if (pg->Symmetric) {
				J[ixc][nmax - izc] = J[ixc][izc];
				J[ixc + 1][nmax - izc] = J[ixc + 1][izc];
				J[ixc][nmax - izc - 1] = J[ixc][izc + 1];
				J[ixc + 1][nmax - izc - 1] = J[ixc + 1][izc + 1];
			}
		}
	}
}

void          AddShellJ(TOKAMAK * td)
{
	int           i;

	for (i = 0; i < td->NumShells; i++)
		AddSubShellJ(td->PsiGrid, td->Shells[i]);
}
