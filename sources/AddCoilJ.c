/*
** TokaMac v2.0
**
** AddCoilJ.c
**
** This file contains routines which add any coil currents
** located within the computational domain.
**
**
** File:		AddCoilJ.c
** Date:		March 20, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <math.h>
#include "psigrid.h"
#include "coil.h"
#include "AddCoilJ.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

#define PT_FLOOR(x, xmin, xmax, nsize)	floor((nsize)*((x)-(xmin))/((xmax)-(xmin)))

/*
**	AddSubCoilJ
**
**
*/

void          AddSubCoilJ(PSIGRID * pg, COIL * c)
{
	int           isc, nmax;
	int           ixc, izc;		/* the closest grid points for the subcoil location */
	SUBCOIL      *sc;
	double        xc, zc;
	double        Jcoil;		/* the current density of the subcoil */
	double        dxn, dzn;
	double      **J;

	J = pg->Current;
	nmax = pg->Nsize;

	for (isc = 0; isc < c->NumSubCoils; isc++) {
		sc = c->SubCoils[isc];
		xc = sc->X;
		zc = sc->Z;
		Jcoil = c->CoilCurrent * sc->CurrentFraction / (pg->dx * pg->dz);
		ixc = (int) PT_FLOOR(xc, pg->Xmin, pg->Xmax, nmax);
		izc = (int) PT_FLOOR(zc, pg->Zmin, pg->Zmax, nmax);
		if ((ixc > 0) && (ixc < nmax) && (izc > 0) && (izc < nmax)) {
			dxn = (xc - pg->X[ixc]) / pg->dx;
			dzn = (zc - pg->Z[izc]) / pg->dz;
			J[ixc][izc] += Jcoil * (1.0 - dxn) * (1.0 - dzn);
			J[ixc + 1][izc] += Jcoil * dxn * (1.0 - dzn);
			J[ixc][izc + 1] += Jcoil * (1.0 - dxn) * dzn;
			J[ixc + 1][izc + 1] += Jcoil * dxn * dzn;
			if (pg->Symmetric) {
				J[ixc][nmax - izc] = J[ixc][izc];
				J[ixc + 1][nmax - izc] = J[ixc + 1][izc];
				J[ixc][nmax - izc - 1] = J[ixc][izc + 1];
				J[ixc + 1][nmax - izc - 1] = J[ixc + 1][izc + 1];
			}
		}
	}
}

void          AddCoilJ(TOKAMAK * td)
{
	int           i;

	for (i = 0; i < td->NumCoils; i++)
		if (td->Coils[i]->Enabled == Coil_On)
			AddSubCoilJ(td->PsiGrid, td->Coils[i]);
}
