/*
** TokaMac v2.0
**
** meas_plasmacur.c
**
** Source file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_plasmacur_Green(TOKAMAK *td, MEAS *m)
**
**				meas_plasmacur_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_plasmacur_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_plasmacur.c
** Date:		March 24, 1993
**
** Revisions:
**
**		August 6, 1993		Added xxx_Now
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdio.h>
#include <math.h>
#include "green.h"
#include "nrutil.h"
#include "measurement.h"
#include "psigrid.h"
#include "plasma.h"
#include "tokamak.h"
#include "FindJ.h"
#include "meas_plasmacur.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

extern double ***meas_dJdy;		/* common pointer to storage */

/*
**	meas_plasmacur_Fit
**
**
*/
void          meas_plasmacur_Fit(TOKAMAK * td, MEAS * m)
{
	double        sum = 0.0;
	int           nmax;
	int           ix, iz, **ip;
	double      **J;

	nmax = td->PsiGrid->Nsize;
	J = td->PsiGrid->Current;
	ip = td->PsiGrid->IsPlasma;

	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (ip[ix][iz])
				sum += J[ix][iz];
	sum *= td->PsiGrid->dx * td->PsiGrid->dz;

	m->Fit = sum / MU0;
}

/*
**	meas_plasmacur_Now
**
**
*/
void          meas_plasmacur_Now(TOKAMAK * td, MEAS * m)
{
	double        sum = 0.0;
	int           nmax;
	int           ix, iz, **ip;

	nmax = td->PsiGrid->Nsize;
	ip = td->PsiGrid->IsPlasma;

	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			sum += FindJ_Loc(td, ix, iz);
	sum *= td->PsiGrid->dx * td->PsiGrid->dz;

	m->Now = sum / MU0;
}

/*
**	meas_plasmacur_L
**
**	Note:	The L array has dimensions [1..NumUnkns].
**
*/
void          meas_plasmacur_L(TOKAMAK * td, MEAS *dummy, double *L)
{
	int           iu, iup;
	int           ix, iz, nmax;

	nmax = td->PsiGrid->Nsize;
	iup = td->NumUnkns - td->NumCoils;	/* iup is the number of plasma unknowns */

	for (iu = 1; iu <= td->NumUnkns; iu++)
		L[iu] = 0.0;

	/* P L A S M A   U N K N O W N S */
	for (iu = 1; iu <= iup; iu++) {
		for (ix = 1; ix < nmax; ix++)
			for (iz = 1; iz < nmax; iz++)
				L[iu] += meas_dJdy[iu][ix][iz];
		L[iu] *= td->PsiGrid->dx * td->PsiGrid->dz;
		L[iu] /= MU0;
	}
}
