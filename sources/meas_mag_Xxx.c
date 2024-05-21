/*
** TokaMac v2.0
**
** meas_mag_Xxx.c
**
** Source file for local poloidal field measurements.
**
** Every magnetic measurement has a common subroutine:
**
**				meas_mag_Fit(TOKAMAK *td, MEAS *m)
**
**
** File:		meas_mag_Xxx.c
** Date:		March 24, 1993
**
** Revisions:
**
**		August 6, 1993		Added perfectly conducting shells
**		August 6, 1993		Added xxx_Now
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <math.h>
#include "measurement.h"
#include "coil.h"
#include "shell.h"
#include "psigrid.h"
#include "tokamak.h"
#include "FindJ.h"
#include "meas_mag_Xxx.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307
#define DEG_RAD		1.74532925e-2

#define TRUE		1
#define FALSE		0

extern double ***meas_dJdy;		/* common pointer to storage */

/*
**	meas_mag_Fit
**
*/
void          meas_mag_Fit(TOKAMAK * td, MEAS * m)
{
	double        sum = 0.0;
	int           nmax;
	int           ix, iz, i, is, **ip;
	double      **PlGrn, *CGrn, *SGrn, **J;
	SHELL        *shell;
	COIL		*c;

	nmax = td->PsiGrid->Nsize;
	J = td->PsiGrid->Current;
	ip = td->PsiGrid->IsPlasma;
	PlGrn = m->PlasmaGreen;
	SGrn = m->ShellGreen;
	CGrn = m->CoilGreen;

	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (ip[ix][iz])
				sum += J[ix][iz] * PlGrn[ix][iz];
	sum = sum * td->PsiGrid->dx * td->PsiGrid->dz;

	for (i = 0; i < td->NumCoils; i++) {
		c = td->Coils[i];
		if (c->Enabled)
			sum += CGrn[i] * c->CoilCurrent;
	}

	for (i = 0; i < td->NumShells; i++) {
		shell = td->Shells[i];
		for (is = 0; is < shell->NumSubShells; is++)
			sum += SGrn[i] * shell->SubShells[i]->Current;
	}

	m->Fit = sum;
}

/*
**	meas_mag_Now
**
*/
void          meas_mag_Now(TOKAMAK * td, MEAS * m)
{
	double        sum = 0.0;
	int           nmax;
	int           ix, iz, i, is, **ip;
	double      **PlGrn, *CGrn, *SGrn;
	SHELL        *shell;
	COIL		*c;

	nmax = td->PsiGrid->Nsize;
	ip = td->PsiGrid->IsPlasma;
	PlGrn = m->PlasmaGreen;
	SGrn = m->ShellGreen;
	CGrn = m->CoilGreen;

	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (ip[ix][iz])
				sum += FindJ_Loc(td, ix, iz) * PlGrn[ix][iz];
	sum = sum * td->PsiGrid->dx * td->PsiGrid->dz;

	for (i = 0; i < td->NumCoils; i++) {
		c = td->Coils[i];
		if (c->Enabled)
			sum += CGrn[i] * c->CoilCurrent;
	}

	for (i = 0; i < td->NumShells; i++) {
		shell = td->Shells[i];
		for (is = 0; is < shell->NumSubShells; is++)
			sum += SGrn[i] * shell->SubShells[i]->Current;
	}

	m->Now = sum;
}

/*
**	meas_mag_L
**
**	Note:	The L array has dimensions [1..NumUnkns].
**
*/
void          meas_mag_L(TOKAMAK * td, MEAS * m, double *L)
{
	int           iu, iup;
	int           ix, iz, nmax;
	double      **PlGrn, *CGrn;

	nmax = td->PsiGrid->Nsize;
	iup = td->NumUnkns - td->NumCoils;	/* iup is the number of plasma unknowns */

	PlGrn = m->PlasmaGreen;
	CGrn = m->CoilGreen;

	for (iu = 1; iu <= td->NumUnkns; iu++)
		L[iu] = 0.0;

	/* P L A S M A   U N K N O W N S */
	for (iu = 1; iu <= iup; iu++) {
		for (ix = 1; ix < nmax; ix++)
			for (iz = 1; iz < nmax; iz++)
				L[iu] += PlGrn[ix][iz] * meas_dJdy[iu][ix][iz];
		L[iu] *= td->PsiGrid->dx * td->PsiGrid->dz;
	}

	/* C O I L   U N K N O W N S */
	for (iu = iup + 1; iu <= td->NumUnkns; iu++)
		L[iu] = CGrn[iu - iup - 1];
}
