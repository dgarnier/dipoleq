/*
** TokaMac v2.0
**
** J_Std.c
**
**
**
** File:		J_Std.c
** Date:		March 22, 1993
**
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <math.h>
#include "nrutil.h"
#include "plasma.h"
#include "psigrid.h"
#include "tokamak.h"
#include "J_Std.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

/***
* 	Remember: the plasma current will be positive.
*	We're defining it exactly as in Eq. 7 of Johnson, et al.
*
* In other words,
*
*  		Æ*Psi = 2piX J(x,z)
*
* and
*
*  		J(x,z) = -2pi*(X ¶p/¶Psi + ((R0*B0)^2/2X)*¶g^2/¶Psi)
*
***/

void          J_Std(TOKAMAK * td, double **J, double p0, double g0)
{
	PSIGRID      *pg;
	PLASMA       *pl;
	int           ix, iz, nmax;
	double       *X, **Psi;
	double        PsiAxis, DelPsi;
	double        Ptemp, Gtemp;
	double        dPress;
	double        dG2;
	double        PsiX;

	pg = td->PsiGrid;
	X = pg->X;
	Psi = pg->Psi;
	PsiAxis = pg->PsiAxis;
	DelPsi = pg->PsiLim - PsiAxis;
	pl = td->Plasma;
	nmax = pg->Nsize;

	for (ix = 0; ix <= nmax; ix++)
		J[ix][0] = J[ix][nmax] = 0.0;
	for (iz = 0; iz <= nmax; iz++)
		J[0][iz] = J[nmax][iz] = 0.0;

	/* Calculate the current... */
	for (ix = 1; ix < nmax; ix++) {
		Gtemp = -PI * DSQR(pl->B0R0) / X[ix];
		Ptemp = -TWOPI * X[ix];
		for (iz = 1; iz < nmax; iz++) {
			if (pg->IsPlasma[ix][iz]) {
				PsiX = (Psi[ix][iz] - PsiAxis) / DelPsi;
				dPress = p0 * pow(1.0 - PsiX, pl->StndP - 1.0);
				dG2 = g0 * pow(1.0 - PsiX, pl->StndG - 1.0);
				J[ix][iz] = Ptemp * dPress + Gtemp * dG2;
			} else
				J[ix][iz] = 0.0;
		}
	}
}

double        J_Std_Loc(TOKAMAK * td, int ix, int iz, double p0, double g0)
{
	PSIGRID      *pg;
	PLASMA       *pl;
	double       *X, **Psi;
	double        PsiAxis, DelPsi;
	double        Ptemp, Gtemp;
	double        dPress;
	double        dG2;
	double        PsiX;
	double        J = 0.0;

	pg = td->PsiGrid;
	X = pg->X;
	Psi = pg->Psi;
	PsiAxis = pg->PsiAxis;
	DelPsi = pg->PsiLim - PsiAxis;
	pl = td->Plasma;

	if (pg->IsPlasma[ix][iz]) {
		Gtemp = -PI * DSQR(pl->B0R0) / X[ix];
		Ptemp = -TWOPI * X[ix];
		PsiX = (Psi[ix][iz] - PsiAxis) / DelPsi;
		dPress = p0 * pow(1.0 - PsiX, pl->StndP - 1.0);
		dG2 = g0 * pow(1.0 - PsiX, pl->StndG - 1.0);
		J = (Ptemp * dPress + Gtemp * dG2);
	}
	return J;
}
