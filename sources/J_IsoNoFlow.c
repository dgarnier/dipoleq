/*
** TokaMac v2.0
**
** J_IsoNoFlow.c
**
**
**
** File:		J_IsoNoFlow.c
** Date:		March 22, 1993
**
** Revisions:
**
**		August 6, 1993		Added local current density computation
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <math.h>
#include "nrutil.h"
#include "plasma.h"
#include "psigrid.h"
#include "tokamak.h"
#include "fpoly.h"
#include "J_IsoNoFlow.h"

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

void          J_IsoNoFlow(TOKAMAK * td, double **J, double *Pp, double *G2p)
{
	PSIGRID      *pg;
	PLASMA       *pl;
	int           ix, iz, nmax;
	int           PpTerms, G2pTerms;
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
	PpTerms = pl->PpTerms;
	G2pTerms = pl->G2pTerms;
	nmax = pg->Nsize;

	for (ix = 0; ix <= nmax; ix++)
		J[ix][0] = J[ix][nmax] = 0.0;
	for (iz = 0; iz <= nmax; iz++)
		J[0][iz] = J[nmax][iz] = 0.0;

	/* Calculate the current... */
	for (ix = 1; ix < nmax; ix++) {
		Gtemp = -PI * DSQR(pl->B0R0) / X[ix];
		Ptemp = -TWOPI * X[ix];
		for (iz = 1; iz < pg->Nsize; iz++) {
			if (pg->IsPlasma[ix][iz]) {
				PsiX = (Psi[ix][iz] - PsiAxis) / DelPsi;
				dPress = fpoly(Pp, PsiX, PpTerms);
				dG2 = fpoly(G2p, PsiX, G2pTerms);
				J[ix][iz] = Ptemp * dPress + Gtemp * dG2;
			} else
				J[ix][iz] = 0.0;
		}
	}
}

double        J_IsoNoFlow_Loc(TOKAMAK * td, int ix, int iz, double *Pp, double *G2p)
{
	PSIGRID      *pg;
	PLASMA       *pl;
	int           PpTerms, G2pTerms;
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
	PpTerms = pl->PpTerms;
	G2pTerms = pl->G2pTerms;

	if (pg->IsPlasma[ix][iz]) {
		Gtemp = -PI * DSQR(pl->B0R0) / X[ix];
		Ptemp = -TWOPI * X[ix];
		PsiX = (Psi[ix][iz] - PsiAxis) / DelPsi;
		dPress = fpoly(Pp, PsiX, PpTerms);
		dG2 = fpoly(G2p, PsiX, G2pTerms);
		J = (Ptemp * dPress + Gtemp * dG2);
	}
	return J;
}
