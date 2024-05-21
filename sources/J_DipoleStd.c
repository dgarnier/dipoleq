/*
** Dipole v0.9
**
** J_DipoleStd.c
**
**
**
** File:		J_DipoleStd.c
** Date:		December 11, 1999
**
** Revisions:
**
**
**
**
** (c) D. Garnier -- Columbia University
*/

#include <math.h>
#include "nrutil.h"
#include "plasma.h"
#include "psigrid.h"
#include "tokamak.h"
#include "J_DipoleStd.h"

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

double		Pp_DipoleStd_Loc(TOKAMAK *td, double Psi)
{
	double		dT_dPsiX, Psi_0, theta;
	double		dPress;
	double      PsiDipole, PsiWall;
	PSIGRID     *pg;
	double      *Pp;

	pg = td->PsiGrid;
	PsiDipole = pg->PsiAxis;
	PsiWall   = pg->PsiAxis;
	Pp        = td->Plasma->Pp;
	Psi_0     = GetPsi(pg, Pp[0], Pp[1]);

	if (Psi <= Psi_0) {
		dT_dPsiX = (PI/2) / (Psi_0 - PsiDipole);
		theta = dT_dPsiX * (Psi - PsiDipole);
		dPress = Pp[2]*dT_dPsiX*2*sin(theta)*cos(theta);
	} else {
		dPress = Pp[2]/Psi_0*Pp[3]*pow( Psi/Psi_0 , Pp[3] - 1);
	}

	return dPress;
}

double		P_DipoleStd_Loc(TOKAMAK *td, double Psi)
{
	double 		theta;
	double		Psi_0;
	double		Press;
	double      PsiDipole, PsiWall;
	PSIGRID     *pg;
	double      *Pp;

	pg = td->PsiGrid;
	PsiDipole = pg->PsiAxis;
	PsiWall   = pg->PsiAxis;
	Pp        = td->Plasma->Pp;
	Psi_0     = GetPsi(pg, Pp[0], Pp[1]);

	if (Psi <= Psi_0) {
		theta = (PI/2) * (Psi - PsiDipole) / (Psi_0 - PsiDipole);
		Press = Pp[2]*sin(theta)*sin(theta);
	} else {
		Press = Pp[2]*pow( Psi/Psi_0 , Pp[3]);
	}

	return Press;

}

double		G2p_DipoleStd_Loc(TOKAMAK *td, double Psi)
{
	return 0.0;
}

double		G2_DipoleStd_Loc(TOKAMAK *td, double Psi)
{
	return 1.0;
}



void          J_DipoleStd(TOKAMAK * td, double **J)
{
	PSIGRID      *pg;
	PLASMA		 *pl;
	int           ix, iz, nmax;
	double       *X, **Psi;
	double        PsiAxis, PsiLim;
	double        Ptemp, Gtemp;
	double        dPress;
	double        dG2;

	pg = td->PsiGrid;
	pl = td->Plasma;
	X = pg->X;
	Psi = pg->Psi;
	PsiAxis = pg->PsiAxis;
	PsiLim = pg->PsiLim;

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
				dPress = Pp_DipoleStd_Loc(td,Psi[ix][iz]);
				dG2    = G2p_DipoleStd_Loc(td,Psi[ix][iz]);
				J[ix][iz] = Ptemp * dPress + Gtemp * dG2;
			} else
				J[ix][iz] = 0.0;
		}
	}
}

double       J_DipoleStd_Loc(TOKAMAK * td, int ix, int iz)
{
	PSIGRID      *pg;
	PLASMA		 *pl;
	double       *X, **Psi;
	double        PsiAxis, PsiLim;
	double        Ptemp, Gtemp;
	double        dPress;
	double        dG2;
	double		  J=0.0;

	pg = td->PsiGrid;
	pl = td->Plasma;
	X = pg->X;
	Psi = pg->Psi;
	PsiAxis = pg->PsiAxis;
	PsiLim = pg->PsiLim;


	if (pg->IsPlasma[ix][iz]) {
		Gtemp = -PI * DSQR(pl->B0R0) / X[ix];
		Ptemp = -TWOPI * X[ix];
		dPress = Pp_DipoleStd_Loc(td,Psi[ix][iz]);
		dG2    = G2p_DipoleStd_Loc(td,Psi[ix][iz]);
		J = (Ptemp * dPress + Gtemp * dG2);
	}
	return J;
}
