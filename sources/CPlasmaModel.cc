#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "CPlasmaModel.h"
#include "CDipoleStd.h"
#include "CDipoleIntStable.h"
#include "CDipoleStablePsiN.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

CPlasmaModel * CPlasmaModel::CreateModel(PLASMA *p)
{
	CPlasmaModel  *it = NULL;

	switch (p->ModelType) {
		case Plasma_DipoleStd :
			it = new CDipoleStd(p);
		break;
		case Plasma_DipoleIntStable :
			it = new CDipoleIntStable(p);
		break;
		case Plasma_DipoleStablePsiN :
			it = new CDipoleStablePsiN(p);
	}
	return it;
}

void CPlasmaModel::UpdateModel(TOKAMAK *td) {
	if (td->VacuumOnly != 0) isVacuum = 1;
}

void CPlasmaModel::FindJ(TOKAMAK *td, double **J)
{
	PSIGRID      *pg;
	PLASMA		 *pl;
	int           ix, iz, nmax;
	double       *X, **Psi;
	double        Ptemp, Gtemp;
	double        dPress;
	double        dG2;

	UpdateModel(td);

	pg = td->PsiGrid;
	pl = td->Plasma;
	X = pg->X;
	Psi = pg->Psi;

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
			if (pg->IsPlasma[ix][iz] && !isVacuum) {
				dPress = Pp(Psi[ix][iz]);
				dG2    = G2p(Psi[ix][iz]);
				J[ix][iz] = Ptemp * dPress + Gtemp * dG2;
			} else
				J[ix][iz] = 0.0;
		}
	}
}

double CPlasmaModel::FindJ_Loc(TOKAMAK * td, int ix, int iz)
{
	PSIGRID      *pg;
	PLASMA		 *pl;
	double       *X, **Psi;
	double        Ptemp, Gtemp;
	double        dPress;
	double        dG2;
	double		  J=0.0;

	pg = td->PsiGrid;
	pl = td->Plasma;
	X = pg->X;
	Psi = pg->Psi;

	if (isVacuum) return 0;

	if (pg->IsPlasma[ix][iz] && !isVacuum) {
		Gtemp = -PI * DSQR(pl->B0R0) / X[ix];
		Ptemp = -TWOPI * X[ix];
		dPress = Pp(Psi[ix][iz]);
		dG2    = G2p(Psi[ix][iz]);
		J = (Ptemp * dPress + Gtemp * dG2);
	}
	return J;
}

/*
** G E T P P A R M
**
**	Compute:
**			P		Pressure
**			G		Normalized toroidal flux
**			B2		Square mod-B
*/

void CPlasmaModel::GetPParam(TOKAMAK *td)
{
	PLASMA       *pl;
	PSIGRID      *pg;
	int           nmax, ix, iz;
	double      **Psi, **Press, **gPsiX, **gPsiZ, **gPsi2, **G, **Bt, **B2;
	double        dx, dz, PsiD, PsiW;
	int         **ip;

	pg = td->PsiGrid;
	pl = td->Plasma;
	nmax = pg->Nsize;
	dx = pg->dx;
	dz = pg->dz;
	PsiW    = pg->PsiLim;
	PsiD    = pg->PsiAxis;
	Psi = pg->Psi;
	ip = pg->IsPlasma;
	Press = pl->Piso;
	gPsiX = pl->GradPsiX;
	gPsiZ = pl->GradPsiZ;
	gPsi2 = pl->GradPsi2;
	G = pl->G;
	Bt = pl->Bt;
	B2 = pl->B2;

	/*  P R E S S U R E,  G,  E T C . . . .  */
	for (ix = 0; ix <= nmax; ix++)
		for (iz = 0; iz <= nmax; iz++) {
			if ((ip[ix][iz]) && (!isVacuum)) {
				Press[ix][iz] = P(Psi[ix][iz]);
				G[ix][iz] = sqrt(G2(Psi[ix][iz]));
			} else {
				Press[ix][iz] = 0.0;
				G[ix][iz] = 1.0;
			}
			Bt[ix][iz] = G[ix][iz] * pl->B0R0 / pg->X[ix];
			B2[ix][iz] = gPsi2[ix][iz] / TWOPI / pg->X[ix] / TWOPI / pg->X[ix];
			B2[ix][iz] = B2[ix][iz] + DSQR(Bt[ix][iz]);
		}
}
