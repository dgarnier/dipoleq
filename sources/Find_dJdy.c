
/*
** TokaMac v2.0
**
** Find_dJdy.c
**
**
**
** File:		Find_dJdy.c
** Date:		March 22, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdlib.h>
#include "VAX_Alloc.h"
#include <math.h>
#include "nrutil.h"
#include "psigrid.h"
#include "plasma.h"
#include "tokamak.h"
#include "Find_dJdy.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

/*
**	new_dJdy
**
**	Allocates memory for a new dJdy array.
**	This is a three-dimensional array.
**	ynum is the number of unknowns, and the dimensions range
**	from [1..NumUnkns].
**
*/
double     ***new_dJdy(int ynum, int nmax)
{
	double     ***a;

	a = d3tensor(1, ynum, 0, nmax, 0, nmax);

	return a;
}

void          free_dJdy(double ***a, int ynum, int nmax)
{

	if (a) {
		free_d3tensor(a, 1, ynum, 0, nmax, 0, nmax);
		a = NULL;
	}
}

void          Find_dJdy(TOKAMAK * td, double ***dJdy)
{
	PSIGRID      *pg;
	PLASMA       *pl;
	double       *X, **Psi, PsiAxis, DelPsi, PsiX;
	double        c1, c2, c3, c4;	/* various constants */
	int         **ip, i, idx, ix, iz, nmax;

	pg = td->PsiGrid;
	PsiAxis = pg->PsiAxis;
	DelPsi = pg->PsiLim - PsiAxis;
	ip = pg->IsPlasma;
	Psi = pg->Psi;
	X = pg->X;
	nmax = pg->Nsize;
	pl = td->Plasma;

	for (i = 1; i <= td->NumUnkns; i++) {
		for (ix = 0; ix <= nmax; ix++)
			dJdy[i][ix][0] = dJdy[i][ix][nmax] = 0.0;
		for (iz = 0; iz <= nmax; iz++)
			dJdy[i][0][iz] = dJdy[i][nmax][iz] = 0.0;
	}

	switch (pl->ModelType) {

	  case Plasma_Std:
		  c1 = pl->B0R0;
		  c1 = c1 * c1 / 2.0;
		  for (ix = 1; ix < nmax; ix++) {
			  c2 = -TWOPI * X[ix];
			  c3 = -TWOPI * c1 / X[ix];
			  for (iz = 1; iz < nmax; iz++) {
				  if (ip[ix][iz]) {
					  PsiX = (Psi[ix][iz] - PsiAxis) / DelPsi;
					  /* --------  G2p[1] ------- */
					  dJdy[1][ix][iz] = c3 * pow(1.0 - PsiX, pl->StndG - 1.0);
					  /* --------   Pp[1] ------- */
					  dJdy[2][ix][iz] = c2 * pow(1.0 - PsiX, pl->StndP - 1.0);
				  } else {
					  /* --------  G2p[1] ------- */
					  dJdy[1][ix][iz] = 0.0;
					  /* --------   Pp[1] ------- */
					  dJdy[2][ix][iz] = 0.0;
				  }
			  }
		  }
		  break;

	  case Plasma_IsoNoFlow:
		  c1 = pl->B0R0;
		  c1 = c1 * c1 / 2.0;
		  for (ix = 1; ix < nmax; ix++) {
			  c2 = -TWOPI * X[ix];
			  c3 = -TWOPI * c1 / X[ix];
			  for (iz = 1; iz < nmax; iz++) {
				  if (ip[ix][iz]) {
					  PsiX = (Psi[ix][iz] - PsiAxis) / DelPsi;
					  /* --------  G2p[i] ------- */
					  idx = 0;
					  c4 = pow(PsiX, pl->G2pTerms);
					  for (i = 1; i <= pl->G2pTerms; i++)
						  dJdy[i + idx][ix][iz] = c3 * (pow(PsiX, i - 1) - c4);
					  /* --------   Pp[i] ------- */
					  idx = pl->G2pTerms;
					  c4 = pow(PsiX, pl->PpTerms);
					  for (i = 1; i <= pl->PpTerms; i++)
						  dJdy[i + idx][ix][iz] = c2 * (pow(PsiX, i - 1) - c4);
				  } else {
					  /* --------  G2p[i] ------- */
					  idx = 0;
					  for (i = 1; i <= pl->G2pTerms; i++)
						  dJdy[i + idx][ix][iz] = 0.0;
					  /* --------   Pp[i] ------- */
					  idx = pl->G2pTerms;
					  for (i = 1; i <= pl->PpTerms; i++)
						  dJdy[i + idx][ix][iz] = 0.0;
				  }
			  }
		  }
		  break;

	  default:
	  	nrinfo("Called Find_dJdy with unsupported plasma model\n");
	}
}
