/*
** TokaMac v2.0
**
** meas_ppsix.c
**
** Source file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_ppsix_Green(TOKAMAK *td, MEAS *m)
**
**				meas_ppsix_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_ppsix_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_ppsix.c
** Date:		March 24, 1993
**
** Revisions:
**
**		August 3, 1993		Only allowed pressure when IsPlasma == TRUE
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdio.h>
#include <math.h>
#include "fpoly.h"
#include "measurement.h"
#include "psigrid.h"
#include "plasma.h"
#include "tokamak.h"
#include "meas_ppsix.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

#define P_EDGE		0.0

/*
** Local variables to this file.
*/

/*
**	meas_ppsix_Fit
**
**
*/
void          meas_ppsix_Fit(TOKAMAK * td, MEAS * m)
{
	double        PsiX, DelPsi;
	double        press;
	int         **ip;
	PLASMA       *pl;
	PSIGRID      *pg;

	pl = td->Plasma;
	pg = td->PsiGrid;
	ip = pg->IsPlasma;

	switch (pl->ModelType) {
	  case Plasma_Std:
		  DelPsi = pg->PsiLim - pg->PsiAxis;
		  PsiX = m->X*m->X;
		  if (PsiX < 1.0)
			  press = -(DelPsi / pl->StndP) * pow(1.0 - PsiX, pl->StndP) * pl->Pp[1];
		  else
			  press = 0.0;
		  break;
	  case Plasma_IsoNoFlow:
		  DelPsi = pg->PsiLim - pg->PsiAxis;
		  PsiX = m->X*m->X;
		  if (PsiX < 1.0)
			  press = fpoly_int(pl->Pp, PsiX, pl->PpTerms, DelPsi, P_EDGE);
		  else
			  press = P_EDGE;
		  break;
	  case Plasma_IsoFlow:
		  break;
	  case Plasma_AnisoNoFlow:
		  break;
	  case Plasma_AnisoFlow:
		  break;
	}

	m->Fit = press / MU0;		/* pressure in Pascals */
}

/*
**	meas_ppsix_Now
**
**
*/
void          meas_ppsix_Now(TOKAMAK * td, MEAS * m)
{
	meas_ppsix_Fit(td, m);
	m->Now = m->Fit;
}

/*
**	meas_ppsix_L
**
**	Note:	The L array has dimensions [1..NumUnkns].
**
*/
void          meas_ppsix_L(TOKAMAK * td, MEAS * m, double *L)
{
	double        PsiX;
	int           idx, iu, pt, **ip;
	double        DelPsi, c4;
	PLASMA       *pl;
	PSIGRID      *pg;

	pl = td->Plasma;
	pg = td->PsiGrid;
	ip = pg->IsPlasma;
	pt = pl->PpTerms;

	for (iu = 1; iu <= td->NumUnkns; iu++)
		L[iu] = 0.0;

	switch (pl->ModelType) {
	  case Plasma_Std:
		  /* --------   L[2] = Pp[1]  ------- */
		  DelPsi = pg->PsiLim - pg->PsiAxis;
		  PsiX = m->X*m->X;
		  if (PsiX < 1.0)
			  L[2] = -(DelPsi / pl->StndP) * pow(1.0 - PsiX, pl->StndP);
		  break;
	  case Plasma_IsoNoFlow:
		  DelPsi = pg->PsiLim - pg->PsiAxis;
		  PsiX = m->X*m->X;
		  if (PsiX < 1.0) {
			  /* --------   Pp[i] ------- */
			  idx = pl->G2pTerms;
			  c4 = pow(PsiX, pt);
			  for (iu = 1; iu <= pt; iu++)
				  L[iu + idx] = -DelPsi * ((1.0 - pow(PsiX, iu)) / iu -
								   (1.0 - pow(PsiX, pt + 1)) / (pt + 1));
		  }
		  break;
	  case Plasma_IsoFlow:
		  break;
	  case Plasma_AnisoNoFlow:
		  break;
	  case Plasma_AnisoFlow:
		  break;
	}

	for (iu = 1; iu <= td->NumUnkns; iu++)
		L[iu] /= MU0;
}
