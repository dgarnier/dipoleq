/*
** TokaMac v2.0
**
** meas_pnorm.c
**
** Source file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_pnorm_Green(TOKAMAK *td, MEAS *m)
**
**				meas_pnorm_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_pnorm_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_pnorm.c
** Date:		November 15, 1993
**
** Revisions:
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
#include "meas_pnorm.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

#define P_EDGE			0.0
#define PRESS_0_MIN		1.0
#define CHECK_PRESS_0 	if (press_0 < PRESS_0_MIN) press_0 = PRESS_0_MIN

/*
** Local variables to this file.
*/

/*
**	meas_pnorm_Fit
**
**
*/
void          meas_pnorm_Fit(TOKAMAK * td, MEAS * m)
{
	double        PsiX, DelPsi;
	double        press, press_0;
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
		  press_0 = -(DelPsi / pl->StndP) * pl->Pp[1];
		  break;
	  case Plasma_IsoNoFlow:
		  DelPsi = pg->PsiLim - pg->PsiAxis;
		  PsiX = m->X*m->X;
		  if (PsiX < 1.0)
			  press = fpoly_int(pl->Pp, PsiX, pl->PpTerms, DelPsi, P_EDGE);
		  else
			  press = P_EDGE;
		  press_0 = fpoly_int(pl->Pp, 0.0, pl->PpTerms, DelPsi, P_EDGE);
		  break;
	  case Plasma_IsoFlow:
		  break;
	  case Plasma_AnisoNoFlow:
		  break;
	  case Plasma_AnisoFlow:
		  break;
	}

	CHECK_PRESS_0;

	m->Fit = press / press_0;		/* normalized pressure */
}

/*
**	meas_pnorm_Now
**
**
*/
void          meas_pnorm_Now(TOKAMAK * td, MEAS * m)
{
	meas_pnorm_Fit(td, m);
	m->Now = m->Fit;
}

/*
**	meas_pnorm_L
**
**	Note:	The L array has dimensions [1..NumUnkns].
**
*/
void          meas_pnorm_L(TOKAMAK * td, MEAS * m, double *L)
{
	double        PsiX;
	int           idx, iu, pt, **ip;
	double        DelPsi, c4, pnorm, press_0;
	PLASMA       *pl;
	PSIGRID      *pg;

	pl = td->Plasma;
	pg = td->PsiGrid;
	ip = pg->IsPlasma;
	pt = pl->PpTerms;

	for (iu = 1; iu <= td->NumUnkns; iu++)
		L[iu] = 0.0;

	/* Use previously calculated m->Now, in order to find present pnorm */
	pnorm = m->Now;

	switch (pl->ModelType) {
	  case Plasma_Std:
		  /* --------   L[2] = Pp[1]  ------- */
		  DelPsi = pg->PsiLim - pg->PsiAxis;
		  press_0 = -(DelPsi / pl->StndP) * pl->Pp[1];
		  PsiX = m->X*m->X;
		  if (PsiX < 1.0)
			  L[2] = -(DelPsi / pl->StndP) * pow(1.0 - PsiX, pl->StndP);
		  L[2] += pnorm * (DelPsi / pl->StndP);
		  break;
	  case Plasma_IsoNoFlow:
		  DelPsi = pg->PsiLim - pg->PsiAxis;
		  press_0 = fpoly_int(pl->Pp, 0.0, pl->PpTerms, DelPsi, P_EDGE);
		  PsiX = m->X*m->X;
		  if (PsiX < 1.0) {
			  /* --------   Pp[i] ------- */
			  idx = pl->G2pTerms;
			  c4 = pow(PsiX, pt);
			  for (iu = 1; iu <= pt; iu++) {
				  L[iu + idx] = -DelPsi * ((1.0 - pow(PsiX, iu)) / iu -
								   (1.0 - pow(PsiX, pt + 1)) / (pt + 1));
				  L[iu + idx] += pnorm * DelPsi * (1.0 / iu - 1.0 / (pt + 1));
			  }
		  }
		  break;
	  case Plasma_IsoFlow:
		  break;
	  case Plasma_AnisoNoFlow:
		  break;
	  case Plasma_AnisoFlow:
		  break;
	}

	CHECK_PRESS_0;

	for (iu = 1; iu <= td->NumUnkns; iu++)
		L[iu] /= press_0;
}

#undef P_EDGE
#undef PRESS_0_MIN
#undef CHECK_PRESS_0
