/*
** TokaMac v2.0
**
** meas_J0.c
**
** Source file for on-axis plasma current "measurements".
**
** Every measurement must define the following subroutines:
**
**				meas_J0_Green(TOKAMAK *td, MEAS *m)
**
**				meas_J0_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_J0_L(TOKAMAK *td, MEAS *m, double *L)
**
**
**	WARNING:  What happens when Jedge <> 0.0?
**
** File:		meas_J0.c
** Date:		February 21, 1996
**
** Revisions:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "fpoly.h"
#include "measurement.h"
#include "psigrid.h"
#include "plasma.h"
#include "tokamak.h"
#include "meas_J0.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

#define P_EDGE		0.0

/*
** Local variables to this file.
*/

/*
**	meas_J0_Fit
**
**
*/
void          meas_J0_Fit(TOKAMAK * td, MEAS * m)
{
	double        dPress,dG2,Ptemp,Gtemp;
	PLASMA       *pl;
	PSIGRID      *pg;

	pl = td->Plasma;
	pg = td->PsiGrid;

	Gtemp = -PI * DSQR(pl->B0R0) / pl->XMagAxis;
	Ptemp = -TWOPI * pl->XMagAxis;

	switch (pl->ModelType) {

	  case Plasma_Std:
			dPress = pl->Pp[1];
			dG2 = pl->G2p[1];
			break;

	  case Plasma_IsoNoFlow:
			dPress = fpoly(pl->Pp, 0.0, pl->PpTerms);
			dG2 = fpoly(pl->G2p, 0.0, pl->G2pTerms);
			break;

	  case Plasma_IsoFlow:
		  	break;

	  case Plasma_AnisoNoFlow:
		  	break;

	  case Plasma_AnisoFlow:
		  	break;

	}

	m->Fit = (Ptemp * dPress + Gtemp * dG2)/MU0;	/* on-axis current in A/m^2 */
}

/*
**	meas_J0_Now
**
**
*/
void          meas_J0_Now(TOKAMAK * td, MEAS * m)
{
	meas_J0_Fit(td, m);
	m->Now = m->Fit;
}

/*
**	meas_J0_L
**
**	Note:	The L array has dimensions [1..NumUnkns].
**
*/
void          meas_J0_L(TOKAMAK * td, MEAS *dummy, double *L)
{
	int           iu, pt, gt;
	double		 Gtemp,Ptemp;
	PLASMA       *pl;
	PSIGRID      *pg;

	pl = td->Plasma;
	pg = td->PsiGrid;
	pt = pl->PpTerms;
	gt = pl->G2pTerms;

	for (iu = 1; iu <= td->NumUnkns; iu++)
		L[iu] = 0.0;

	Gtemp = -PI * DSQR(pl->B0R0) / pl->XMagAxis;
	Ptemp = -TWOPI * pl->XMagAxis;

	switch (pl->ModelType) {

	  case Plasma_Std:
			/* --------  L[1] = G2p[1] ------- */
			L[1] = Gtemp;
			/* --------   L[2] = Pp[1]  ------- */
			L[2] = Ptemp;
			break;

	  case Plasma_IsoNoFlow:
			/* --------  G2p[i] ------- */
			L[1] = Gtemp;
			/* --------   Pp[i] ------- */
			L[1+gt] = Ptemp;
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
