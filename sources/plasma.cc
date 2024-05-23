/*
** TokaMac v2.0
**
** plasma.c
**
** This file contains routines which define the plasma
** data structure.
**
** File:		plasma.c
** Date:		February 10, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdlib.h>
#include <string.h>
#include "VAX_Alloc.h"
#include <math.h>
#include "nrutil.h"
#include "plasma.h"
#include "CPlasmaModel.h"
#include "fpoly.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

#define IF_free_dmatrix(x)		if (x) free_dmatrix(x, 0, nmax, 0, nmax)
#define IF_free_dvector(x)		if (x) free_dvector(x, 0, MaxPolyTerms)
#define IF_free_dvectorg(x,n)	if (x) free_dvector(x, 0, n)

/*
**
**	new_Plasma
**
*/
PLASMA       *new_Plasma()
{
	PLASMA       *p;

	p = (PLASMA *) malloc((unsigned) sizeof(PLASMA));
	memset(p, 0, sizeof(PLASMA));  // zero out the structure

	if (!p)
		nrerror("ERROR: Allocation error in new_Plasma.");

	/* Defaults */
	p->Nsize = 16;
	p->ModelType = Plasma_Std;
	p->Model = NULL;

	p->PpTerms = 3;
	p->G2pTerms = 1;
	p->SisoTerms = 0;
	p->SperTerms = 0;
	p->SparTerms = 0;
	p->RotTerms = 0;
	p->HTerms = 0;

	p->Pp = NULL;
	p->G2p = NULL;
	p->Rot = NULL;
	p->Siso = NULL;
	p->Sper = NULL;
	p->Spar = NULL;
	p->H = NULL;

	p->B2 = NULL;
	p->GradPsiX = NULL;
	p->GradPsiZ = NULL;
	p->GradPsi2 = NULL;
	p->Bt = NULL;
	p->G = NULL;
	p->Rho = NULL;
	p->Piso = NULL;
	p->Ppar = NULL;
	p->Pper = NULL;
	p->Alpha = NULL;

	p->StndP = 1.5;
	p->StndG = 1.5;
	p->Jedge = 0.0;
	p->ChiSqr = 1.0;

	p->NumBndMomts = 5;			/* moments from 0 to 4 */
	p->NumPsiPts = 10;			/* ten points to plot */
	p->PsiXmax = 0.98;			/* 98% flux surface */
	p->q_pr = NULL;				/* safety factor profile */
	p->Volp_pr = NULL;			/* flux derivative of flux surface volume */
	p->Vol_pr = NULL;			/* flux surface volume */
	p->S_pr = NULL;				/* global shear of q */
	p->B2_pr = NULL;			/* flux surface average of B2 */
	p->Well_pr = NULL;			/* magnetic well */
	p->J_pr = NULL;				/* flux-surface averaged toroidal current */
	p->Beta_pr = NULL;			/* flux-tube averaged beta */

	p->BBetaMax_pr = NULL;		/* beta at BetaMax on a flux surface */
	p->BMax_pr = NULL;			/* maximum field strength on a flux surface */


	p->Psi_pr = NULL;
	p->PsiX_pr = NULL;
	p->P_pr = NULL;
	p->G_pr = NULL;
	p->Pp_pr = NULL;
	p->G2p_pr = NULL;

	p->R0 = 0.2350;
	p->Z0 = 0.0;
	p->a0 = 0.5;
	p->Ip0 = 6.5e3;
	p->B0 = 4.0;
	p->B0R0 = p->B0 * p->R0;

	p->Ip = p->Ip0;
	p->beta0 = 0.00;
	p->beta = 0.000;
	p->betap = 1.0;
	p->li = 1.0;
	p->Ltotal = MU0 * p->R0;
	p->mu = 0.0;

	p->Volume = 2.0 * PI * PI * p->R0 * p->a0 * p->a0;
	p->Diamag = 0.0;
	p->q0 = 1.0;
	p->XMagAxis = p->R0;
	p->ZMagAxis = p->Z0;
	p->PsiAxis = 0.0;
	p->PsiLim = 0.1;
	p->RSurfaceAvg = p->R0;
	p->HalfWidth = p->a0;
	p->Elongation = 1.0;
	p->RCentroid = p->R0;
	p->qCircular = 2.0;

	p->TotKinEnergy = 1.0e6;
	p->TotMagEnergy = 1.0e6;
	p->RStar = p->R0;

	return p;
}

/*
**
**	free_Plasma
**
*/

void          free_Plasma(PLASMA * p)
{
	int           nmax = p->Nsize;

	if (p) {
		IF_free_dvector(p->Pp);
		IF_free_dvector(p->G2p);

		IF_free_dvector(p->Rot);
		IF_free_dvector(p->Siso);
		IF_free_dvector(p->Sper);
		IF_free_dvector(p->Spar);
		IF_free_dvector(p->H);

		IF_free_dvectorg(p->q_pr, p->NumPsiPts - 1);
		IF_free_dvectorg(p->Volp_pr, p->NumPsiPts - 1);
		IF_free_dvectorg(p->Vol_pr, p->NumPsiPts - 1);
		IF_free_dvectorg(p->S_pr, p->NumPsiPts - 1);
		IF_free_dvectorg(p->B2_pr, p->NumPsiPts - 1);
		IF_free_dvectorg(p->Well_pr, p->NumPsiPts - 1);
		IF_free_dvectorg(p->J_pr, p->NumPsiPts - 1);
		IF_free_dvectorg(p->Beta_pr, p->NumPsiPts - 1);

		IF_free_dvectorg(p->Psi_pr, p->NumPsiPts - 1);
		IF_free_dvectorg(p->PsiX_pr, p->NumPsiPts - 1);
		IF_free_dvectorg(p->P_pr, p->NumPsiPts - 1);
		IF_free_dvectorg(p->G_pr, p->NumPsiPts - 1);
		IF_free_dvectorg(p->Pp_pr, p->NumPsiPts - 1);
		IF_free_dvectorg(p->G2p_pr, p->NumPsiPts - 1);

		IF_free_dmatrix(p->B2);
		IF_free_dmatrix(p->GradPsiX);
		IF_free_dmatrix(p->GradPsiZ);
		IF_free_dmatrix(p->GradPsi2);
		IF_free_dmatrix(p->Bt);

		IF_free_dmatrix(p->G);
		IF_free_dmatrix(p->Rho);
		IF_free_dmatrix(p->Piso);
		IF_free_dmatrix(p->Ppar);
		IF_free_dmatrix(p->Pper);
		IF_free_dmatrix(p->Alpha);

		free(p);
	}
}

/*
**
**	init_Plasma
**
*/
void          init_Plasma(PLASMA * p)
{
	int           nmax = p->Nsize;

	p->B0R0 = p->B0 * p->R0;

	/***
	* The plasma model
	***/
	p->G2p = dvector0(0, MaxPolyTerms);
	p->H = dvector0(0, MaxPolyTerms);
	p->Rot = dvector0(0, MaxPolyTerms);
	p->Siso = dvector0(0, MaxPolyTerms);
	p->Sper = dvector0(0, MaxPolyTerms);
	p->Spar = dvector0(0, MaxPolyTerms);
	p->Pp = dvector0(0, MaxPolyTerms);

	/*
	** The computed plasma parameters
	**
	*/
	p->B2 = dmatrix0(0, nmax, 0, nmax);
	p->GradPsiX = dmatrix0(0, nmax, 0, nmax);
	p->GradPsiZ = dmatrix0(0, nmax, 0, nmax);
	p->GradPsi2 = dmatrix0(0, nmax, 0, nmax);

	p->Bt = dmatrix0(0, nmax, 0, nmax);
	p->G = dmatrix0(0, nmax, 0, nmax);

	switch (p->ModelType) {
	  case Plasma_Std:
	  case Plasma_IsoNoFlow:
		  p->Piso = dmatrix(0, nmax, 0, nmax);
		  break;
	  case Plasma_IsoFlow:
		  p->Piso = dmatrix(0, nmax, 0, nmax);
		  p->Rho = dmatrix(0, nmax, 0, nmax);
		  break;
	  case Plasma_AnisoNoFlow:
		  p->Rho = dmatrix(0, nmax, 0, nmax);
		  p->Ppar = dmatrix(0, nmax, 0, nmax);
		  p->Pper = dmatrix(0, nmax, 0, nmax);
		  p->Alpha = dmatrix0(0, nmax, 0, nmax);
		  break;
	  case Plasma_AnisoFlow:
		  p->Rho = dmatrix(0, nmax, 0, nmax);
		  p->Ppar = dmatrix(0, nmax, 0, nmax);
		  p->Pper = dmatrix(0, nmax, 0, nmax);
		  p->Alpha = dmatrix0(0, nmax, 0, nmax);
		  break;
	  default:
		  p->Piso = dmatrix(0, nmax, 0, nmax);
		  p->Model = CPlasmaModel::CreateModel(p);
	}
}


#define P_EDGE		0.0


double PlasmaP(PLASMA *pl, double Psi)
{
	double PsiX, DelPsi;

	switch (pl->ModelType) {
		case Plasma_Std:
			DelPsi =  pl->PsiLim - pl->PsiAxis;
		    PsiX = (Psi - pl->PsiAxis)/DelPsi;
			return (-DelPsi * pl->Pp[1] * pow(1.0 - PsiX, pl->StndP) / pl->StndP);
		break;

		case  Plasma_IsoNoFlow:
			DelPsi =  pl->PsiLim - pl->PsiAxis;
		    PsiX = (Psi - pl->PsiAxis)/DelPsi;
			return ( fpoly_int(pl->Pp, PsiX, pl->PpTerms, DelPsi, P_EDGE) );
		break;

		default :
		    if (pl->Model) return ( pl->Model->P(Psi) );
	}
	return (0);
}

double PlasmaPp(PLASMA *pl, double Psi)
{
	double PsiX, DelPsi;

	switch (pl->ModelType) {
		case Plasma_Std:
			DelPsi =  pl->PsiLim - pl->PsiAxis;
		    PsiX = (Psi - pl->PsiAxis)/DelPsi;
			return ( pl->Pp[1] * pow(1.0 - PsiX, pl->StndP - 1.0));
		break;

		case  Plasma_IsoNoFlow:
			DelPsi =  pl->PsiLim - pl->PsiAxis;
		    PsiX = (Psi - pl->PsiAxis)/DelPsi;
			return ( fpoly(pl->Pp, PsiX, pl->PpTerms) );
		break;

		default :
		    if (pl->Model) return ( pl->Model->Pp(Psi) );
	}
	return (0);
}

double PlasmaG(PLASMA *pl, double Psi)
{
	double PsiX, DelPsi;

	switch (pl->ModelType) {
		case Plasma_Std:
			DelPsi =  pl->PsiLim - pl->PsiAxis;
		    PsiX = (Psi - pl->PsiAxis)/DelPsi;
			return (1.0 - DelPsi * pl->G2p[1] * pow(1.0 - PsiX, pl->StndG) / pl->StndG);
		break;

		case  Plasma_IsoNoFlow:
			DelPsi =  pl->PsiLim - pl->PsiAxis;
		    PsiX = (Psi - pl->PsiAxis)/DelPsi;
			return ( fpoly_int(pl->G2p, PsiX, pl->G2pTerms, DelPsi, 1.0));
		break;

		default :
		    if (pl->Model) return ( pl->Model->G2(Psi) );
	}
	return (0);
}

double PlasmaG2p(PLASMA *pl, double Psi)
{
	double PsiX, DelPsi;

	switch (pl->ModelType) {
		case Plasma_Std:
			DelPsi =  pl->PsiLim - pl->PsiAxis;
		    PsiX = (Psi - pl->PsiAxis)/DelPsi;
			return (pl->G2p[1] * pow(1.0 - PsiX, pl->StndG - 1.0));
		break;

		case  Plasma_IsoNoFlow:
			DelPsi =  pl->PsiLim - pl->PsiAxis;
		    PsiX = (Psi - pl->PsiAxis)/DelPsi;
			return ( fpoly(pl->G2p, PsiX, pl->G2pTerms) );
		break;

		default :
		    if (pl->Model) return ( pl->Model->G2p(Psi) );
	}
	return (0);
}
