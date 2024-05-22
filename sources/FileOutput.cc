/*
** 	TokaMac v2.0
**
** 	FileOutput.c
**
** 	This file defines subroutines to write
** 	various output quantities.
**
** 	File:		FileOutput.c
** 	Date:		April 7, 1993
**
**	Modifications:
**
**		April 10, 1993		Added BndMomOutput
**		May 3, 1993			Iterated to find (X0,Z0) in BndMomOutput
**		May 4, 1993			Added EQGRUMOutput
** 		May 10, 1993		Formated EQGRUMOutput as per Sabbagh
**		July 9, 1993		Fixed typo for mag axis in EQGRUMOutput
**		August 8, 1993		Fixed typo for Plasma_Std
**		Sept. 18, 1993		Added virial output
**		Feb. 16, 1996		Added DCON output file (with help from A. Glasser)
**		Feb. 1, 2002            Added GS2 output
**
** 	Routine list:
**
**		InValuesOutput		Reports the values input by the user.
**		PsiGridOutput		Reports values about the solution and IsPlasma.
**		PlasmaOutput		Reports many plasma parameters about solution.
**		FluxProfileOutput	Reports flux integrals like q(psi).
**		MeasOutput			Summaries the fit for each measurement.
**		BndMomentsOutput	Computes and reports the moments of outer flux surface.
**
** 	(c) L. Bai and M. Mauel -- Columbia University
*/


#include <stdlib.h>
#include "VAX_Alloc.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "DelChiSqr.h"
#include "fpoly.h"
#include "coil.h"
#include "shell.h"
#include "limiter.h"
#include "separatrix.h"
#include "measurement.h"
#include "plasma.h"
#include "psigrid.h"
#include "tokamak.h"
#include "CPlasmaModel.h"
#ifdef HDFOUTPUT
#include "HDFOutput.h"
#endif /* HDFOUTPUT */
#include "GetFluxMoments.h"
#include "FileOutput.h"
#include "J_DipoleStd.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

#define P_EDGE		0.0

/*
**	InValuesOutput
**
**
*/
void          InValuesOutput(TOKAMAK * td)
{
	FILE         *fi;
	PLASMA       *pl;
	PSIGRID      *pg;
	LIMITER      *lm;
	SEPARATRIX   *sp;
	COIL         *c;
	SUBCOIL      *sc;
	SHELL        *s;
	SUBSHELL     *ss;
	MEAS         *m;
	int           i, isc;
	char          fname[32] = "";

	pl = td->Plasma;
	pg = td->PsiGrid;

	strncat(fname, td->Oname, 18);	/* take 1st 18 characters */
	strcat(fname, "_InValues.out");

	fi = fopen(fname, "w");
	if (!fi)
		nrerror("ERROR:	Could not open file for writing in InValuesOut.");

	/*  H E A D E R */
	fprintf(fi, "TokaMac Values Output. From Input FileName: %s\n", td->Iname);
	fprintf(fi, "    RunName: %s. Info: %s\n", td->Name, td->Info);
	fprintf(fi, "    Run started at %s\n", td->Start);
	if (td->RestartStatus) {
		fprintf(fi, "    Run was retarted from file %s.\n", td->RSname);
		if (td->RestartUnkns)
			fprintf(fi, "    The previously saved unknowns were restored.\n");
		else
			fprintf(fi, "    The previously saved unknowns were NOT restored.\n");
	} else
		fprintf(fi, "    Run was initialized internally.\n\n");

	/* I N P U T   V A L U E S */
	fprintf(fi, "Input Values for run:\n");
	fprintf(fi, "    Fixed, Free Iterations = %d, %d\n", td->MaxIterFixed, td->MaxIterFree);
	fprintf(fi, "    Number of Coils, Shells, Limiters, Seps, Meas = %d, %d, %d, %d, %d\n",
			td->NumCoils, td->NumShells, td->NumLimiters, td->NumSeps, td->NumMeasures);
	fprintf(fi, "    LHGreens functions placed in %s\n", td->LHname);
	fprintf(fi, "    MeasGreens functions placed in %s\n", td->MGname);
	fprintf(fi, "    ShellGreens functions placed in %s\n", td->SGname);
	fprintf(fi, "    Shell inductance matrix placed in %s\n", td->SMname);
	fprintf(fi, "    Restart file placed in %s\n", td->RSname);
	fprintf(fi, "    Initial Plasma Current = %g (A), B0 = %.3f (T)\n",
			pl->Ip0, pl->B0);
	fprintf(fi, "    Initial Current Geometry R0 = %.3f Z0 = %.3f a0 = %.3f (m)\n",
			pl->R0, pl->Z0, pl->a0);
	switch (pl->ModelType) {
	  case Plasma_Std:
		  fprintf(fi, "    Plasma_Std with StdP = %.3f StdG = %.3f\n",
				  pl->StndP, pl->StndG);
		  break;
	  case Plasma_IsoNoFlow:
		  fprintf(fi, "    Plasma_IsoNoFlow with %d PpTerms %d G2Terms\n",
				  pl->PpTerms, pl->G2pTerms);
		  break;
	  case Plasma_IsoFlow:
		  fprintf(fi, "    Plasma_IsoFlow with H, G2, S, Rot Terms = %d, %d, %d, %d\n",
				  pl->HTerms, pl->G2pTerms, pl->SisoTerms, pl->RotTerms);
		  break;
	  case Plasma_AnisoNoFlow:
		  fprintf(fi, "    Plasma_AnisoNoFlow with H, G2, Sper, Spar Terms = %d, %d, %d, %d\n",
				  pl->HTerms, pl->G2pTerms, pl->SperTerms, pl->SparTerms);
		  break;
	  case Plasma_AnisoFlow:
		  fprintf(fi, "    Plasma_AnisoFlow with H, G2, Sper, Spar, Rot Terms = %d, %d, %d, %d, %d\n",
				  pl->HTerms, pl->G2pTerms, pl->SperTerms, pl->SparTerms, pl->RotTerms);
		  break;
	  default :
	  	  pl->Model->ModelDescription(fi);
//	  case Plasma_DipoleStd:
//		  fprintf(fi, "    Plasma_DipoleStd with %d PpTerms %d G2Terms\n",
//				  pl->PpTerms, pl->G2pTerms);
//		  break;

	}
	fprintf(fi, "    Approximate edge current density = %g (A/m2)\n",
			pl->Jedge/MU0);
	fprintf(fi, "    PsiGrid %d x %d  (%.3f < X < %.3f) (%.3f < Z < %.3f)\n",
			pg->Nsize, pg->Nsize, pg->Xmin, pg->Xmax, pg->Zmin, pg->Zmax);
	if (pg->Symmetric)
		fprintf(fi, "    PsiGrid is up-down symmetric.\n");
	else
		fprintf(fi, "    PsiGrid is not up-down symmetric.\n");
	fprintf(fi, "    Thresholds for boundary = %g, residual = %g\n",
			pg->BoundThreshold, pg->ResThreshold);
	fprintf(fi, "    The under relaxation parameters are %g, %g\n\n",
		pg->UnderRelax1, pg->UnderRelax2);

	/* L I M I T E R S */
	/*          1234567890123456789012345678901234567890123456789012345678901234567890 */
	fprintf(fi, "Limiter                     X1        Z1        X2        Z2   On/Off\n");
	for (i = 0; i < td->NumLimiters; i++) {
		lm = td->Limiters[i];
		fprintf(fi, "%20s %9.3f %9.3f %9.3f %9.3f   %d\n",
				lm->Name, lm->X1, lm->Z1, lm->X2, lm->Z2, lm->Enabled);
	}
	fprintf(fi, "\n");

	/* S E P A R A T R I E S */
	if (td->NumSeps > 0) {
		/*          1234567890123456789012345678901234567890123456789012345678901234567890 */
		fprintf(fi, "Separatrix   X1        Z1        X2        Z2        XC        ZC   On/Off\n");
		for (i = 0; i < td->NumSeps; i++) {
			sp = td->Seps[i];
			fprintf(fi, "%5s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f   %d\n",
				  sp->Name, sp->X1, sp->Z1, sp->X2, sp->Z2, sp->XC, sp->ZC, sp->Enabled);
		}
		fprintf(fi, "\n");
	}
	/* C O I L   S E T S */
	for (i = 0; i < td->NumCoils; i++) {
		c = td->Coils[i];
		fprintf(fi, "Coil Set %s On/Off = %d Coil Current = %g (A)\n",
				c->Name, c->Enabled, c->CoilCurrent / MU0);
		/*          1234567890123456789012345678901234567890123456789012345678901234567890 */
		fprintf(fi, "    SubCoil                  X         Z          Fract\n");
		for (isc = 0; isc < c->NumSubCoils; isc++) {
			sc = c->SubCoils[isc];
			fprintf(fi, "    %16s %9.3f %9.3f     %10.4f\n", sc->Name, sc->X, sc->Z, sc->CurrentFraction);
		}
		fprintf(fi, "\n");
	}

	/* S H E L L S */
	for (i = 0; i < td->NumShells; i++) {
		s = td->Shells[i];
		fprintf(fi, "Shell %s On/Off = %d\n", s->Name, s->Enabled);
		/*          1234567890123456789012345678901234567890123456789012345678901234567890 */
		fprintf(fi, "    SubShell                 X         Z         Radius\n");
		for (isc = 0; isc < s->NumSubShells; isc++) {
			ss = s->SubShells[isc];
			fprintf(fi, "    %16s %9.3f %9.3f     %10.4f\n", ss->Name, ss->X, ss->Z, ss->Radius);
		}
		fprintf(fi, "\n");
	}

	/* M E A S U R E M E N T S */
	/*          1234567890123456789012345678901234567890123456789012345678901234567890 */
	fprintf(fi, "Measurements         mType     Value    StdDev         X        Z\n");
	for (i = 0; i < td->NumMeasures; i++) {
		m = td->Measures[i];
		switch (m->mType) {
		  case meas_bp:
			  fprintf(fi, "    %16s   %3d %9.2g %9.2g %9.3f %9.3f %9.3f\n",
					  m->Name, m->mType, m->Value, m->StdDev, m->X, m->Z, m->parm.bp.Angle);
			  break;
		  case meas_bpangle:
			  fprintf(fi, "    %16s   %3d %9.2g %9.2g %9.3f %9.3f\n",
					  m->Name, m->mType, m->Value, m->StdDev, m->X, m->Z);
			  break;
		  case meas_flux:
			  fprintf(fi, "    %16s   %3d %9.2g %9.2g %9.3f %9.3f\n",
					  m->Name, m->mType, m->Value, m->StdDev, m->X, m->Z);
			  break;
		  case meas_saddle:
			  fprintf(fi, "    %16s   %3d %9.2g %9.2g %9.3f %9.3f %9.3f %9.3f\n",
					  m->Name, m->mType, m->Value, m->StdDev,
					  m->parm.saddle.X1, m->parm.saddle.Z1, m->parm.saddle.X2, m->parm.saddle.Z2);
			  break;
		  case meas_circle:
			  fprintf(fi, "    %16s   %3d %9.2g %9.2g %9.3f %9.3f %9.3f %d Ctype: %d\n",
					  m->Name, m->mType, m->Value, m->StdDev, m->X, m->Z, m->parm.circle.Radius,
					  m->parm.circle.Number, m->parm.circle.CircleType);
			  break;
		  case meas_press:
			  fprintf(fi, "    %16s   %3d %9.2g %9.2g %9.3f %9.3f\n",
					  m->Name, m->mType, m->Value, m->StdDev, m->X, m->Z);
			  break;
		  case meas_ppsix:   case meas_pnorm:
			  fprintf(fi, "    %16s   %3d %9.2g %9.2g %9.3f\n",
					  m->Name, m->mType, m->Value, m->StdDev, m->X);
			  break;
		  default:
			  fprintf(fi, "    %16s   %3d %9.2g %9.2g\n", m->Name, m->mType, m->Value, m->StdDev);
			  break;
		}
	}
	fprintf(fi, "\n");

	fclose(fi);
}

/*
**	ConductorsOutput
**
**
*/
void          ConductorsOutput(TOKAMAK * td)
{
	FILE         *fi;
	COIL         *c;
	SUBCOIL      *sc;
	SHELL        *s;
	SUBSHELL     *ss;
	int           i, isc;
	char          fname[32] = "";

	strncat(fname, td->Oname, 18);	/* take 1st 18 characters */
	strcat(fname, "_Conductors.out");

	fi = fopen(fname, "w");
	if (!fi)
		nrerror("ERROR:	Could not open file for writing in ConductorsOutput.");

	/*  H E A D E R */
	fprintf(fi, "TokaMac Conductors Output. From Input FileName: %s\n", td->Iname);
	fprintf(fi, "    RunName: %s. Info: %s\n", td->Name, td->Info);
	fprintf(fi, "    Run started at %s\n", td->Start);
	if (td->RestartStatus)
		fprintf(fi, "    Run was retarted from file %s.\n\n", td->RSname);
	else
		fprintf(fi, "    Run was initialized internally.\n\n");

	/* C O I L   S E T S */
	for (i = 0; i < td->NumCoils; i++) {
		c = td->Coils[i];
		fprintf(fi, "Coil Set %s On/Off = %d Coil Current = %g (A)\n",
				c->Name, c->Enabled, c->CoilCurrent / MU0);
		/*          1234567890123456789012345678901234567890123456789012345678901234567890 */
		fprintf(fi, "    SubCoil                  X         Z        Current\n");
		for (isc = 0; isc < c->NumSubCoils; isc++) {
			sc = c->SubCoils[isc];
			fprintf(fi, "    %16s %9.3f %9.3f     %10.4f\n", sc->Name, sc->X, sc->Z,
					sc->CurrentFraction * c->CoilCurrent / MU0);
		}
		fprintf(fi, "\n");
	}

	/* S H E L L S */
	for (i = 0; i < td->NumShells; i++) {
		s = td->Shells[i];
		fprintf(fi, "Shell %s On/Off = %d\n", s->Name, s->Enabled);
		/*          1234567890123456789012345678901234567890123456789012345678901234567890 */
		fprintf(fi, "    SubShell                 X         Z        Current\n");
		for (isc = 0; isc < s->NumSubShells; isc++) {
			ss = s->SubShells[isc];
			fprintf(fi, "    %16s %9.3f %9.3f     %10.4f\n", ss->Name, ss->X, ss->Z, ss->Current / MU0);
		}
		fprintf(fi, "\n");
	}

	fclose(fi);
}

/*
**	PsiGridOutput
*/
void          PsiGridOutput(TOKAMAK * td)
{
	PSIGRID      *pg;
	FILE         *fi;
	double      **J, **Psi, **Res;
	int         **IP;
	char          fname[32] = "";
	int           ix, iz, nmax;

	pg = td->PsiGrid;
	nmax = pg->Nsize;
	J = pg->Current;
	Psi = pg->Psi;
	Res = pg->Residual;
	IP = pg->IsPlasma;

	strncat(fname, td->Oname, 20);	/* take 1st 20 characters */
	strcat(fname, "_PsiGrid.out");

	fi = fopen(fname, "w");
	if (!fi)
		nrerror("ERROR:	Could not open file for writing in PsiGridOut.");

	/* H E A D E R */
	fprintf(fi, "PsiGrid Output. From Input FileName: %s\n", td->Iname);
	fprintf(fi, "    RunName: %s. Info: %s\n", td->Name, td->Info);
	fprintf(fi, "    Run started at %s\n", td->Start);
	if (td->RestartStatus)
		fprintf(fi, "    Run was retarted from file %s.\n", td->RSname);
	else
		fprintf(fi, "    Run was initialized internally.\n");
	fprintf(fi, "    Run ended at   %s\n\n", td->Stop);

	fprintf(fi, "    PsiLim = %g, PsiAxis = %g, DelPsi = %g\n",
			pg->PsiLim, pg->PsiAxis, pg->DelPsi);
	fprintf(fi, "    Symmetric = %d [1/0], UnderRelax1,2 = %g, %g\n",
			pg->Symmetric, pg->UnderRelax1, pg->UnderRelax2);
	fprintf(fi, "    MaxRes = %g, BoundError = %g\n\n",
			pg->MaxRes, pg->BoundError);

	/* I S P L A S M A */
	fprintf(fi, "IsPlasma\n");
	for (iz = nmax; iz >= 0; iz--) {
		for (ix = 0; ix <= nmax; ix++)
			fprintf(fi, "%1d", IP[ix][iz]);
		fprintf(fi, "\n");
	}
	fprintf(fi, "\n");

	fclose(fi);

#ifdef HDFOUTPUT
	HDFPsiGrid(pg, td->Oname);
#endif

}

/*
**	PlasmaOutput
*/
void          PlasmaOutput(TOKAMAK * td)
{
	PSIGRID      *pg;
	PLASMA       *pl;
	FILE         *fi;
	char          fname[32] = "";
	int            nmax, i;
	double        temp1;

//	if (td->VacuumOnly) return;
	pg = td->PsiGrid;
	pl = td->Plasma;
	nmax = pg->Nsize;

	strncat(fname, td->Oname, 20);	/* take 1st 20 characters */
	strcat(fname, "_Plasma.out");

	fi = fopen(fname, "w");
	if (!fi)
		nrerror("ERROR:	Could not open file for writing in PlasmaOut.");

	/* H E A D E R */
	fprintf(fi, "Plasma Output. From Input FileName: %s\n", td->Iname);
	fprintf(fi, "    RunName: %s. Info: %s\n", td->Name, td->Info);
	fprintf(fi, "    Run started at %s\n", td->Start);
	if (td->RestartStatus)
		fprintf(fi, "    Run was retarted from file %s.\n", td->RSname);
	else
		fprintf(fi, "    Run was initialized internally.\n");
	fprintf(fi, "    Run ended at   %s\n\n", td->Stop);

	/* G E O M E T R Y */
	fprintf(fi, "Geometry:\n");
	fprintf(fi, "    Magnetic axis (x, z) = (%g, %g)\n", pl->XMagAxis, pl->ZMagAxis);
	fprintf(fi, "    Current centroid R = %g (m), RCenter = %g (m)\n",
			pl->RCentroid, pl->RCenter);
	fprintf(fi, "    RStar = %g, RSurfaceAvg = %g, MinorHalfWidth = %g\n",
			pl->RStar, pl->RSurfaceAvg, pl->HalfWidth);
	fprintf(fi, "    Volume = %g (m3), Perimeter = %g (m).\n\n", pl->Volume, pl->Perimeter);

	/* E N E R G Y   &   B E T A */
	temp1 = pl->B0R0 / pl->RCenter;
	fprintf(fi, "Energy and beta:\n");
	fprintf(fi, "    Plasma current %g (A), Diamagnetic flux %g (weber)\n",
			pl->Ip, pl->Diamag);
	fprintf(fi, "    beta poloidal = %g, li = %g, mu = %g\n",
			pl->betap, pl->li, pl->mu);
	fprintf(fi, "    Volume-averaged beta = %g. BetaN = %g\n",
			pl->beta, pl->beta * 100.0 / (pl->Ip / 1.0e6 / pl->HalfWidth / temp1));
	fprintf(fi, "    Safety factor: q0 is %g; q at the %5.3f flux surface is %g\n",
			pl->q0, pl->PsiXmax, pl->q_pr[pl->NumPsiPts - 1]);
	fprintf(fi, "    Circular q is %g, and qStar is %g.\n", pl->qCircular, pl->qStar);
	fprintf(fi, "    Total stored kinetic energy = %g (J)\n", pl->TotKinEnergy);
	fprintf(fi, "    Total stored magnetic energy = %g (J)\n\n", pl->TotMagEnergy);

	/* V I R I A L   O U T P U T  */
	fprintf(fi, "Virial Outputs:\n");
	fprintf(fi, "    R_vr = %g (m), Alpha_vr = %g\n", pl->R_vr, pl->Alpha_vr);
	fprintf(fi, "    S1 = %g, S2 = %g, S3 = %g\n", pl->S1_vr, pl->S2_vr, pl->S3_vr);
	fprintf(fi, "    Lambda = %g, betap(perp) = %g\n",
		 0.25 * (pl->S1_vr + pl->S2_vr * (1.0 + pl->R_vr / pl->RCenter)),
			0.5 * (pl->S1_vr + pl->S2_vr * (1.0 - pl->R_vr / pl->RCenter)) + pl->mu);
	if (pl->Alpha_vr == 1.0)
		fprintf(fi, "    li(mu) = %g, li(alpha) = inf\n\n",
				-(0.5 * pl->S1_vr + 0.5 * pl->S2_vr * (1.0 - 3.0 * pl->R_vr / pl->RCenter) + 2.0 * pl->mu));
	else
		fprintf(fi, "    li(mu) = %g, li(alpha) = %g\n\n",
				-(0.5 * pl->S1_vr + 0.5 * pl->S2_vr * (1.0 - 3.0 * pl->R_vr / pl->RCenter) + 2.0 * pl->mu),
				(0.5 * pl->S1_vr + 0.5 * pl->S2_vr * (1.0 - pl->R_vr / pl->RCenter) - pl->S3_vr) / (pl->Alpha_vr - 1.0));

	/* R E C O N S T R U C T E D    M O D E L */
	fprintf(fi, "Plasma model:\n");
	switch (pl->ModelType) {
	  case Plasma_Std:
		  fprintf(fi, "    Plasma_Std with StdP = %.3f StdG = %.3f\n",
				  pl->StndP, pl->StndG);
		  break;
	  case Plasma_IsoNoFlow:
		  fprintf(fi, "    Plasma_IsoNoFlow with %d PpTerms %d G2Terms\n",
				  pl->PpTerms, pl->G2pTerms);
		  break;
	  case Plasma_IsoFlow:
		  fprintf(fi, "    Plasma_IsoFlow with H, G2, S, Rot Terms = %d, %d, %d, %d\n",
				  pl->HTerms, pl->G2pTerms, pl->SisoTerms, pl->RotTerms);
		  break;
	  case Plasma_AnisoNoFlow:
		  fprintf(fi, "    Plasma_AnisoNoFlow with H, G2, Sper, Spar Terms = %d, %d, %d, %d\n",
				  pl->HTerms, pl->G2pTerms, pl->SperTerms, pl->SparTerms);
		  break;
	  case Plasma_AnisoFlow:
		  fprintf(fi, "    Plasma_AnisoFlow with H, G2, Sper, Spar, Rot Terms = %d, %d, %d, %d, %d\n",
				  pl->HTerms, pl->G2pTerms, pl->SperTerms, pl->SparTerms, pl->RotTerms);
		  break;
	  default :
	  	  pl->Model->ModelDescription(fi);
//	  case Plasma_DipoleStd:
//		  fprintf(fi, "    Plasma_DipoleStd with %d PpTerms %d G2Terms\n",
//				  pl->PpTerms, pl->G2pTerms);
//		  break;

	}
	fprintf(fi, "\n");

	/*          1234567890123456789012345678901234567890123456789012345678901234567890 */
	fprintf(fi, "FluxF           0         1         2         3         4         5         6         7\n");
	fprintf(fi, "%7s", "G2p");
	for (i = 0; i <= MaxPolyTerms; i++)
		fprintf(fi, "%10.3g", pl->G2p[i]);
	fprintf(fi, "\n");
	if ((pl->ModelType == Plasma_IsoFlow) ||
		(pl->ModelType == Plasma_AnisoNoFlow) ||
		(pl->ModelType == Plasma_AnisoFlow)) {
		fprintf(fi, "%7s", "H");
		for (i = 0; i <= MaxPolyTerms; i++)
			fprintf(fi, "%10.3g", pl->H[i]);
		fprintf(fi, "\n");
	}
	if ((pl->ModelType == Plasma_IsoFlow) ||
		(pl->ModelType == Plasma_AnisoFlow)) {
		fprintf(fi, "%7s", "Rot");
		for (i = 0; i <= MaxPolyTerms; i++)
			fprintf(fi, "%10.3g", pl->Rot[i]);
		fprintf(fi, "\n");
	}
	if ((pl->ModelType == Plasma_AnisoNoFlow) ||
		(pl->ModelType == Plasma_AnisoFlow)) {
		fprintf(fi, "%7s", "Spar");
		for (i = 0; i <= MaxPolyTerms; i++)
			fprintf(fi, "%10.3g", pl->Spar[i]);
		fprintf(fi, "\n");
		fprintf(fi, "%7s", "Sper");
		for (i = 0; i <= MaxPolyTerms; i++)
			fprintf(fi, "%10.3g", pl->Sper[i]);
		fprintf(fi, "\n");
	}
	if (pl->ModelType == Plasma_IsoFlow) {
		fprintf(fi, "%7s", "Siso");
		for (i = 0; i <= MaxPolyTerms; i++)
			fprintf(fi, "%10.3g", pl->Siso[i]);
		fprintf(fi, "\n");
	}
	if ((pl->ModelType == Plasma_IsoNoFlow) ||
//	    (pl->ModelType == Plasma_DipoleStd) ||
		(pl->ModelType == Plasma_Std)) {
		fprintf(fi, "%7s", "Pp");
		for (i = 0; i <= MaxPolyTerms; i++)
			fprintf(fi, "%10.3g", pl->Pp[i]);
		fprintf(fi, "\n");
	}

	if (pl->Model)
		pl->Model->ModelOutput(fi);

	fclose(fi);
#ifdef HDFOUTPUT
	HDFPlasma(pl, pg, td->Oname);
#endif /* HDFOUTPUT */
}


/*
**	FluxProfileOutput
*/
void          FluxProfileOutput(TOKAMAK * td)
{
	FILE         *fi;
	PLASMA       *pl;
	int           i, npts;
	double        PsiXmax, PsiX, Psi, DelPsi, P, G, Pp, G2p;
	char          fname[32] = "";
	double		  *PV, *GV, *PpV, *G2V, *PsiV, *PsiXV;

	strncat(fname, td->Oname, 20);	/* take 1st 20 characters */
	strcat(fname, "_FluxPr.out");

	fi = fopen(fname, "w");
	if (!fi)
		nrerror("ERROR:	Could not open file for writing in FluxProfileOutput.");


	/*  H E A D E R */
	fprintf(fi, "TokaMac Flux Profile Output. From Input FileName: %s\n", td->Iname);
	fprintf(fi, "    RunName: %s. Info: %s\n", td->Name, td->Info);
	fprintf(fi, "    Run started at %s\n", td->Start);
	if (td->RestartStatus)
		fprintf(fi, "    Run was retarted from file %s.\n", td->RSname);
	else
		fprintf(fi, "    Run was initialized internally.\n");
	fprintf(fi, "    Run ended at   %s\n\n", td->Stop);


	/* F L U X   P R O F I L E S */
	pl = td->Plasma;
	npts = pl->NumPsiPts;
	PsiXmax = pl->PsiXmax;
	DelPsi = pl->PsiLim - pl->PsiAxis;

	PsiV  = pl->Psi_pr ;
	PsiXV = pl->PsiX_pr;
	PV    = pl->P_pr   ;
	GV    = pl->G_pr   ;
	PpV   = pl->Pp_pr  ;
	G2V   = pl->G2p_pr ;

#ifndef DIPOLE
	/*          12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890 */
	fprintf(fi, "i         PsiX         Psi           q           P           G          Pp         G2p        dVol         Vol          S         Well         <J>        <B2>\n");
#else
	fprintf(fi, "i         PsiX         Psi           P           Pp        dVol         Vol        <J>        <B2>        <Beta>     BetaMax     RBetaMax    ZBetaMax    BBetaMax	BMax         RBMax       ZBMax\n");
#endif
	for (i = 0; i < npts; i++) {
		PsiX = PsiXV[i];
		Psi = PsiV[i];
		P = PV[i];
		G = GV[i];
		Pp = PpV[i];
		G2p = G2V[i];
#ifndef DIPOLE
		fprintf(fi, "%02d %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g\n",
				i, PsiX, Psi, pl->q_pr[i], P, G, Pp, G2p, pl->Volp_pr[i], pl->Vol_pr[i],
				pl->S_pr[i], pl->Well_pr[i], pl->J_pr[i], pl->B2_pr[i]);
#else
		fprintf(fi, "%02d %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g %11.4g\n",
				i, PsiX, Psi, P, Pp, pl->Volp_pr[i], pl->Vol_pr[i],
			 pl->J_pr[i], pl->B2_pr[i], pl->Beta_pr[i],
			 pl->BetaMax_pr[i], pl->XBetaMax_pr[i], pl->ZBetaMax_pr[i], pl->BBetaMax_pr[i],
                         pl->BMax_pr[i],pl->XBMax_pr[i],pl->ZBMax_pr[i]);

#endif
	}
	fprintf(fi, "\n");

	fclose(fi);
#ifdef HDFOUTPUT
	HDFFluxFuncs(td->Oname,npts,PsiXV,PsiV,PV,GV,PpV,G2V,pl->q_pr,pl->Volp_pr, pl->Vol_pr,
			pl->S_pr, pl->Well_pr, pl->J_pr, pl->B2_pr, pl->Beta_pr,
			pl->BetaMax_pr, pl->XBetaMax_pr, pl->ZBetaMax_pr, pl->BBetaMax_pr,
                         pl->BMax_pr,pl->XBMax_pr,pl->ZBMax_pr);
#endif /* HDFOUTPUT */

}

/*
**	MeasOutput
*/
void          MeasOutput(TOKAMAK * td)
{
	FILE         *fi;
	MEAS         *m;
	int           i, deg_freedom;
	char          fname[32] = "";
	double        chi, chisum = 0.0;

	strncat(fname, td->Oname, 20);	/* take 1st 20 characters */
	strcat(fname, "_Meas.out");

	fi = fopen(fname, "w");
	if (!fi)
		nrerror("ERROR:	Could not open file for writing in MeasOutput.");

	deg_freedom = td->NumMeasures - td->NumUnkns;

	/*  H E A D E R */
	fprintf(fi, "TokaMac Meas Output. From Input FileName: %s\n", td->Iname);
	fprintf(fi, "    RunName: %s. Info: %s\n", td->Name, td->Info);
	fprintf(fi, "    Run started at %s\n", td->Start);
	if (td->RestartStatus)
		fprintf(fi, "    Run was retarted from file %s.\n", td->RSname);
	else
		fprintf(fi, "    Run was initialized internally.\n");
	fprintf(fi, "    The final value of chisq = %g; BoundError = %g\n",
			td->Plasma->ChiSqr, td->PsiGrid->BoundError);
        if (deg_freedom >= 1)
	fprintf(fi, "    The chi-square probability, Q = %g.\n",
		gammq(0.5*deg_freedom,0.5*td->Plasma->ChiSqr));
	fprintf(fi, "    Run ended at   %s\n\n", td->Stop);

	/* M E A S U R E M E N T S */
	/*          1234567890123456789012345678901234567890123456789012345678901234567890 */
	fprintf(fi, "Measurements         mType     Value    StdDev       Fit     |Chi|\n");
	for (i = 0; i < td->NumMeasures; i++) {
		m = td->Measures[i];
		chi = fabs((m->Value - m->Fit) / m->StdDev);
		chisum += chi * chi;
		fprintf(fi, "    %16s   %3d %9.2g %9.2g %9.2g %9.3f\n",
				m->Name, m->mType, m->Value, m->StdDev, m->Fit, chi);
	}
	fprintf(fi, "                                                         ---------\n");
	fprintf(fi, "Sum                                                      %9.3f\n",
			chisum);
	fprintf(fi, "\n");

	fclose(fi);
}

/*
**	BndMomentsOutput
*/
void          BndMomentsOutput(TOKAMAK * td)
{
	FILE         *fi;
	PSIGRID      *pg;
	int           m, nmomts;
	double       *Xm, *Zm;
	double        PsiXmax;
	char          fname[32] = "";

	strncat(fname, td->Oname, 20);	/* take 1st 20 characters */
	strcat(fname, "_BndMomts.out");

	fi = fopen(fname, "w");
	if (!fi)
		nrerror("ERROR:	Could not open file for writing in BndMomentsOutput.");

	/*  H E A D E R */
	fprintf(fi, "TokaMac Boundary Moments Output. From Input FileName: %s\n", td->Iname);
	fprintf(fi, "    RunName: %s. Info: %s\n", td->Name, td->Info);
	fprintf(fi, "    Run started at %s\n", td->Start);
	if (td->RestartStatus)
		fprintf(fi, "    Run was retarted from file %s.\n", td->RSname);
	else
		fprintf(fi, "    Run was initialized internally.\n");
	fprintf(fi, "    Run ended at   %s\n\n", td->Stop);

	/* B O U N D A R Y   M O M E N T S */
	pg = td->PsiGrid;
	nmomts = td->Plasma->NumBndMomts;
	PsiXmax = td->Plasma->PsiXmax;
	fprintf(fi, "First %d moments of the %g flux surface:\n\n", nmomts, PsiXmax);

	/*          12345678901234567890123456789012345678901234567890 */
	fprintf(fi, "m           Xm          Zm  \n");

	/* Find the moments */
	Xm = dvector(0, nmomts);
	Zm = dvector(0, nmomts);

	GetFluxMoments(pg, PsiXmax, Xm, Zm, nmomts);

	/* Print the moments */
	for (m = 0; m <= nmomts; m++)
		fprintf(fi, "%02d %11.4g %11.4g\n", m, Xm[m], Zm[m]);

	fprintf(fi, "\n");
        fclose(fi);


#ifdef HDFOUTPUT
        {
            int len;
            double *X,*Z;

#ifdef DIPOLE
            GetFluxContour(pg,0.0,&X,&Z,&len);
			if (X != NULL) {
				HDFBoundary(td->Oname,FCFS_NAME,0,X,Z,len);
				free_dvector(X,0,len);
				free_dvector(Z,0,len);
			} else {
				fprintf(stderr,"Warning: Could not find the flux contour for FCFS.\n");
			}
#endif /* DIPOLE */
            GetFluxContour(pg,PsiXmax,&X,&Z,&len);
            if (X != NULL) {
                HDFBoundary(td->Oname,LCFS_NAME,PsiXmax,X,Z,len);
                free_dvector(X,0,len);
                free_dvector(Z,0,len);
            } else {
				fprintf(stderr,"Warning: Could not find the flux contour for LCFS.\n");
			}
			HDFLimiters(td->Oname, td->Limiters, td->NumLimiters);
        }
#endif /* HDFOUTPUT */

	free_dvector(Zm, 0, nmomts);
	free_dvector(Xm, 0, nmomts);

}

/*
**	EQGRUMOutput
*/
void          EQGRUMOutput(TOKAMAK * td)
{
	FILE         *fi;
	PSIGRID      *pg;
	PLASMA       *pl;
	double        PsiXmax, PsiX, Psi, DelPsi, P, Pp;
	double        PsiAxis;
	int           m, nmomts, i, npts;
	double       *Xm, *Zm;
	char          fname[32] = "";

	strncat(fname, td->Oname, 20);	/* take 1st 20 characters */
	strcat(fname, "_EQGRUM.stb");

	fi = fopen(fname, "w");
	if (!fi)
		nrerror("ERROR:	Could not open file for writing in EQGRUMOutput.");

	/*  H E A D E R */
	fprintf(fi, " equilibrium    1 at time  0.0000001E+00 nstep    0\n");
	fprintf(fi, " +++ TOKAMAC    00000 at 0.0000E+00 seconds.  +++\n");
	fprintf(fi, " +++ TRANSP file:     Interpolation Order:  1 +++\n");
	fprintf(fi, " TokaMac RunName: %s from Input FileName: %s\n", td->Name, td->Iname);
	fprintf(fi, " +++ Smoothing (%%): Pressure     q      psi   +++\n");
	fprintf(fi, " +++                   0.00    0.00    0.00   +++\n");

	pg = td->PsiGrid;
	pl = td->Plasma;
	nmomts = td->Plasma->NumBndMomts;
	npts = pl->NumPsiPts;
	PsiXmax = pl->PsiXmax;
	DelPsi = pl->PsiLim - pl->PsiAxis;
	PsiAxis = pl->PsiAxis;

	/* V A L U E   N U M B E R S */
	fprintf(fi, " mjbal= %4d mombnd= %4d mom=    0\n", npts + 2, nmomts);

	/* M O M E N T S */
	Xm = dvector(0, nmomts);
	Zm = dvector(0, nmomts);

	GetFluxMoments(pg, PsiXmax, Xm, Zm, nmomts);

	fprintf(fi, " r0b= %14.7e", Xm[0]);
	for (m = 1; m <= nmomts; m++) {
		if (m % 5 == 1)
			fprintf(fi, "\n rmb=");
		fprintf(fi, " %14.7e", Xm[m]);
	}

	for (m = 1; m <= nmomts; m++) {
		if (m % 5 == 1)
			fprintf(fi, "\n ymb=");
		fprintf(fi, " %14.7e", Zm[m]);
	}
	fprintf(fi, "\n");

	/* Toroidal Field */
	fprintf(fi, " btor= %14.7e rtor= %14.7e eqcamp= %14.7e\n",
			pl->B0, pl->R0, pl->Ip);

	/* flux values */
	fprintf(fi, "  j      psi            q              d q / d psi    p              d p / d psi\n");

	/* inner guard point */
	PsiX = 1 * PsiXmax / (npts - 1) + 0.0001;	/* start at the 0.01% flux surface */
	Psi = pl->PsiAxis + PsiX * DelPsi;
	Pp = P = 0.0;
	if (pl->ModelType == Plasma_Std) {
		P = -DelPsi * pl->Pp[1] * pow(1.0 - PsiX, pl->StndP) / pl->StndP;
		Pp = pl->Pp[1] * pow(1.0 - PsiX, pl->StndP - 1.0);
	} else if (pl->ModelType == Plasma_IsoNoFlow) {
		P = fpoly_int(pl->Pp, PsiX, pl->PpTerms, DelPsi, P_EDGE);
		Pp = fpoly(pl->Pp, PsiX, pl->PpTerms);
	}
	P = P / MU0;
	Pp = Pp / MU0;
	fprintf(fi, " %4d %14.7e %14.7e %14.7e %14.7e %14.7e\n", 1,
			(Psi - PsiAxis) / TWOPI, pl->q_pr[1], PI * (pl->q_pr[2] - pl->q_pr[0]) * (npts - 1) / DelPsi,
			P, TWOPI * Pp);

	/* magnetic axis */
	PsiX = 0.0001;				/* start at the 0.01% flux surface */
	Psi = pl->PsiAxis + PsiX * DelPsi;
	Pp = P = 0.0;
	if (pl->ModelType == Plasma_Std) {
		P = -DelPsi * pl->Pp[1] * pow(1.0 - PsiX, pl->StndP) / pl->StndP;
		Pp = pl->Pp[1] * pow(1.0 - PsiX, pl->StndP - 1.0);
	} else if (pl->ModelType == Plasma_IsoNoFlow) {
		P = fpoly_int(pl->Pp, PsiX, pl->PpTerms, DelPsi, P_EDGE);
		Pp = fpoly(pl->Pp, PsiX, pl->PpTerms);
	}
	P = P / MU0;
	Pp = Pp / MU0;
	fprintf(fi, " %4d %14.7e %14.7e %14.7e %14.7e %14.7e\n", 2,
			0.0, pl->q_pr[0], TWOPI * (pl->q_pr[1] - pl->q_pr[0]) * (npts - 1) / DelPsi, P, TWOPI * Pp);

	/* the flux profiles */
	for (i = 1; i < npts - 1; i++) {
		PsiX = i * PsiXmax / (npts - 1) + 0.0001;	/* start at the 0.01% flux surface */
		Psi = pl->PsiAxis + PsiX * DelPsi;
		Pp = P = 0.0;
		if (pl->ModelType == Plasma_Std) {
			P = -DelPsi * pl->Pp[1] * pow(1.0 - PsiX, pl->StndP) / pl->StndP;
			Pp = pl->Pp[1] * pow(1.0 - PsiX, pl->StndP - 1.0);
		} else if (pl->ModelType == Plasma_IsoNoFlow) {
			P = fpoly_int(pl->Pp, PsiX, pl->PpTerms, DelPsi, P_EDGE);
			Pp = fpoly(pl->Pp, PsiX, pl->PpTerms);
		}
		P = P / MU0;
		Pp = Pp / MU0;
		fprintf(fi, " %4d %14.7e %14.7e %14.7e %14.7e %14.7e\n", 2 + i,
				(Psi - PsiAxis) / TWOPI, pl->q_pr[i], PI * (pl->q_pr[i + 1] - pl->q_pr[i - 1]) * (npts - 1) / DelPsi,
				P, TWOPI * Pp);
	}

	/* the last flux surface */
	PsiX = PsiXmax + 0.0001;	/* start at the 0.01% flux surface */
	Psi = pl->PsiAxis + PsiX * DelPsi;
	Pp = P = 0.0;
	if (pl->ModelType == Plasma_Std) {
		P = -DelPsi * pl->Pp[1] * pow(1.0 - PsiX, pl->StndP) / pl->StndP;
		Pp = pl->Pp[1] * pow(1.0 - PsiX, pl->StndP - 1.0);
	} else if (pl->ModelType == Plasma_IsoNoFlow) {
		P = fpoly_int(pl->Pp, PsiX, pl->PpTerms, DelPsi, P_EDGE);
		Pp = fpoly(pl->Pp, PsiX, pl->PpTerms);
	}
	P = P / MU0;
	Pp = Pp / MU0;
	fprintf(fi, " %4d %14.7e %14.7e %14.7e %14.7e %14.7e\n", npts + 1,
			(Psi - PsiAxis) / TWOPI, pl->q_pr[npts - 1], TWOPI * (pl->q_pr[npts - 1] - pl->q_pr[npts - 2]) * (npts - 1) / DelPsi,
			P, TWOPI * Pp);

	/* the outer guard point */
	fprintf(fi, " %4d %14.7e %14.7e %14.7e %14.7e %14.7e\n", npts + 2, 0.0, 0.0, 0.0, 0.0, 0.0);

	free_dvector(Zm, 0, nmomts);
	free_dvector(Xm, 0, nmomts);
	fclose(fi);
}

/*
**	DCONOutput
**      this code is obsolete and should be changed to reflect newest plasma models..
*/
void          DCONOutput(TOKAMAK * td)
{
	FILE         *fi;
	PSIGRID      *pg;
	PLASMA       *pl;
	double        PsiXmax, PsiX, Psi, DelPsi, P, G;
	double        PsiAxis;
	int            i, npts;
	int			   ix, iz, Nsize;
	char          fname[32] = "";

	strncat(fname, td->Oname, 20);	/* take 1st 20 characters */
	strcat(fname, "_DCON.dat");

	fi = fopen(fname, "w");
	if (!fi)
		nrerror("ERROR:	Could not open file for writing in DCONOutput.");

	/*  H E A D E R */
	fprintf(fi, " ** TokaMac Equilibrium Output prepared for DCON\n");
	fprintf(fi, " ** TokaMac RunName: %s from Input FileName: %s\n", td->Name, td->Iname);
	fprintf(fi, " ** ");
	fprintf(fi, " ** Output below consists of following values:\n");
	fprintf(fi, " **     mr, mz, mpsi   ! grid size for (r,z) and Psi\n");
	fprintf(fi, " **     rmin, rmax, zmin, zmax\n");
	fprintf(fi, " **     Volumn Averaged Beta\n");
	fprintf(fi, " **     Global Beta Poloidal\n");
	fprintf(fi, " **     -- blank line-----------------------------------\n");
	fprintf(fi, " **     PsiGrid(ir,iz) ! ir varies first, iz varies next\n");
	fprintf(fi, " **     -- blank line-----------------------------------\n");
	fprintf(fi, " **     Psi(ipsi)      ! Poloidal flux (MKS)\n");
	fprintf(fi, " **     -- blank line-----------------------------------\n");
	fprintf(fi, " **     F(ipsi)        ! Toroidal flux (MKS)\n");
	fprintf(fi, " **     -- blank line-----------------------------------\n");
	fprintf(fi, " **     P(ipsi)        ! Pressure (MKS)\n");
	fprintf(fi, " **     -- blank line-----------------------------------\n");
	fprintf(fi, " **     q(ipsi)        ! safety factor\n");
	fprintf(fi, " **     -- blank line-----------------------------------\n");
	fprintf(fi, " ** File format dated February 17, 1996.\n\n");

	pg = td->PsiGrid;
	pl = td->Plasma;
	npts = pl->NumPsiPts;
	PsiXmax = pl->PsiXmax;
	DelPsi = pl->PsiLim - pl->PsiAxis;
	PsiAxis = pl->PsiAxis;

	Nsize = pg->Nsize;

	/* Grid Size */
	fprintf(fi, " %6d %6d %6d\n",Nsize+1,Nsize+1,npts);
	fprintf(fi, " %14.7e %14.7e %14.7e %14.7e\n",pg->Xmin,pg->Xmax,pg->Zmin,pg->Zmax);
	fprintf(fi, " %14.7e\n",pl->beta);
	fprintf(fi, " %14.7e\n\n",pl->betap);

	/* PsiGrid */
	for (iz = 0; iz <= Nsize; iz++)
		for (ix = 0; ix <= Nsize; ix++)
			fprintf(fi, " %14.7e\n", pg->Psi[ix][iz]);
	fprintf(fi, "\n");

	/* Psi(ipsi) */
	for (i = 0; i < npts; i++) {
		PsiX = i * PsiXmax / (npts - 1);
		Psi = pl->PsiAxis + PsiX * DelPsi;
		fprintf(fi, " %14.7e\n", Psi);
	}
	fprintf(fi, "\n");

	/* F(ipsi) */
	for (i = 0; i < npts; i++) {
		PsiX = i * PsiXmax / (npts - 1);
		Psi = pl->PsiAxis + PsiX * DelPsi;
		if (pl->ModelType == Plasma_Std) {
			G = 1.0 - DelPsi * pl->G2p[1] * pow(1.0 - PsiX, pl->StndG) / pl->StndG;
		} else if (pl->ModelType == Plasma_IsoNoFlow) {
			G = fpoly_int(pl->G2p, PsiX, pl->G2pTerms, DelPsi, 1.0);
		}
		G = sqrt(G)*pl->B0*pl->R0;
		fprintf(fi, " %14.7e\n", G);
	}
	fprintf(fi, "\n");

	/* P(ipsi) */
	for (i = 0; i < npts; i++) {
		PsiX = i * PsiXmax / (npts - 1);
		Psi = pl->PsiAxis + PsiX * DelPsi;
		if (pl->ModelType == Plasma_Std) {
			P = -DelPsi * pl->Pp[1] * pow(1.0 - PsiX, pl->StndP) / pl->StndP;
		} else if (pl->ModelType == Plasma_IsoNoFlow) {
			P = fpoly_int(pl->Pp, PsiX, pl->PpTerms, DelPsi, P_EDGE);
		}
		P = P / MU0;
		fprintf(fi, " %14.7e\n", P);
	}
	fprintf(fi, "\n");

	/* q(ipsi) */
	for (i = 0; i < npts; i++) {
		fprintf(fi, " %14.7e\n", pl->q_pr[i]);
	}
	fprintf(fi, "\n");

	fclose(fi);
}

void GS2Output(TOKAMAK * td)
{
        FILE         *fi;
	PSIGRID      *pg;
	PLASMA       *pl;
	double        Psi, DelPsi, P;
	double        PsiFCFS, PsiLCFS, PsiNorm;
	int           i,j,k,q,np,npts;
	int	      ix, iz, Nsize;
	char          fname[32] = "";

	strncat(fname, td->Oname, 20);	/* take 1st 20 characters */
	strcat(fname, "_gs2.out");

	fi = fopen(fname, "w");
	if (!fi)
		nrerror("ERROR:	Could not open file for writing in GS2Output.");

	pg = td->PsiGrid;
	pl = td->Plasma;
	npts = pl->NumPsiPts;
	DelPsi = pl->PsiLim - pl->PsiAxis;
	PsiFCFS = pl->PsiAxis;
        PsiLCFS = pl->PsiLim;
        PsiNorm = 1.0;  /* what should this be?  2*PI*Current... what for? */

	Nsize = pg->Nsize;


	/* Grid Size */
        /* write(7,200) nr, nz, nrho */
	fprintf(fi, " %6d %6d %6d\n",Nsize+1,Nsize+1,npts);
        /* write(7,210) R_max, Z_max */
        fprintf(fi, " %14.7e %14.7e\n",pg->Xmax,pg->Zmax);
        /* Rmin, Zmin */
        /* write(7,210) 1., 0. */
        fprintf(fi, " %14.7e %14.7e\n",pg->Xmin,pg->Zmin);
        /* Psi_FCFS, Psi_LCFS, not normalized and 2*pi*current */
        fprintf(fi, " %14.7e %14.7e %14.7e\n", PsiFCFS, PsiLCFS, PsiNorm);
        /* Pressure profile p(psi_norm)  n_psi_norm = nrho ?? or nr ??? */
        /* priyanka has this as nr... but with no good reason. */
        np = Nsize+1;
        for (i=0, k=0; i<np; i++) {
            Psi = PsiFCFS + i*DelPsi/(np-1);
            P = PlasmaP(pl,Psi);
            /* Calculate p of psi at normalized rho point */
            fprintf(fi, " %16.8e", P);
            if (!(++k % 5)) fprintf(fi,"\n");
        }
        if (k % 5) fprintf(fi,"\n");
        /* Psi (norm?) profile as function of R & Z */
        /* careful with the i, j direction... */
        /*   write(7,210) ((Psi(i,j), i=1,nr), j=1,nz) */
        for (j=0, k=0; j<=Nsize; j++)
            for(i=0; i<=Nsize; i++) {
            fprintf(fi, " %16.8e", pg->Psi[i][j]);
            if (!(++k % 5)) fprintf(fi,"\n");
        }
        if (k % 5) fprintf(fi,"\n");

        /* write(7,210) (rho_mid(i), i=1,nrho) */

        /* write(7,210) (psi_mid(i), i=1,nrho) */
	/*  100 format (4(1x,e12.5))
            200 format (5(1x,i5))
            210 format (5(e16.8)) */

	fclose(fi);


}
