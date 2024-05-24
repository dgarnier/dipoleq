/*
** TokaMac v2.0
**
** GetPlasmaParameters.c
**
**
**
** File:		GetPlasmaParameters.c
** Date:		March 31, 1993
**
** Revisions:
**
**		Sept. 16, 1993		Added virial terms (from Lao, Nuc Fusion, 1985)
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdlib.h>
#include "VAX_Alloc.h"
#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "contour.h"
#include "fpoly.h"
#include "psigrid.h"
#include "plasma.h"
#include "tokamak.h"
#include "GetFluxParameters.h"
#include "GetPlasmaParameters.h"
#include "CPlasmaModel.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

#define P_EDGE		0.0


extern "C" FILE  *LogFile;

/*
** Local variables
**
*/
double        gPerimeter;		/* a variable to store the perimeter during contour trace */

/*
** Local prototypes
**
*/
void          GetPParm_Std(TOKAMAK *);
void          GetPParm_IsoNoFlow(TOKAMAK *);
void          GetPParm_IsoFlow(TOKAMAK *);
void          GetPParm_AnisoNoFlow(TOKAMAK *);
void          GetPParm_AnisoFlow(TOKAMAK *);
void          TracePerimeter(double x, double z, double p, int flag);
void          GetGeometry(TOKAMAK * td);
void          GetBeta(TOKAMAK *td);
void          GetVirial_Vol(TOKAMAK * td);

/*
** G E T G R A D P S I
**
*/
void          GetGradPsi(TOKAMAK * td)
{
	PLASMA       *pl;
	PSIGRID      *pg;
	int  		n, m1, m2, ix, iz;
	double      **Psi, **gPsiX, **gPsiZ, **gPsi2;
	double        dx, dz;

	pg = td->PsiGrid;
	pl = td->Plasma;
	n = pg->Nsize;
	dx = pg->dx;
	dz = pg->dz;
	Psi = pg->Psi;
	gPsiX = pl->GradPsiX;
	gPsiZ = pl->GradPsiZ;
	gPsi2 = pl->GradPsi2;

	m1 = n-1;
	m2 = n-2;

	/*  G R A D  P S I,  E T C . . . .  */
	/* DTG 2/18/99 - changed definitions of edge derivatives to be good to 2nd order */
	for (ix = 1; ix < n; ix++)
		for (iz = 1; iz < n; iz++) {  /* 1st deriv, 2nd order */
			gPsiX[ix][iz] = (Psi[ix + 1][iz] - Psi[ix - 1][iz]) / 2.0 / dx;
			gPsiZ[ix][iz] = (Psi[ix][iz + 1] - Psi[ix][iz - 1]) / 2.0 / dz;
		}
	for (ix = 1; ix < n; ix++) {	/* along top and bottom */
		gPsiX[ix][0] = (Psi[ix + 1][0] - Psi[ix - 1][0]) / 2.0 / dx;
		gPsiX[ix][n] = (Psi[ix + 1][n] - Psi[ix - 1][n]) / 2.0 / dx;
		gPsiZ[ix][0] = (-3 * Psi[ix][0] + 4 * Psi[ix][ 1] - Psi[ix][ 2]) / 2.0 / dz;
		gPsiZ[ix][n] = ( 3 * Psi[ix][n] - 4 * Psi[ix][m1] + Psi[ix][m2]) / 2.0 / dz;
	}
	for (iz = 1; iz < n; iz++) {	/* along inside and outside */
		gPsiX[0][iz] = (-3 * Psi[0][iz] + 4 * Psi[ 1][iz] - Psi[ 2][iz]) / 2.0 / dx;
		gPsiX[n][iz] = ( 3 * Psi[n][iz] - 4 * Psi[m1][iz] + Psi[m2][iz]) / 2.0 / dx;
		gPsiZ[0][iz] = (Psi[0][iz + 1] - Psi[0][iz - 1]) / 2.0 / dz;
		gPsiZ[n][iz] = (Psi[n][iz + 1] - Psi[n][iz - 1]) / 2.0 / dz;
	}
	/* Corners */
	gPsiX[0][0] = (-3 * Psi[0][0] + 4 * Psi[ 1][ 0] - Psi[ 2][ 0]) / 2.0 / dx;
	gPsiZ[0][0] = (-3 * Psi[0][0] + 4 * Psi[ 0][ 1] - Psi[ 0][ 2]) / 2.0 / dz;
	gPsiX[n][0] = ( 3 * Psi[n][0] - 4 * Psi[m1][ 0] + Psi[m2][ 0]) / 2.0 / dx;
	gPsiZ[n][0] = (-3 * Psi[n][0] + 4 * Psi[ n][ 1] - Psi[ n][ 2]) / 2.0 / dz;
	gPsiX[0][n] = (-3 * Psi[0][n] + 4 * Psi[ 1][ n] - Psi[ 2][ n]) / 2.0 / dx;
	gPsiZ[0][n] = ( 3 * Psi[0][n] - 4 * Psi[ 0][m1] + Psi[ 0][m2]) / 2.0 / dz;
	gPsiX[n][n] = ( 3 * Psi[n][n] - 4 * Psi[m1][ n] + Psi[m2][ n]) / 2.0 / dx;
	gPsiZ[n][n] = ( 3 * Psi[n][n] - 4 * Psi[ n][m1] + Psi[ n][m2]) / 2.0 / dz;
	/* GradPsi2 */
	for (ix = 0; ix <= n; ix++)
		for (iz = 0; iz <= n; iz++)
			gPsi2[ix][iz] = gPsiX[ix][iz] * gPsiX[ix][iz] + gPsiZ[ix][iz] * gPsiZ[ix][iz];
}

/*
**	TracePerimeter
*/
void          TracePerimeter(double x, double z, double , int flag)
{
	static double lastx, lastz;
	double        ds;

	switch (flag) {
	  case CONTOUR_START:
		  gPerimeter = 0.0;
		  break;
	  case CONTOUR_TRACE:
		  ds = (x - lastx) * (x - lastx) + (z - lastz) * (z - lastz);
		  ds = sqrt(ds);
		  gPerimeter += ds;
		  break;
	}
	lastx = x;
	lastz = z;
}

/*
** G E T G E O M E T R Y
**
*/
void          GetGeometry(TOKAMAK * td)
{
	PLASMA       *pl;
	PSIGRID      *pg;
	int           nmax, ix, iz, i1, i2, i3;
	double      **Psi, **J;
	double       *X, *Z, dx, dz, sum1, sum2;

	pg = td->PsiGrid;
	pl = td->Plasma;
	nmax = pg->Nsize;
	dx = pg->dx;
	dz = pg->dz;
	Psi = pg->Psi;
	J = pg->Current;
	X = pg->X;
	Z = pg->Z;

	/* F I N D   M I N O R   H A L F W I D T H */
	iz = (int) floor((pl->ZMagAxis - pg->Zmin) / dz);	/* iz index for magnetic axis */
	i1 = (int) floor((pl->XMagAxis - pg->Xmin) / dx);	/* ix index for magnetic axis */
	for (ix = i1; ix <= nmax; ix++)
		if (Psi[ix][iz] >= pg->PsiLim) {
			i2 = ix;
			break;
		}
	for (ix = i1; ix >= 0; ix--)
		if (Psi[ix][iz] >= pg->PsiLim) {
			i3 = ix;
			break;
		}
	pl->HalfWidth = (X[i2] - X[i3]) / 2.0;
	pl->RCenter = X[i3] + pl->HalfWidth;	/* Center of outer flux surface */

	/* C U R R E N T   C E N T R O I D,   V O L U M E,  &  A R E A */
	sum1 = sum2 = 0.0;
	i1 = 0;
	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++) {
			sum1 += X[ix] * J[ix][iz];
			if (pg->IsPlasma[ix][iz]) {
				sum2 += X[ix];
				i1++;
			}
		}
	pl->CrossSection = i1 * dx * dz;
	pl->Volume = TWOPI * dx * dz * sum2;
	pl->RSurfaceAvg = sum2 / i1;
	pl->RCentroid = sum1 * dx * dz / pl->Ip / MU0;
	pl->Elongation = pl->CrossSection / PI / pl->HalfWidth / pl->HalfWidth;

	/* P E R I M E T E R    A N D    R S T A R  */
	contour(X, Z, Psi, 0, nmax, 0, nmax, pg->PsiLim, CONTOUR_ONLY_CLOSED, CONTOUR_NO_MIDPOINT, TracePerimeter);
	pl->Perimeter = gPerimeter;
	pl->RStar = 2.0 * pl->Volume / gPerimeter / gPerimeter;

}

/*
** G E T B E T A
**
*/
void          GetBeta(TOKAMAK * td)
{
	PLASMA       *pl;
	PSIGRID      *pg;
	int           nmax, ix, iz, i1;
	double      **P;
	double       *X, *Z, dx, dz, Pavg, B_R0_Vac, Wnorm, Wpol, dWtor, diamag;

	pg = td->PsiGrid;
	pl = td->Plasma;
	nmax = pg->Nsize;
	dx = pg->dx;
	dz = pg->dz;
	X = pg->X;
	Z = pg->Z;

	/* F I N D   V O L U M N   I N T E G R A L   O F   P R E S S U R E */
	switch (pl->ModelType) {
	  default:
		  P = pl->Piso;
		  Pavg = 0.0;
		  for (ix = 1; ix < nmax; ix++)
			  for (iz = 1; iz < nmax; iz++)
				  Pavg += X[ix] * P[ix][iz];
		  Pavg = TWOPI * Pavg * dx * dz;
		  break;
	  case Plasma_AnisoNoFlow:

		  break;
	  case Plasma_AnisoFlow:

		  break;
	}

	pl->TotKinEnergy = 1.5 * Pavg / MU0;

	/* V O L   A V G   B E T A */
	i1 = (int) floor((pl->RCenter - pg->Xmin) / dx);
	B_R0_Vac = pl->B0R0 / X[i1];
	pl->beta = Pavg / MU0 / pl->Volume / (B_R0_Vac * B_R0_Vac / 2.0 / MU0);

	/* P O L O I D A L    B E T A */
	Wnorm = 0.25 * MU0 * pl->RStar * pl->Ip * pl->Ip;
	pl->betap = Pavg / Wnorm / MU0;

	/* I N T   I N D U C T A N C E   &   M U */
	diamag = dWtor = Wpol = 0.0;
	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (pg->IsPlasma[ix][iz]) {
				Wpol += pl->GradPsi2[ix][iz] / X[ix];
				dWtor += X[ix] * (DSQR(pl->B0R0 / X[ix]));
				dWtor -= X[ix] * (DSQR(pl->Bt[ix][iz]));
				diamag += pl->Bt[ix][iz] - pl->B0R0 / X[ix];
			}
	Wpol = Wpol * dx * dz / 2.0 / MU0 / TWOPI;
	pl->TotMagEnergy = Wpol;
	dWtor = TWOPI * dWtor * dx * dz / 2.0 / MU0;
	pl->Diamag = diamag * dx * dz;
	pl->li = Wpol / Wnorm;
	pl->mu = dWtor / Wnorm;
}

/*
**	G E T V I R I A L _ V O L
**
*/
void          GetVirial_Vol(TOKAMAK * td)
{
	PLASMA       *pl;
	PSIGRID      *pg;
	int           nmax, ix, iz;
	int         **ip;
	double        g, sum1 = 0.0, sum2 = 0.0;
	double       *X, *Z, dx, dz;

	pg = td->PsiGrid;
	pl = td->Plasma;
	nmax = pg->Nsize;
	dx = pg->dx;
	dz = pg->dz;
	X = pg->X;
	Z = pg->Z;
	ip = pg->IsPlasma;

	/* F I N D   R _ V R   */
	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (ip[ix][iz] == YesPlasma) {
				g = pl->B0R0 * pl->G[ix][iz] / X[ix];
				g = 0.5 * (pl->B2[ix][iz] - 2.0 * pl->Bt[ix][iz] * pl->Bt[ix][iz] + g * g);
				switch (pl->ModelType) {
				  case Plasma_IsoFlow:
					  g = pl->Piso[ix][iz] + g;	/* not complete */
					  break;
				  case Plasma_AnisoNoFlow:
					  g = pl->Ppar[ix][iz] + g;	/* not complete */
					  break;
				  case Plasma_AnisoFlow:
					  g = pl->Ppar[ix][iz] + g;	/* not complete */
					  break;
				  default :
				  	  g = pl->Piso[ix][iz] + g;
				}
				sum1 += X[ix] * g;
				sum2 += g;
			}
	pl->R_vr = sum1 / sum2;

	/* F I N D    A L P H A _ V R  */
	sum1 = sum2 = 0.0;
	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (ip[ix][iz] == YesPlasma) {
				sum1 += pl->GradPsiX[ix][iz] * pl->GradPsiX[ix][iz];
				sum2 += pl->GradPsi2[ix][iz];
			}
	pl->Alpha_vr = 2.0 * sum1 / sum2;

}

/*
** G E T   P L A S M A   P A R A M E T E R S
**
*/
void          GetPlasmaParameters(TOKAMAK * td)
{
	PLASMA       *pl;

	printf("INFO:	Getting plasma parameters.\n");
	fprintf(LogFile, "INFO:	Getting plasma parameters.\n");

	pl = td->Plasma;


	printf("		[Gradient Psi]\n");
	fprintf(LogFile, "		[Gradient Psi]\n");

	GetGradPsi(td);

	printf("		[Plasma Geometry]\n");
	fprintf(LogFile, "		[Plasma Geometry]\n");

	GetGeometry(td);

	printf("		[Pressure and magnetic field]\n");
	fprintf(LogFile, "		[Pressure and magnetic field]\n");

	switch (pl->ModelType) {
	  case Plasma_Std:
		  GetPParm_Std(td);
		  break;
	  case Plasma_IsoNoFlow:
		  GetPParm_IsoNoFlow(td);
		  break;
	  case Plasma_IsoFlow:
		  GetPParm_IsoFlow(td);
		  break;
	  case Plasma_AnisoNoFlow:
		  GetPParm_AnisoNoFlow(td);
		  break;
	  case Plasma_AnisoFlow:
		  GetPParm_AnisoFlow(td);
		  break;
	 default:
	 	  pl->Model->GetPParam(td);
	}

	printf("		[Plasma beta]\n");
	fprintf(LogFile, "		[Plasma beta]\n");

	GetBeta(td);

	printf("		[Virial Integrals]\n");
	fprintf(LogFile, "		[Virial Integrals]\n");

	GetVirial_Vol(td);

 	printf("		[Flux profile parameters]\n");
	fprintf(LogFile, "		[Flux profile parameters]\n");

	GetFluxParameters(td);
}

/*
** 	G E T P P A R M _ S T D
**
**	Compute:
**			P		Pressure
**			G		Normalized toroidal flux
**			B2		Square mod-B
*/
void          GetPParm_Std(TOKAMAK * td)
{
	PLASMA       *pl;
	PSIGRID      *pg;
	int           nmax, ix, iz;
	double      **Psi, **P, **gPsiX, **gPsiZ, **gPsi2, **G, **Bt, **B2;
	double        dx, dz, DelPsi, PsiAxis, PsiX;
	int         **ip;

	pg = td->PsiGrid;
	pl = td->Plasma;
	nmax = pg->Nsize;
	dx = pg->dx;
	dz = pg->dz;
	DelPsi = pg->DelPsi;
	PsiAxis = pg->PsiAxis;
	Psi = pg->Psi;
	ip = pg->IsPlasma;
	P = pl->Piso;
	gPsiX = pl->GradPsiX;
	gPsiZ = pl->GradPsiZ;
	gPsi2 = pl->GradPsi2;
	G = pl->G;
	Bt = pl->Bt;
	B2 = pl->B2;

	/*  P R E S S U R E,  G,  E T C . . . .  */
	for (ix = 0; ix <= nmax; ix++)
		for (iz = 0; iz <= nmax; iz++) {
			if (ip[ix][iz]) {
				PsiX = (Psi[ix][iz] - PsiAxis) / DelPsi;
				P[ix][iz] = -DelPsi * (pl->Pp[1] / pl->StndP) * pow(1.0 - PsiX, pl->StndP);
				G[ix][iz] = 1.0 - DelPsi * (pl->G2p[1] / pl->StndG) * pow(1.0 - PsiX, pl->StndG);
				G[ix][iz] = sqrt(G[ix][iz]);
			} else {
				P[ix][iz] = 0.0;
				G[ix][iz] = 1.0;
			}
			Bt[ix][iz] = G[ix][iz] * pl->B0R0 / pg->X[ix];
			B2[ix][iz] = gPsi2[ix][iz] / TWOPI / pg->X[ix] / TWOPI / pg->X[ix];
			B2[ix][iz] = B2[ix][iz] + DSQR(Bt[ix][iz]);
		}
}

/*
** G E T P P A R M _ I S O N O F L O W
**
**	Compute:
**			P		Pressure
**			G		Normalized toroidal flux
**			B2		Square mod-B
*/
void          GetPParm_IsoNoFlow(TOKAMAK * td)
{
	PLASMA       *pl;
	PSIGRID      *pg;
	int           nmax, ix, iz;
	double      **Psi, **P, **gPsiX, **gPsiZ, **gPsi2, **G, **Bt, **B2;
	double        dx, dz, DelPsi, PsiAxis, PsiX;
	int         **ip;

	pg = td->PsiGrid;
	pl = td->Plasma;
	nmax = pg->Nsize;
	dx = pg->dx;
	dz = pg->dz;
	DelPsi = pg->DelPsi;
	PsiAxis = pg->PsiAxis;
	Psi = pg->Psi;
	ip = pg->IsPlasma;
	P = pl->Piso;
	gPsiX = pl->GradPsiX;
	gPsiZ = pl->GradPsiZ;
	gPsi2 = pl->GradPsi2;
	G = pl->G;
	Bt = pl->Bt;
	B2 = pl->B2;

	/*  P R E S S U R E,  G,  E T C . . . .  */
	for (ix = 0; ix <= nmax; ix++)
		for (iz = 0; iz <= nmax; iz++) {
			if (ip[ix][iz]) {
				PsiX = (Psi[ix][iz] - PsiAxis) / DelPsi;
				P[ix][iz] = fpoly_int(pl->Pp, PsiX, pl->PpTerms, DelPsi, P_EDGE);
				G[ix][iz] = fpoly_int(pl->G2p, PsiX, pl->G2pTerms, DelPsi, 1.0);
				G[ix][iz] = sqrt(G[ix][iz]);
			} else {
				P[ix][iz] = 0.0;
				G[ix][iz] = 1.0;
			}
			Bt[ix][iz] = G[ix][iz] * pl->B0R0 / pg->X[ix];
			B2[ix][iz] = gPsi2[ix][iz] / TWOPI / pg->X[ix] / TWOPI / pg->X[ix];
			B2[ix][iz] = B2[ix][iz] + DSQR(Bt[ix][iz]);
		}
}

void          GetPParm_IsoFlow(TOKAMAK *)
{

}

void          GetPParm_AnisoNoFlow(TOKAMAK *)
{

}

void          GetPParm_AnisoFlow(TOKAMAK *)
{

}
