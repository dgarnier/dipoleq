/*
** TokaMac v2.0
**
** psigrid.c
**
** This file contains routines which solve the grad-shafranov
** equation for a given plasma current.
**
** File:		psigrid.c
** Date:		June 19, 1992
**
** Revisions:
**
**		August 3, 1993		Added GetIsPlasma
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdlib.h>
#include "VAX_Alloc.h"
#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "interpolate.h"
#include "psigrid.h"
#include "multitask.h"

#ifdef THINK_C
#pragma options(!assign_registers)
#endif

#define PI          		3.14159265358979323
#define MU0					12.56637061E-7
#define TWOPI				6.283185307
#define DEG_RAD 			0.01745329252

#define NUM_RELAXATIONS 		3
#define GOPDE_MAX_ITERATIONS 	6
#define MAX_RESTHRESHOLD		1.0e-4

#ifndef IDL
extern FILE  *LogFile;
#endif

/***************************************************************
 *	The Private Functions Within This File
 ***************************************************************/

void          DoRelaxation(PSIGRID *);
void          NewSolution(PSIGRID *);
void          DoNewTwoSolution(PSIGRID *);
void          MakePsiSymmetric(PSIGRID *);

void          InitializeFromParent(PSIGRID *, PSIGRID *);
void          NewMSolution(PSIGRID *);
void          RefineSol(PSIGRID *, PSIGRID *);

/***************************************************************
 * new_PsiGrid
 *
 * allocates an empty PSIGRID structure
 ***************************************************************/
PSIGRID      *new_PsiGrid()
{

	PSIGRID      *pg;

	pg = (PSIGRID *) malloc((unsigned) sizeof(PSIGRID));
	if (!pg)
		nrerror("ERROR: Allocation error in new_PsiGrid.");

	pg->Psi = NULL;
	pg->Current = NULL;
	pg->Residual = NULL;
	pg->X = NULL;
	pg->Z = NULL;
	pg->IsPlasma = NULL;

	pg->Nsize = 16;				/* defaults */
	pg->Symmetric = UpDownSymmetric;
	pg->MaxRes = 1.0;
	pg->PastMaxRes = 1.0;

	pg->Xmin = 1.0;				/* defaults */
	pg->Xmax = 2.0;
	pg->Zmin = -0.5;
	pg->Zmax = 0.5;

	pg->BoundError = 1.0;
	pg->BoundThreshold = 1.0e-4;
	pg->ResThreshold = 1.0e-4;
	pg->UnderRelax1 = 0.4;
	pg->UnderRelax2 = 0.4;

	pg->PsiAxis = -1.0;
	pg->PsiLim = 0.0;
	pg->DelPsi = pg->PsiLim - pg->PsiAxis;

	return pg;
}

/***************************************************************
 * init_PsiGrid
 *
 * allocates and fills space within a PSIGRID structure
 ***************************************************************/
void          init_PsiGrid(PSIGRID * pg)
{
	int           ix, iz;
	int           nsize;

	nsize = pg->Nsize;

	pg->Psi = dmatrix(0, nsize, 0, nsize);
	pg->Current = dmatrix(0, nsize, 0, nsize);
	pg->Residual = dmatrix(0, nsize, 0, nsize);
	pg->X = dvector(0, nsize);
	pg->Z = dvector(0, nsize);
	pg->IsPlasma = imatrix(0, nsize, 0, nsize);

	pg->dx = (pg->Xmax - pg->Xmin) / nsize;
	pg->dz = (pg->Zmax - pg->Zmin) / nsize;

	for (ix = 0; ix <= nsize; ix++)
		pg->X[ix] = pg->Xmin + ix * pg->dx;
	for (iz = 0; iz <= nsize; iz++)
		pg->Z[iz] = pg->Zmin + iz * pg->dz;

	for (ix = 0; ix <= nsize; ix++)
		for (iz = 0; iz <= nsize; iz++) {
			pg->IsPlasma[ix][iz] = 0;
			pg->Psi[ix][iz] = 0.0;
			pg->Current[ix][iz] = 0.0;
			pg->Residual[ix][iz] = 0.0;
		}
}

/***************************************************************
 * free_PsiGrid
 *
 *
 ***************************************************************/
void          free_PsiGrid(PSIGRID * pg)
{
	int           n;

	n = pg->Nsize;

	free_dvector(pg->X, 0, n);
	free_dvector(pg->Z, 0, n);
	free_imatrix(pg->IsPlasma, 0, n, 0, n);
	free_dmatrix(pg->Psi, 0, n, 0, n);
	free_dmatrix(pg->Current, 0, n, 0, n);
	free_dmatrix(pg->Residual, 0, n, 0, n);

	free(pg);
}

/***************************************************************
 * GetPsi
 *
 * This function interpolates/extrapolates PsiGrid in order to
 * return the value of Psi at the requested point, (x, z).
 *
 * We use a simple bilinear interpolation.
 *
 ***************************************************************/
double        GetPsi(PSIGRID * pg, double x, double z)
{
	return interpolate(pg, pg->Psi, x, z);
}

/***************************************************************
 * GetIsPlasma
 *
 * This function interpolates/extrapolates IsPlasma in order to
 * return a value between 0 and 1 at the requested point, (x, z).
 *
 * We use a simple bilinear interpolation.
 *
 ***************************************************************/
double        GetIsPlasma(PSIGRID * pg, double x, double z)
{
	return interpolate_int(pg, pg->IsPlasma, x, z);
}

/***************************************************************
|# Method: GetNewResidual
|#
|# The residual is the difference between
|#
|# Æ*^2(Psi) = + 2pi x ( J(x,z) + Res(x,z) ).
|# or
|# Æ*^2(Psi)/2piX - J(x,z) = Res(x,z).
|#
***************************************************************/
void          GetNewResidual(PSIGRID * pg)
{
	double        dx2, dz2;
	double      **P;
	double      **J;
	double      **R;
	double       *x;

	int           ix, iz, NMax;
	double        MaxCurDen = 0.0;
	double        temp1, temp2;
	double        del2;

    MULTI;

	pg->PastMaxRes = pg->MaxRes;
	pg->MaxRes = 0.0;

	dx2 = pg->dx * pg->dx;
	dz2 = pg->dz * pg->dz;

	P = pg->Psi;
	J = pg->Current;
	R = pg->Residual;
	x = pg->X;

	NMax = pg->Nsize;

	for (ix = 0; ix <= NMax; ix++)
		for (iz = 0; iz <= NMax; iz++)
			R[ix][iz] = 0.0;

	for (ix = 1; ix < NMax; ix++)
		for (iz = 1; iz < NMax; iz++)
			if (fabs(J[ix][iz]) > MaxCurDen)
				MaxCurDen = fabs(J[ix][iz]);

	for (ix = 1; ix < NMax; ix++) {
		temp1 = 1.0 / (2.0 * pg->dx * x[ix]);
		temp2 = TWOPI * x[ix];
		for (iz = 1; iz < NMax; iz++) {
			del2 = (P[ix + 1][iz] - 2.0 * P[ix][iz] + P[ix - 1][iz]) / dx2;
			del2 -= (P[ix + 1][iz] - P[ix - 1][iz]) * temp1;
			del2 += (P[ix][iz + 1] - 2.0 * P[ix][iz] + P[ix][iz - 1]) / dz2;
			R[ix][iz] = del2 / temp2 - J[ix][iz];
			if (pg->MaxRes < fabs(R[ix][iz]))
				pg->MaxRes = fabs(R[ix][iz]);
		}
	}
	pg->MaxRes = pg->MaxRes / MaxCurDen;	/* Normalize the residual.*/
}

/***************************************************************
 * MakePsiSymmetric
 *
 *
 ***************************************************************/
void          MakePsiSymmetric(PSIGRID * pg)
{
	int           ix, iz;
	int           NMax;
	double      **P;

	P = pg->Psi;
	NMax = pg->Nsize;
    MULTI;
	for (ix = 0; ix <= NMax; ix++)
		for (iz = 0; iz < (NMax / 2); iz++) {
			P[ix][iz] = 0.5 * (P[ix][iz] + P[ix][NMax - iz]);
			P[ix][NMax - iz] = P[ix][iz];
		}
}

/***************************************************************
|# Method: DoRelaxation
|#
|#	RED-BLACK GAUSS SEIDEL
|#
|#	This uses an algorithm that USES A "RED" "Black" relaxation technique
|#	which is much faster than the tridiagonal technique shown above.	This
|#	really only works with Multi-Grid since it is terrible for long-wavelength
|#	relaxation...but great for the short.
|#
***************************************************************/
void          DoRelaxation(PSIGRID * pg)
{
	int           ix, iz;
	double        Temp;
	double      **p, **j, *x;
	double        TwoPidx2, dxdz2, dxo2, bx;

    MULTI;

	p = pg->Psi;				/* p is a pointer to the poloidal flux */
	j = pg->Current;			/* J is the toroidal current */
	x = pg->X;					/* X is a vector of major radial coordinate */

	TwoPidx2 = TWOPI * pg->dx * pg->dx;
	dxdz2 = pg->dx / pg->dz;
	dxdz2 = dxdz2 * dxdz2;
	dxo2 = pg->dx / 2.0;
	bx = -2.0 * (1.0 + dxdz2);

	/* RED SWEEP...*/
	for (iz = 1; iz < pg->Nsize; iz++) {
		ix = (iz % 2) ? 1 : 2;	/* ix = 2, if iz is even */
		do {
			Temp = TwoPidx2 * x[ix] * j[ix][iz]
				+ (dxo2 / x[ix]) * (p[ix + 1][iz] - p[ix - 1][iz])
				- (p[ix + 1][iz] + p[ix - 1][iz])
				- dxdz2 * (p[ix][iz + 1] + p[ix][iz - 1]);
			p[ix][iz] = Temp / bx;
			ix += 2;
		} while (ix < pg->Nsize);
	};
	/* BLACK SWEEP... */
	for (iz = 1; iz < pg->Nsize; iz++) {
		ix = (iz % 2) ? 2 : 1;
		do {
			Temp = TwoPidx2 * x[ix] * j[ix][iz]
				+ (dxo2 / x[ix]) * (p[ix + 1][iz] - p[ix - 1][iz])
				- (p[ix + 1][iz] + p[ix - 1][iz])
				- dxdz2 * (p[ix][iz + 1] + p[ix][iz - 1]);
			p[ix][iz] = Temp / bx;
			ix += 2;
		} while (ix < pg->Nsize);
	};
}

/***************************************************************
|# Method: InitializeFromParent
|#
|# Here, -Res(x,y) is used as the source for the smaller grid.
|#
|# I have added "2D Full Weighting" of the residual into the coarser
|# grid.	This should work better especially with the Red-Black Gauss-Seidel
|# solution procedure.
|#
|# mg is the multigrid, and
|# pg is the parent.
|#
***************************************************************/
void          InitializeFromParent(PSIGRID * mg, PSIGRID * pg)
{
	int           ix, iz, ix2, iz2;
	int           SmallerN;
	double      **R;

	R = pg->Residual;

	SmallerN = mg->Nsize;
	for (ix = 0; ix <= SmallerN; ix++)
		for (iz = 0; iz <= SmallerN; iz++)
			mg->Current[ix][iz] = 0.0;
    MULTI;
	for (ix = 1; ix < SmallerN; ix++)
		for (iz = 1; iz < SmallerN; iz++) {
			ix2 = 2 * ix;
			iz2 = 2 * iz;
			mg->Current[ix][iz] = R[ix2 - 1][iz2 - 1] + R[ix2 - 1][iz2 + 1] + R[ix2 + 1][iz2 - 1] + R[ix2 + 1][iz2 + 1];
			mg->Current[ix][iz] += 2.0 * (R[ix2][iz2 - 1] + R[ix2][iz2 + 1] + R[ix2 - 1][iz2] + R[ix2 + 1][iz2]);
			mg->Current[ix][iz] += 4.0 * R[ix2][iz2];
			mg->Current[ix][iz] = -mg->Current[ix][iz] / 16.0;
		}
}

/***************************************************************
|# Method: RefineSol
|#
|# Refine the solution to a finer grid using linear interpolation.
|#
|# mg is the multigrid and pg is the parent psigrid.
|#
***************************************************************/
void          RefineSol(PSIGRID * mg, PSIGRID * pg)
{
	int           ix, iz;		/* indicies for mg's, smaller Solution */
	int           nx, nz;		/* indicies for pg's, bigger Solution */
	int           xodd, zodd;
	double        add;
	double      **mgPsi;
	double      **pgPsi;

	mgPsi = mg->Psi;
	pgPsi = pg->Psi;

    MULTI;

	for (nx = 1; nx < pg->Nsize; nx++)
		for (nz = 1; nz < pg->Nsize; nz++) {
			ix = nx / 2;
			iz = nz / 2;
			xodd = nx % 2;		/* xodd == 1, if odd */
			zodd = nz % 2;		/* zodd == 1, if odd */
			if (!xodd && !zodd)
				add = mgPsi[ix][iz];
			if (!xodd && zodd)
				add = 0.5 * (mgPsi[ix][iz] + mgPsi[ix][iz + 1]);
			if (xodd && !zodd)
				add = 0.5 * (mgPsi[ix][iz] + mgPsi[ix + 1][iz]);
			if (xodd && zodd)
				add = 0.25 * (mgPsi[ix][iz] + mgPsi[ix + 1][iz] +
							  mgPsi[ix][iz + 1] + mgPsi[ix + 1][iz + 1]);
			pgPsi[nx][nz] += add;
		}

}

/***************************************************************
|# Method: DoNewTwoSolution
|#
|# This method performs a complete multi-grid iteration for a 2x2 grid.
|# It should be called for a MGrid which is 4x4. This routine will
|#
|#  1. Prepare the course grid
|#  2. Solve the 2x2
|#  3. Refine the calling 4x4
|#
|#
***************************************************************/
void          DoNewTwoSolution(PSIGRID * pg)
{
	double      **P;
	double        Psi22, Cur22;
	int           nx, nz;
	double        bx;			/* -2(1+(dx/dz)^2) and -2(1+(dz/dx)^2) */
	double        dxdz2;		/* (dx/dz)^2 */
	int           xodd, zodd;
	double        add;

	/* The following are the variables of the reduced 2x2 grid...*/
	dxdz2 = pg->dx / pg->dz;
	dxdz2 = dxdz2 * dxdz2;

	P = pg->Psi;

	/* Using the residual to find the current for the 2x2 grid...*/
	Cur22 = -pg->Residual[1 * 2][1 * 2];

	/* Finding the solution to the 2x2 grid...*/
	bx = -2.0 * (1.0 + dxdz2);
	Psi22 = TWOPI * (pg->dx * pg->dx / 4.0) * pg->X[1 * 2] * Cur22;
	Psi22 = Psi22 / bx;

	/* Refining the 2x2 back up to the 4x4...*/
	for (nx = 1; nx < pg->Nsize; nx++)
		for (nz = 1; nz < pg->Nsize; nz++) {
			xodd = nx % 2;		/* if (xodd==1), then xodd is odd */
			zodd = nz % 2;		/* if (zodd==1), then zodd is odd */
			if (!xodd && !zodd)
				add = Psi22;
			if ((!xodd && zodd) || (xodd && !zodd))
				add = 0.5 * Psi22;
			if (xodd && zodd)
				add = 0.25 * Psi22;
			P[nx][nz] += add;
		}
}

/***************************************************************
|# Method: NewMSolution
|#
|# This solution technique recursively calls courser grids to reduce
|# the size of the residual.
|# It uses a simple smoothing routine to create the courser grids.
|#
|# We choose to do 2 ADI solutions per MULTIGRID solution.
|# We check for NMax = 4, so that we can do the simplier problem...
|#
|# **NOTE** GetNewResidual must always be called before NewMSolution.
|#
***************************************************************/
void          NewMSolution(PSIGRID * pg)
{
	PSIGRID      *aMGrid;

    MULTI;
	/** 1. Find residual to be used as a source for the new grid...**/
	/**    This should have been called already.                   **/
	/** GetNewResidual(pg);                                        **/

	/* CHECK if NMax = 4...*/
	switch (pg->Nsize) {
	  case 4:
		  /* Do complete 2x2 multi-grid solution*/
		  DoNewTwoSolution(pg);
		  break;
	  default:
		  /***2. Set up the smaller grid...*/
		  aMGrid = new_PsiGrid();
		  aMGrid->Nsize = pg->Nsize / 2;
		  aMGrid->Xmin = pg->Xmin;
		  aMGrid->Xmax = pg->Xmax;
		  aMGrid->Zmin = pg->Zmin;
		  aMGrid->Zmax = pg->Zmax;
		  init_PsiGrid(aMGrid);
		  InitializeFromParent(aMGrid, pg);
		  break;
	}

	/** 3. Remove residual, since we no longer need it.*/

	/** 4. Solve the auxillary problem...*/
	switch (pg->Nsize) {
	  case 4:
		  /* Do nothing */
		  break;
	  default:
		  /* a "W" cycle...*/
		  NewSolution(aMGrid);
		  GetNewResidual(aMGrid);
		  NewMSolution(aMGrid);	/* A recursive call! */
		  GetNewResidual(aMGrid);
		  NewMSolution(aMGrid);	/* A recursive call! */
		  break;
	}

	if (pg->Nsize > 4) {
		/** 5. Refine the solution back up to the larger grid...*/
		RefineSol(aMGrid, pg);
		free_PsiGrid(aMGrid);
	}
	/** 6...and try again.*/
	NewSolution(pg);			/* Do we need to do this if pg->Nsize == 4 ? */
}

/***************************************************************
 * NewSolution
 ***************************************************************/
void          NewSolution(PSIGRID * pg)
{
	int           i;
	for (i = 1; i <= NUM_RELAXATIONS; i++)
		DoRelaxation(pg);
}

/*
**
**	GoPDE
**
*/
void          GoPDE(PSIGRID * pg)
{
	int           i = 0;
	double        ResThr;

	ResThr = sqrt(pg->ResThreshold * pg->BoundError);

#ifdef IDL
    IDL_Message(IDL_M_GENERIC,IDL_MSG_INFO, "INFO:	Solving Grad-Shafranov Equation.\n");
#else
	printf("INFO:	Solving Grad-Shafranov Equation.\n");
	fprintf(LogFile, "INFO:	Solving Grad-Shafranov Equation.\n");
#endif

	ResThr = sqrt(pg->ResThreshold * pg->BoundError);
	if (ResThr > MAX_RESTHRESHOLD)
		ResThr = MAX_RESTHRESHOLD;

	NewSolution(pg);
	GetNewResidual(pg);
	while ((i < GOPDE_MAX_ITERATIONS) && (pg->MaxRes > ResThr)) {
		NewMSolution(pg);
		if (pg->Symmetric)
			MakePsiSymmetric(pg);
		GetNewResidual(pg);
		i++;
	}
	NewSolution(pg);
	if (pg->Symmetric)
		MakePsiSymmetric(pg);
	GetNewResidual(pg);

#ifndef IDL
	printf("		[After %d iterations, MaxRes=%g]\n", i, pg->MaxRes);
	fprintf(LogFile, "		[After %d iterations, MaxRes=%g]\n", i, pg->MaxRes);
#endif
}
