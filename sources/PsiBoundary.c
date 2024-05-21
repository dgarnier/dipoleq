/*
** TokaMac v2.0
**
** PsiBoundary.c
**
** Computes the value of psi around the computational domain
** for an axisymmetric tokamak.
**
** The plasma current and a coil set are computed separately.
**
** NOTE: You must run PlasmaBoundary BEFORE running CoilBoundary.
**
**
** File:		PsiBoundary.c
** Date:		March 20, 1993
**
** Revisions:
**
**		August 5, 1993		Added perfectly conducting shells
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <math.h>
#include <stdlib.h>
#include "VAX_Alloc.h"
#include <stdio.h>
#include "nrutil.h"
#include "coil.h"
#include "shell.h"
#include "psigrid.h"
#include "green.h"
#include "tokgreen.h"
#include "PsiBoundary.h"
#include "multitask.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

extern FILE  *LogFile;

#ifdef __cplusplus
extern "C" {
#endif

double        LHLoop(LHARY * lha, PSIGRID * pg, LHVEC * dudn, int ix, int iz);
void          PsiShellBoundary(PSIGRID * pg, SHELL * c);

#ifdef __cplusplus
}
#endif


/*
**
** This performs a Lackner / Von Hagenow integral around the outside of the
** computational mesh. We use Simpson's Rule.
** lhv contains (dx/x)*Green or (dz/x)*Green depending on the side.
**
** The arguments (ix,iz) correspond to a location on the boundary.
**
*/
double        LHLoop(LHARY * lha, PSIGRID * pg, LHVEC * dudn, int ix, int iz)
{
	int           ixp, izp, nmax, izs;
	LHVEC        *lhv, *lhvp;
	double        Sum;
	int           flip = 0;		/* if flip = 1, then switch top & bottom */

	nmax = pg->Nsize;

	/* find the appropriate LHVEC for the given boundary location */
	izs = iz;
	if (iz > nmax / 2) {
		flip = 1;
		izs = nmax - iz;
	}
	if (iz == nmax)				/* (ix,iz) is located on the top */
		lhvp = lha->Bot[ix];
	else if (iz == 0)			/* (ix,iz) is located on the bottom */
		lhvp = lha->Bot[ix];
	else if (ix == 0)			/* (ix,iz) is located on the inside */
		lhvp = lha->In[izs];
	else						/* (ix,iz) is located on the outside */
		lhvp = lha->Out[izs];

	if (!flip)
		lhv = lhvp;
	else {
		lhv = new_LHvec(nmax);
		for (ixp = 0; ixp <= nmax; ixp++) {
			lhv->Top[ixp] = lhvp->Bot[ixp];
			lhv->Bot[ixp] = lhvp->Top[ixp];
		}
		for (izp = 0; izp <= nmax; izp++) {
			lhv->In[izp] = lhvp->In[nmax - izp];
			lhv->Out[izp] = lhvp->Out[nmax - izp];
		}
	}

	/* integrate the LHVEC around the computational domain */
	Sum = 0.0;

	/*Top*/
	ixp = 1;
	do {
		Sum += (lhv->Top[ixp - 1] * dudn->Top[ixp - 1] + 4.0 * lhv->Top[ixp] * dudn->Top[ixp]
				+ lhv->Top[ixp + 1] * dudn->Top[ixp + 1]);
		ixp += 2;
	}
	while (ixp < nmax);

	/*Bottom*/
	ixp = 1;
	do {
		Sum += (lhv->Bot[ixp - 1] * dudn->Bot[ixp - 1] + 4.0 * lhv->Bot[ixp] * dudn->Bot[ixp]
				+ lhv->Bot[ixp + 1] * dudn->Bot[ixp + 1]);
		ixp += 2;
	}
	while (ixp < nmax);

	/*Outside*/
	izp = 1;
	do {
		Sum += (lhv->Out[izp - 1] * dudn->Out[izp - 1] + 4.0 * lhv->Out[izp] * dudn->Out[izp]
				+ lhv->Out[izp + 1] * dudn->Out[izp + 1]);
		izp += 2;
	}
	while (izp < nmax);

	/*inside*/
	izp = 1;
	do {
		Sum += (lhv->In[izp - 1] * dudn->In[izp - 1] + 4.0 * lhv->In[izp] * dudn->In[izp]
				+ lhv->In[izp + 1] * dudn->In[izp + 1]);
		izp += 2;
	}
	while (izp < nmax);

	Sum = Sum / 3.0;

	/***
	 * Eq. 11 of Johnson, et al. is not correct since it doesn't have TWOPI.
	 * With our definition of Green's Function...As Ling and Jardin's Eq. 35,
	 * we don't need an extra TWOPI.
	 ***/

	if (flip)
		free_LHvec(lhv, nmax);
	return Sum;
}

/*
**
**
** Finds psi along outer boundary due to the plasma.
**
**
*/
void          PsiPlasmaBoundary(LHARY * lha, PSIGRID * pg)
{
	int           nmax, ix, iz, i;
	double        dx, dz, dxdz;
	LHVEC        *dudn;
	PSIGRID      *USol;			/*This is a temporary solution for plasma's boundary*/
	double      **Upsi;

	dx = pg->dx;
	dz = pg->dz;
	dxdz = dx * dz;
	nmax = pg->Nsize;

	USol = new_PsiGrid();
	USol->Nsize = nmax;
	USol->Xmin = pg->Xmin;
	USol->Xmax = pg->Xmax;
	USol->Zmin = pg->Zmin;
	USol->Zmax = pg->Zmax;
	USol->Symmetric = pg->Symmetric;
	USol->BoundError = pg->ResThreshold;	/* to set tolerance for solution */
	init_PsiGrid(USol);			/* allocate arrays */

	/* First, some the contribution due to the plasma current...*/
	/* Copy the toroidal current into USol...*/
	/* This only contains the current from the plasma.*/
	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			USol->Current[ix][iz] = pg->Current[ix][iz];

	/* Next, solve the fixed boundary problem...*/
	GoPDE(USol);

	Upsi = USol->Psi;
	dudn = new_LHvec(nmax);
	for (i = 0; i <= nmax; i++) {
		if ((i == 0) || (i == nmax)) {
			dudn->Top[i] = dudn->Bot[i] = 0.0;
			dudn->In[i] = dudn->Out[i] = 0.0;
		} else {
			dudn->Top[i] = (3.0 * Upsi[i][nmax] - 4.0 * Upsi[i][nmax - 1] + Upsi[i][nmax - 2]) / (2.0 * dz);
			dudn->Bot[i] = (3.0 * Upsi[i][0] - 4.0 * Upsi[i][1] + Upsi[i][2]) / (2.0 * dz);
			dudn->In[i] = (3.0 * Upsi[0][i] - 4.0 * Upsi[1][i] + Upsi[2][i]) / (2.0 * dx);
			dudn->Out[i] = (3.0 * Upsi[nmax][i] - 4.0 * Upsi[nmax - 1][i] + Upsi[nmax - 2][i]) / (2.0 * dx);
		}
	}

	free_PsiGrid(USol);

	/*
	 * Finally, loop around the boundary in order to find each boundary point...
	 * using the Lackner/Von Hagenow Formulation...
	 */

	/*Bottom*/
	for (ix = 0; ix <= nmax; ix++)
		pg->Psi[ix][0] += LHLoop(lha, pg, dudn, ix, 0);

	/*Top*/
	for (ix = 0; ix <= nmax; ix++) {
		if (pg->Symmetric)
			pg->Psi[ix][nmax] = pg->Psi[ix][0];
		else
			pg->Psi[ix][nmax] += LHLoop(lha, pg, dudn, ix, nmax);
	}

	/*Inside*/
	if (pg->Symmetric)
		for (iz = 1; iz <= nmax / 2; iz++) {
			pg->Psi[0][iz] += LHLoop(lha, pg, dudn, 0, iz);
			pg->Psi[0][nmax - iz] = pg->Psi[0][iz];
	} else {
		for (iz = 1; iz < nmax; iz++)
			pg->Psi[0][iz] += LHLoop(lha, pg, dudn, 0, iz);
	}

	/*Outside*/
	if (pg->Symmetric)
		for (iz = 1; iz <= nmax / 2; iz++) {
			pg->Psi[nmax][iz] += LHLoop(lha, pg, dudn, nmax, iz);
			pg->Psi[nmax][nmax - iz] = pg->Psi[nmax][iz];
	} else {
		for (iz = 1; iz < nmax; iz++)
			pg->Psi[nmax][iz] += LHLoop(lha, pg, dudn, nmax, iz);
	}

	free_LHvec(dudn, nmax);
}

/*
**
**
** Finds psi along outer boundary due to a coil set.
**
**
*/
void          PsiCoilBoundary(PSIGRID * pg, COIL * c)
{
	int           nmax, ix, iz;
	double        cur;
	COILGREEN    *cg;

	if (c->Enabled==0) return;

	nmax = pg->Nsize;

	cg = c->CoilGreen;
	cur = c->CoilCurrent;

	/*Top and Bottom*/
	for (ix = 0; ix <= nmax; ix++) {
		pg->Psi[ix][nmax] += cg->Top[ix] * cur;
		pg->Psi[ix][0] += cg->Bot[ix] * cur;
	}
	/*In and Out*/
	for (iz = 1; iz < nmax; iz++) {
		pg->Psi[0][iz] += cg->In[iz] * cur;
		pg->Psi[nmax][iz] += cg->Out[iz] * cur;
	}
}

/*
**
**
** Finds psi along outer boundary due to a conducting shell.
**
**
*/
void          PsiShellBoundary(PSIGRID * pg, SHELL * c)
{
	int           iss, nmax, ix, iz;
	double        cur;
	SUBSHELL     *subshell;
	SHELLGREEN   *sg;

	nmax = pg->Nsize;

	for (iss = 0; iss < c->NumSubShells; iss++) {
		subshell = c->SubShells[iss];
		sg = subshell->ShellGreen;
		cur = subshell->Current;
		/*Top and Bottom*/
		for (ix = 0; ix <= nmax; ix++) {
			pg->Psi[ix][nmax] += sg->Top[ix] * cur;
			pg->Psi[ix][0] += sg->Bot[ix] * cur;
		}
		/*In and Out*/
		for (iz = 1; iz < nmax; iz++) {
			pg->Psi[0][iz] += sg->In[iz] * cur;
			pg->Psi[nmax][iz] += sg->Out[iz] * cur;
		}
	}
}

/*
**
**	PsiBoundary
**
*/
void          PsiBoundary(TOKAMAK * td)
{
	int           i, nmax;
	double        MaxPsiBndry, DelPsiBndry;
	double        op, np;
	double      **Psi;
	PSIGRID      *pg;
	LHVEC        *oldBndry;

	pg = td->PsiGrid;
	nmax = pg->Nsize;
	Psi = pg->Psi;

	printf("INFO:	Computing Psi boundary.\n");
	fprintf(LogFile, "INFO:	Computing Psi boundary.\n");

	/* S A V E   O L D   B O U N D A R Y */
	oldBndry = new_LHvec(nmax);
	for (i = 0; i <= nmax; i++) {
		oldBndry->Top[i] = Psi[i][nmax];
		oldBndry->Bot[i] = Psi[i][0];
		oldBndry->In[i] = Psi[0][i];
		oldBndry->Out[i] = Psi[nmax][i];
	}

    MULTI;


	/* Z E R O   P S I   B O U N D A R Y */
	for (i = 0; i <= nmax; i++)
		Psi[i][nmax] = Psi[i][0] = Psi[0][i] = Psi[nmax][i] = 0.0;

    MULTI;

	/* G E T    P L A S M A    B O U N D A R Y */
	PsiPlasmaBoundary(td->LHPlasmaGreen, pg);

    MULTI;

	/* A D D    C O I L    B O U N D A R Y */
	for (i = 0; i < td->NumCoils; i++)
		PsiCoilBoundary(pg, td->Coils[i]);

    MULTI;

	/* A D D   S H E L L    B O U N D A R Y */
	for (i = 0; i < td->NumShells; i++)
		PsiShellBoundary(pg, td->Shells[i]);

    MULTI;

	/* C O M P A R E   O L D   W I T H   N E W */
	MaxPsiBndry = fabs(oldBndry->Bot[0] + Psi[0][0]);
	DelPsiBndry = fabs(oldBndry->Bot[0] - Psi[0][0]);
	for (i = 0; i <= nmax; i++) {
		op = oldBndry->Top[i];
		np = Psi[i][nmax];
		if (fabs(op - np) > DelPsiBndry)
			DelPsiBndry = fabs(op - np);
		if (fabs(op + np) > MaxPsiBndry)
			MaxPsiBndry = fabs(op + np);
		op = oldBndry->Bot[i];
		np = Psi[i][0];
		if (fabs(op - np) > DelPsiBndry)
			DelPsiBndry = fabs(op - np);
		if (fabs(op + np) > MaxPsiBndry)
			MaxPsiBndry = fabs(op + np);
		op = oldBndry->In[i];
		np = Psi[0][i];
		if (fabs(op - np) > DelPsiBndry)
			DelPsiBndry = fabs(op - np);
		if (fabs(op + np) > MaxPsiBndry)
			MaxPsiBndry = fabs(op + np);
		op = oldBndry->Out[i];
		np = Psi[nmax][i];
		if (fabs(op - np) > DelPsiBndry)
			DelPsiBndry = fabs(op - np);
		if (fabs(op + np) > MaxPsiBndry)
			MaxPsiBndry = fabs(op + np);
	}

	pg->BoundError = DelPsiBndry / MaxPsiBndry;

	printf("		[BoundErr = %g]\n", pg->BoundError);
	fprintf(LogFile, "		[BoundErr = %g]\n", pg->BoundError);

	free_LHvec(oldBndry, nmax);
}
