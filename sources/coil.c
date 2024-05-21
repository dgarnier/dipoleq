
/*
** TokaMac v2.0
**
** coil.c
**
** This file define basic global creation and destruction
** subroutines.
**
** File:		coil.c
** Date:		February 7, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdlib.h>
#include "VAX_Alloc.h"
#include <string.h>
#include <math.h>
#include "nrutil.h"
#include "coil.h"
#include "psigrid.h"

SUBCOIL      *copy_SubCoil(SUBCOIL * sc);

/*
**
**	new_SubCoil
**
*/

SUBCOIL      *new_SubCoil(void)
{
	SUBCOIL      *sc;

	sc = (SUBCOIL *) malloc((size_t) sizeof(SUBCOIL));
	if (sc == NULL)
		nrerror("ERROR: Allocation error in new_SubCoil");

	sc->X = 1.0;
	sc->Z = 0.0;
	sc->CurrentFraction = 1.0;
	strcpy(sc->Name, " - ");

	return sc;
}

SUBCOIL      *copy_SubCoil(SUBCOIL * sc)
{
	SUBCOIL      *new_sc;

	new_sc = new_SubCoil();
	new_sc->X = sc->X;
	new_sc->Z = sc->Z;
	new_sc->CurrentFraction = sc->CurrentFraction;
	strcpy(new_sc->Name, sc->Name);

	return new_sc;
}

/*
**
**	add_SubCoil
**
*/

void          add_SubCoil(COIL * c, SUBCOIL * sc)
{
	SUBCOIL     **newAry;
	int           is, scnum;

	scnum = c->NumSubCoils;

	/*   realloc(..) appears to fail with VAXC !!
	newAry = (SUBCOIL **) realloc(c->SubCoils, (size_t) (scnum + 1) * sizeof(SUBCOIL *));
	if (newAry == NULL)
		nrerror("allocation error in add_SubCoil");
	*/

	newAry = (SUBCOIL **) malloc((size_t) (scnum + 1) * sizeof(SUBCOIL *));
	if (newAry == NULL)
		nrerror("allocation error in add_SubCoil");

	for (is = 0; is < scnum; is++) {
		newAry[is] = copy_SubCoil(c->SubCoils[is]);
		free(c->SubCoils[is]);
	}
	free(c->SubCoils);

	c->SubCoils = newAry;
	c->SubCoils[scnum] = sc;
	c->NumSubCoils = scnum + 1;
}

/*
**
** 	new_Coil
**
*/

COIL         *new_Coil(int NumSubCoils)
{
	COIL         *c;
	int           isc;

	c = (COIL *) malloc((size_t) sizeof(COIL));
	if (c == NULL)
		nrerror("ERROR: Allocation error1 in new_Coil.");

	c->NumSubCoils = NumSubCoils;
	c->CoilCurrent = 0.0;
	c->Enabled = Coil_Off;
	c->CoilGreen = NULL;
	strcpy(c->Name, " - ");
	c->X = c->Z = 0.0;
	c->dX = c->dZ = -1.0;

	if (NumSubCoils == 0)
		c->SubCoils = NULL;
	else {
		c->SubCoils = (SUBCOIL **) malloc((size_t) NumSubCoils * sizeof(SUBCOIL *));
		if (c->SubCoils == NULL)
			nrerror("ERROR: Allocation error2 in new_Coil.");
		for (isc = 0; isc < NumSubCoils; isc++)
			c->SubCoils[isc] = new_SubCoil();
	}
	return c;
}

void          free_Coil(COIL * c, int nmax)
{
	int           isc;

	if (c->CoilGreen)
		free_CoilGreen(c->CoilGreen, nmax);
	for (isc = 0; isc < c->NumSubCoils; isc++)
		if (c->SubCoils[isc])
			free(c->SubCoils[isc]);
	free(c->SubCoils);
	free(c);
}

COILGREEN    *new_CoilGreen(int nmax)
{
	COILGREEN    *cg;
	double       *gvec;

	/* alloc space for new lhvec */
	cg = (COILGREEN *) malloc((size_t) sizeof(COILGREEN));
	if (!cg)
		nrerror("ERROR: Allocation failure in new_CoilGreen.");

	gvec = dvector(0, nmax);
	cg->Top = gvec;

	gvec = dvector(0, nmax);
	cg->Bot = gvec;

	gvec = dvector(0, nmax);
	cg->In = gvec;

	gvec = dvector(0, nmax);
	cg->Out = gvec;

	return cg;
}

void          free_CoilGreen(COILGREEN * cg, int nmax)
{
	free_dvector(cg->Top, 0, nmax);
	free_dvector(cg->Bot, 0, nmax);
	free_dvector(cg->In, 0, nmax);
	free_dvector(cg->Out, 0, nmax);

	free(cg);
}

/*   routine to calculate a new subcoil set based on the coil information.
**   added because I hated adding these in the input file and I wanted the
**   the best possible use of subcoils on the grid, which is variable.
**
**   if I was really clever I could make a class of coil and get rid of this
**   whole subcoil idea....  oh well.


     Added 10/19/99
*/


void compute_SubCoils(COIL *c, PSIGRID * pg)
{
	double cdx, cdz, cxm, czm;
	double gdx, gdz, gxm, gzm;
	int	   isc, ix, ix1, ix2, iz, iz1, iz2, nSub = 0;
	double fx, fx1, fx2, fxn, fz, fz1, fz2, fzn, jcur, jtot;
	SUBCOIL **sca;

	/* free old subcoil array */

	if ((c->NumSubCoils != 0) && (c->SubCoils != NULL)) {
		for (isc = 0 ; isc > c->NumSubCoils; isc++)
		if (c->SubCoils[isc])
			free(c->SubCoils[isc]);
		free(c->SubCoils);
	}

	gdx = pg->dx;
	gdz = pg->dz;
	gxm = pg->Xmin;
	gzm = pg->Zmin;

	cdx = c->dX;
	cdz = c->dZ;
	cxm = c->X-cdx/2.0;
	czm = c->Z-cdz/2.0;

	jcur = gdx*gdz/(cdx*cdz);

	ix1 = (int) floor((cxm-gxm)/gdx + 0.5);
	ix2 = (int) floor((cxm+cdx-gxm)/gdx + 0.5);
	iz1 = (int) floor((czm-gzm)/gdz + 0.5);
	iz2 = (int) floor((czm+cdz-gzm)/gdz + 0.5);

	nSub = (1+ix2-ix1)*(1+iz2-iz1);

	c->SubCoils = (SUBCOIL **) malloc((size_t) nSub * sizeof(SUBCOIL *));
	if (c->SubCoils == NULL)
		nrerror("ERROR: Allocation error2 in compute_SubCoils.");
	for (sca = c->SubCoils; sca < c->SubCoils+nSub; sca++)
		*sca = new_SubCoil();
	c->NumSubCoils = nSub;
	jtot=0;
	sca = c->SubCoils;

	fx1 = (cxm-gxm)/gdx - ix1;
	fx2 = (cxm+cdx-gxm)/gdx - ix2;
	fxn = fx1;

	fz1 = (czm-gzm)/gdz - iz1;
	fz2 = (czm+cdz-gzm)/gdz - iz2;
	fzn = fz1;

	for (ix=ix1; ix <= ix2; ix++) {
		if (ix < ix2) {
			fx = 0.5 - fxn;
			fxn = - 0.5;
		} else {
			fx = fx2 - fxn;
			fxn = fx1;
		}
		for (iz=iz1; iz <= iz2; iz++) {
			if (iz < iz2) {
				fz = 0.5 - fzn;
				fzn = - 0.5;
			} else {
				fz = fz2 - fzn;
				fzn = fz1;
			}
			(*sca)->X = gxm + ix*gdx;
			(*sca)->Z = gzm + iz*gdz;
			(*sca)->CurrentFraction = fx*fz*jcur;
			sca++;
			jtot+=fx*fz*jcur;   /* if you get this calc right, ftot should be 1 */
		}
	}
}
