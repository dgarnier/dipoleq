/*
** TokaMac v2.0
**
** meas_coilcur.c
**
** Source file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_coilcur_Green(TOKAMAK *td, MEAS *m)
**
**				meas_coilcur_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_coilcur_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_coilcur.c
** Date:		March 24, 1993
**
** Revisions:
**
**		August 6, 1993		Added xxx_Now
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdio.h>
#include <math.h>
#include "measurement.h"
#include "coil.h"
#include "tokamak.h"
#include "meas_coilcur.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

/*
** Local variables to this file.
*/

/*
**	meas_coilcur_Fit
**
**
*/
void          meas_coilcur_Fit(TOKAMAK * td, MEAS * m)
{
	int           cnum;
	COIL         *c;

	cnum = m->parm.coilcur.CoilNum;
	c = td->Coils[cnum];

	m->Fit = c->CoilCurrent / MU0;
}

/*
**	meas_coilcur_Now
**
**
*/
void          meas_coilcur_Now(TOKAMAK * td, MEAS * m)
{
	int           cnum;
	COIL         *c;

	cnum = m->parm.coilcur.CoilNum;
	c = td->Coils[cnum];

	m->Now = c->CoilCurrent / MU0;
}

/*
**	meas_coilcur_L
**
**	Note:	The L array has dimensions [1..NumUnkns].
**
*/
void          meas_coilcur_L(TOKAMAK * td, MEAS * m, double *L)
{
	int           iu, dnum;		/* the "diagonal" term of the L vector */
	int           cnum;

	cnum = m->parm.coilcur.CoilNum;	/* td->Coils[i] is a zero-referenced array */

	for (iu = 1; iu <= td->NumUnkns; iu++)
		L[iu] = 0.0;

	dnum = td->NumUnkns - td->NumCoils + 1 + cnum;

	L[dnum] = 1.0 / MU0;
}
