/*
** TokaMac v2.0
**
** meas_circle.c
**
** Source file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_circle_Green(TOKAMAK *td, MEAS *m)
**
**				meas_circle_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_circle_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_circle.c
** Date:		March 24, 1993
**
** Revised:
**
**		August 6, 1993		Added perfectly conducting shells
**		August 6, 1993		Added xxx_Now
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdio.h>
#include <math.h>
#include "green.h"
#include "nrutil.h"
#include "measurement.h"
#include "coil.h"
#include "shell.h"
#include "psigrid.h"
#include "plasma.h"
#include "tokamak.h"
#include "meas_mag_Xxx.h"
#include "meas_circle.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307
#define DEG_RAD		1.74532925e-2

#define TRUE		1
#define FALSE		0

extern FILE  *LogFile;

/*
**
**	g_circle
**
**	The basic Greens function for a meas_circle.
**
*/
double        g_circle(MEAS * m, int symmetric, double xc, double zc)
{
	double        GG, dGx, dGz;
	double        x, z, rad, sum = 0.0;
	double        thetaLoc, cosLoc, sinLoc, xLoc, zLoc, ang, trig;
	int           num, iloc;

	num = m->parm.circle.Number;
	x = m->X;
	z = m->Z;
	rad = m->parm.circle.Radius;

	for (iloc = 1; iloc <= num; iloc++) {
		thetaLoc = TWOPI * ((iloc - 1.0) / num);
		cosLoc = cos(thetaLoc);
		sinLoc = sin(thetaLoc);
		xLoc = x + rad * cosLoc;
		zLoc = z + rad * sinLoc;
		switch (m->parm.circle.CircleType) {
		  case CircleType_brsin:
		  case CircleType_brcos:
			  ang = (PI / 2.0) - thetaLoc;
			  break;
		  case CircleType_btcos:
			  ang = -thetaLoc;
			  break;
		}
		switch (m->parm.circle.CircleType) {
		  case CircleType_brcos:
		  case CircleType_btcos:
			  trig = cosLoc;
			  break;
		  case CircleType_brsin:
			  trig = sinLoc;
			  break;
		}

		GetdGreen(&GG, &dGx, &dGz, xLoc, zLoc, xc, zc);
		sum += (trig / xLoc) * (-dGx * cos(ang) + dGz * sin(ang));
		if (symmetric) {
			GetdGreen(&GG, &dGx, &dGz, xLoc, zLoc, xc, -zc);
			sum += (trig / xLoc) * (-dGx * cos(ang) + dGz * sin(ang));
		}
	}

	return sum / num;
}

/*
**	meas_circle_Green
**
**
*/
void          meas_circle_Green(TOKAMAK * td, MEAS * m)
{
	PSIGRID      *pg;
	COIL         *c;
	SUBCOIL      *sc;
	SHELL        *s;
	SUBSHELL     *ss;
	int           ncoil, nsubshell = 0, nmax, i, isc, iss, ix, iz;
	double       sum;
	double       *CGrn, *SGrn, **PlGrn;

	pg = td->PsiGrid;

	nmax = pg->Nsize;
	ncoil = td->NumCoils;

	for (i = 0; i < td->NumShells; i++)
		nsubshell += td->Shells[i]->NumSubShells;

	m->CoilGreen = CGrn = dvector(0, ncoil - 1);
	if (nsubshell > 0)
		m->ShellGreen = SGrn = dvector(0, nsubshell - 1);
	m->PlasmaGreen = PlGrn = dmatrix(0, nmax, 0, nmax);

	/* C O I L    G R E E N S */
	for (i = 0; i < ncoil; i++) {
		c = td->Coils[i];
		if (!c->Enabled)
			continue;
		sum = 0.0;
		for (isc = 0; isc < c->NumSubCoils; isc++) {
			sc = c->SubCoils[isc];
			sum += sc->CurrentFraction * g_circle(m, pg->Symmetric, sc->X, sc->Z);
		}
		CGrn[i] = sum / TWOPI;
	}

	/* S H E L L    G R E E N S */
	isc = 0;
	for (i = 0; i < td->NumShells; i++) {
		s = td->Shells[i];
		if (!s->Enabled)
			continue;
		for (iss = 0; iss < s->NumSubShells; iss++) {
			sum = 0.0;
			ss = s->SubShells[iss];
			sum += g_circle(m, pg->Symmetric, ss->X, ss->Z);
			SGrn[isc] = sum / TWOPI;
			isc++;
		}
	}

	/* P L A S M A    G R E E N S */
	for (ix = 0; ix <= nmax; ix++)
		PlGrn[ix][0] = PlGrn[ix][nmax] = 0.0;
	for (iz = 0; iz <= nmax; iz++)
		PlGrn[0][iz] = PlGrn[nmax][iz] = 0.0;

	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			PlGrn[ix][iz] = g_circle(m, FALSE, pg->X[ix], pg->Z[iz]) / TWOPI;

	printf("		[%s]\n", m->Name);
	fprintf(LogFile, "		[%s]\n", m->Name);
}

/*
**	meas_circle_Fit
**
**
*/
void          meas_circle_Fit(TOKAMAK * td, MEAS * m)
{
	meas_mag_Fit(td, m);
}

/*
**	meas_circle_Now
**
**
*/
void          meas_circle_Now(TOKAMAK * td, MEAS * m)
{
	meas_mag_Now(td, m);
}

/*
**	meas_circle_L
**
**
*/
void          meas_circle_L(TOKAMAK * td, MEAS * m, double *L)
{
	meas_mag_L(td, m, L);

}
