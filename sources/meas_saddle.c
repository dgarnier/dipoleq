/*
** TokaMac v2.0
**
** meas_saddle.c
**
** Source file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_saddle_Green(TOKAMAK *td, MEAS *m)
**
**				meas_saddle_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_saddle_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_saddle.c
** Date:		March 24, 1993
**
** Revisions:
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
#include "meas_saddle.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307
#define DEG_RAD		1.74532925e-2

extern FILE  *LogFile;

/*
**
**	g_saddle
**
**	The basic Greens function for a meas_saddle.
**
*/
double        g_saddle(double x1, double z1, double x2, double z2, double xc, double zc)
{

	return Green(x1, z1, xc, zc) - Green(x2, z2, xc, zc);
}

/*
**	meas_saddle_Green
**
**
*/
void          meas_saddle_Green(TOKAMAK * td, MEAS * m)
{
	PSIGRID      *pg;
	COIL         *c;
	SUBCOIL      *sc;
	SHELL        *s;
	SUBSHELL     *ss;
	int           ncoil, nsubshell = 0, nmax, i, isc, iss, ix, iz;
	double        x1, x2, z1, z2, sum;
	double       *CGrn, *SGrn, **PlGrn;

	pg = td->PsiGrid;

	nmax = pg->Nsize;
	ncoil = td->NumCoils;

	for (i = 0; i < td->NumShells; i++)
		nsubshell += td->Shells[i]->NumSubShells;

	x1 = m->parm.saddle.X1;
	z1 = m->parm.saddle.Z1;
	x2 = m->parm.saddle.X2;
	z2 = m->parm.saddle.Z2;

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
			sum += sc->CurrentFraction * g_saddle(x1, z1, x2, z2, sc->X, sc->Z);
			if (pg->Symmetric)
				sum += sc->CurrentFraction * g_saddle(x1, z1, x2, z2, sc->X, -(sc->Z));
		}
		CGrn[i] = sum;
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
			sum += g_saddle(x1, z1, x2, z2, ss->X, ss->Z);
			if (pg->Symmetric)
				sum += g_saddle(x1, z1, x2, z2, ss->X, -(ss->Z));
			SGrn[isc] = sum;
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
			PlGrn[ix][iz] = g_saddle(x1, z1, x2, z2, pg->X[ix], pg->Z[iz]);

	printf("		[%s]\n", m->Name);
	fprintf(LogFile, "		[%s]\n", m->Name);
}

/*
**	meas_saddle_Fit
**
**
*/
void          meas_saddle_Fit(TOKAMAK * td, MEAS * m)
{
	meas_mag_Fit(td, m);
}

/*
**	meas_saddle_Now
**
**
*/
void          meas_saddle_Now(TOKAMAK * td, MEAS * m)
{
	meas_mag_Now(td, m);
}

/*
**	meas_saddle_L
**
**
*/
void          meas_saddle_L(TOKAMAK * td, MEAS * m, double *L)
{
	meas_mag_L(td, m, L);
}
