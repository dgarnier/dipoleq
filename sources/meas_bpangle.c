/*
** TokaMac v2.0
**
** meas_bpangle.c
**
** Source file for local poloidal field measurements as expressed by a
** LINEAR poloidal field angle.  That is the deviation of the toroidal
** field from its vacuum value is IGNORED.  Thus, bpangle = atan(Bz,Bt).
**
** Note, the bpangle depends on both the POLOIDAL field and the TOROIDAL
** field.  The POLOIDAL field part depends upon the toroidal currents.
** The TOROIDAL field part depends only upon those terms contributing to
** Bt.  This is G2p for IsoNoFlow.
**
**
** Every measurement must define the following subroutines:
**
**				meas_bpangle_Green(TOKAMAK *td, MEAS *m)
**
**				meas_bpangle_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_bpangle_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_bpangle.c
** Date:		November 15, 1993
**
** Revisions:
**
**  DTG 12/9/98 --- Found bug??? g_bp in place of g_pbangle!
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
#include "meas_bpangle.h"
#include "fpoly.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307
#define DEG_RAD		1.74532925e-2

#define CHECK_PLASMA	(GetIsPlasma(pg, m->X, m->Z) > 0.5)

extern FILE  	*LogFile;
extern double 	***meas_dJdy;		/* common pointer to storage */

double        g_bpangle(double x, double z, double xc, double zc);


/*
**
**	g_bpangle
**
**	The basic Greens function for a meas_bpangle.
**  Returns a value proportional to the vertical poloidal field.
**
*/
double        g_bpangle(double x, double z, double xc, double zc)
{
	double        GG, dGx, dGz;

	GetdGreen(&GG, &dGx, &dGz, x, z, xc, zc);
	return (-dGx);
}

/*
**	meas_bpangle_Green
**
**  Note:  This simply stores the Green's functions needed to find the
** 			vertical component of the poloidal field.
*/
void          meas_bpangle_Green(TOKAMAK * td, MEAS * m)
{
	PSIGRID      *pg;
	COIL         *c;
	SUBCOIL      *sc;
	SHELL        *s;
	SUBSHELL     *ss;
	int           ncoil, nsubshell = 0, nmax, i, isc, iss, ix, iz;
	double        x, z, sum;
	double       *CGrn, *SGrn, **PlGrn;

	pg = td->PsiGrid;

	nmax = pg->Nsize;
	ncoil = td->NumCoils;

	for (i = 0; i < td->NumShells; i++)
		nsubshell += td->Shells[i]->NumSubShells;

	x = m->X;
	z = m->Z;

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
			sum += sc->CurrentFraction * g_bpangle(x, z, sc->X, sc->Z);
			if (pg->Symmetric)
				sum += sc->CurrentFraction * g_bpangle(x, z, sc->X, -(sc->Z));
		}
		CGrn[i] = sum / TWOPI / x;
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
			/* DTG this looks like a bug! is it? */
			// sum += g_bp(x, z, ss->X, ss->Z);
			sum += g_bpangle(x, z, ss->X, ss->Z);
			if (pg->Symmetric)
				sum += g_bpangle(x, z, ss->X, -(ss->Z));
			SGrn[isc] = sum / TWOPI / x;
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
			PlGrn[ix][iz] = g_bpangle(x, z, pg->X[ix], pg->Z[iz]) / TWOPI / x;

	printf("		[%s]\n", m->Name);
	fprintf(LogFile, "		[%s]\n", m->Name);
}

/*
**	meas_bpangle_Fit
**
**
*/
void          meas_bpangle_Fit(TOKAMAK * td, MEAS * m)
{
	double        Psi, PsiX;
	double        G2, DelPsi, Bt;
	int         **ip;
	PLASMA       *pl;
	PSIGRID      *pg;

	pl = td->Plasma;
	pg = td->PsiGrid;
	ip = pg->IsPlasma;

	meas_mag_Fit(td, m);	/* This computes the poloidal field */

	/* Compute the toroidal magnetic field */

	switch (pl->ModelType) {
	  case Plasma_Std:
		  Psi = GetPsi(pg, m->X, m->Z);
		  DelPsi = pg->PsiLim - pg->PsiAxis;
		  PsiX = (Psi - pg->PsiAxis) / DelPsi;
		  if ((CHECK_PLASMA) && (PsiX < 1.0))
			  G2 = 1.0 - (DelPsi / pl->StndG) * pow(1.0 - PsiX, pl->StndG) * pl->G2p[1];
		  else
			  G2 = 1.0;
		  break;
	  case Plasma_IsoNoFlow:
		  Psi = GetPsi(pg, m->X, m->Z);
		  DelPsi = pg->PsiLim - pg->PsiAxis;
		  PsiX = (Psi - pg->PsiAxis) / DelPsi;
		  if ((CHECK_PLASMA) && (PsiX < 1.0))
			  G2 = fpoly_int(pl->G2p, PsiX, pl->G2pTerms, DelPsi, 1.0);
		  else
			  G2 = 1.0;
		  break;
	  case Plasma_IsoFlow:
		  break;
	  case Plasma_AnisoNoFlow:
		  break;
	  case Plasma_AnisoFlow:
		  break;
	}

	Bt = pl->B0R0*sqrt(G2)/m->X;
	m->Fit = atan2(m->Fit,Bt);		/* atan(Bz/Bt) */
}

/*
**	meas_bpangle_Now
**
**
*/
void          meas_bpangle_Now(TOKAMAK * td, MEAS * m)
{
	double        Psi, PsiX;
	double        G2, DelPsi, Bt;
	int         **ip;
	PLASMA       *pl;
	PSIGRID      *pg;

	pl = td->Plasma;
	pg = td->PsiGrid;
	ip = pg->IsPlasma;

	meas_mag_Now(td, m);	/* This computes the poloidal field */

	/* Compute the toroidal magnetic field */

	switch (pl->ModelType) {
	  case Plasma_Std:
		  Psi = GetPsi(pg, m->X, m->Z);
		  DelPsi = pg->PsiLim - pg->PsiAxis;
		  PsiX = (Psi - pg->PsiAxis) / DelPsi;
		  if ((CHECK_PLASMA) && (PsiX < 1.0))
			  G2 = 1.0 - (DelPsi / pl->StndG) * pow(1.0 - PsiX, pl->StndG) * pl->G2p[1];
		  else
			  G2 = 1.0;
		  break;
	  case Plasma_IsoNoFlow:
		  Psi = GetPsi(pg, m->X, m->Z);
		  DelPsi = pg->PsiLim - pg->PsiAxis;
		  PsiX = (Psi - pg->PsiAxis) / DelPsi;
		  if ((CHECK_PLASMA) && (PsiX < 1.0))
			  G2 = fpoly_int(pl->G2p, PsiX, pl->G2pTerms, DelPsi, 1.0);
		  else
			  G2 = 1.0;
		  break;
	  case Plasma_IsoFlow:
		  break;
	  case Plasma_AnisoNoFlow:
		  break;
	  case Plasma_AnisoFlow:
		  break;
	}

	Bt = pl->B0R0*sqrt(G2)/m->X;
	m->Now = atan2(m->Now,Bt);		/* atan(Bz/Bt) */
}

/*
**	meas_bpangle_L
**
**  Note:  We assume in this routine that meas->Now has already been
**			computed correctly.  We use this value in computing L.
**
**	Note:  The Green's function only specify the vertical component of
**			the poloidal field.
**
*/
void          meas_bpangle_L(TOKAMAK * td, MEAS * m, double *L)
{
	int           iu, iup, gt;
	int           ix, iz, nmax;
	double      **PlGrn, *CGrn;
	double        Psi, PsiX;
	double        G2, DelPsi, Bt, Bz, c4;
	int         **ip;
	PLASMA       *pl;
	PSIGRID      *pg;

	pl = td->Plasma;
	pg = td->PsiGrid;
	ip = pg->IsPlasma;

	gt = pl->G2pTerms;

	for (iu = 1; iu <= td->NumUnkns; iu++)
		L[iu] = 0.0;

	/* Compute the toroidal magnetic field */

	switch (pl->ModelType) {
	  case Plasma_Std:
		  Psi = GetPsi(pg, m->X, m->Z);
		  DelPsi = pg->PsiLim - pg->PsiAxis;
		  PsiX = (Psi - pg->PsiAxis) / DelPsi;
		  if ((CHECK_PLASMA) && (PsiX < 1.0))
			  G2 = 1.0 - (DelPsi / pl->StndG) * pow(1.0 - PsiX, pl->StndG) * pl->G2p[1];
		  else
			  G2 = 1.0;
		  break;
	  case Plasma_IsoNoFlow:
		  Psi = GetPsi(pg, m->X, m->Z);
		  DelPsi = pg->PsiLim - pg->PsiAxis;
		  PsiX = (Psi - pg->PsiAxis) / DelPsi;
		  if ((CHECK_PLASMA) && (PsiX < 1.0))
			  G2 = fpoly_int(pl->G2p, PsiX, gt, DelPsi, 1.0);
		  else
			  G2 = 1.0;
		  break;
	  case Plasma_IsoFlow:
		  break;
	  case Plasma_AnisoNoFlow:
		  break;
	  case Plasma_AnisoFlow:
		  break;
	}

	/* Use previously calculated m->Now, in order to find present Bz */
	Bt = pl->B0R0*sqrt(G2)/m->X;
	Bz = tan(m->Now)*Bt;		/* atan(Bz/Bt) */

	/* Compute dependencies on toroidal currents */
	nmax = td->PsiGrid->Nsize;
	iup = td->NumUnkns - td->NumCoils;	/* iup is the number of plasma unknowns */

	PlGrn = m->PlasmaGreen;
	CGrn = m->CoilGreen;

	/* P L A S M A   U N K N O W N S */
	for (iu = 1; iu <= iup; iu++) {
		for (ix = 1; ix < nmax; ix++)
			for (iz = 1; iz < nmax; iz++)
				L[iu] += PlGrn[ix][iz] * meas_dJdy[iu][ix][iz];
		L[iu] *= td->PsiGrid->dx * td->PsiGrid->dz;
	}

	/* C O I L   U N K N O W N S */
	for (iu = iup + 1; iu <= td->NumUnkns; iu++)
		L[iu] = CGrn[iu - iup - 1];

	/* A D D   T O R O I D A L   F I E L D   D E P E N D E N C I E S */
	switch (pl->ModelType) {
	  case Plasma_Std:
		  /* --------   L[1] = G2p[1]  ------- */
		  if ((CHECK_PLASMA) && (PsiX < 1.0))
			  L[1] += 0.5*(Bz/Bt/Bt)*(DelPsi / pl->StndG) * pow(1.0 - PsiX, pl->StndG);
		  break;
	  case Plasma_IsoNoFlow:
		  if ((CHECK_PLASMA) && (PsiX < 1.0)) {
			  /* --------   G2p[i] ------- */
			  c4 = pow(PsiX, gt);
			  for (iu = 1; iu <= gt; iu++)
				  L[iu] += 0.5*(Bz/Bt/Bt)*DelPsi * ((1.0 - pow(PsiX, iu)) / iu -
								   (1.0 - pow(PsiX, gt + 1)) / (gt + 1));
		  }
		  break;
	  case Plasma_IsoFlow:
		  break;
	  case Plasma_AnisoNoFlow:
		  break;
	  case Plasma_AnisoFlow:
		  break;
	}

	/* A D J U S T   F O R   A R C T A N G E N T */
	for (iu = 1; iu <= td->NumUnkns; iu++)
		L[iu] *= (Bt/(Bt*Bt + Bz*Bz));
}

#undef CHECK_PLASMA
