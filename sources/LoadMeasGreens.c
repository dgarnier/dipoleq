/*
** TokaMac v2.0
**
** LoadMeasGreens.c
**
** Reads the Meas Green's functions from disk.
** Also, it rewrites the file (erasing if it already exists),
** by calling RewriteMeasGreens.
**
** File:		LoadMeasGreens.c
** Date:		March 24, 1993
**
** Revisions:
**
**		August 6, 1993		Added perfectly conducting shells
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdlib.h>
#include "VAX_Alloc.h"
#include <stdio.h>
#include "nrutil.h"
#include "stdio_dmatrix.h"
#include "coil.h"
#include "shell.h"
#include "green.h"
#include "tokamak.h"
#include "measurement.h"
#include "LoadMeasGreens.h"

#define cchk(c,n)	if ((c) != ((size_t) n)) nrerror("ERROR: Could not read/write Measurement Greens.")

extern FILE  *LogFile;

/*
**
**	LoadMeasGreens
**
*/
void          LoadMeasGreens(TOKAMAK * td)
{
	FILE         *fi;
	MEAS         *m;
	size_t        c;
	int           i, nmax, ncoil, nsubshells = 0, ngrid, im;

	nmax = td->PsiGrid->Nsize;
	ncoil = td->NumCoils;
	ngrid = (nmax + 1) * (nmax + 1);

	for (i = 0; i < td->NumShells; i++)
		nsubshells += td->Shells[i]->NumSubShells;

	if (td->MGreenStatus == MGreenOK) {	/* read MeasGreens from file */
		fi = fopen(td->MGname, "rb");
		if (!fi)
			nrerror("ERROR:	Could not open MeasGreen file for reading.");

		/* R E A D    M E A S   G R E E N S */
		for (im = 0; im < td->NumMeasures; im++) {
			m = td->Measures[im];
			if (m->FindGreen) {
				m->CoilGreen = dvector(0, ncoil - 1);
				c = fread(m->CoilGreen, sizeof(double), ncoil, fi);
				cchk(c, ncoil);
				m->PlasmaGreen = dmatrix(0, nmax, 0, nmax);
				c = fread_dmatrix(m->PlasmaGreen, 0, nmax, 0, nmax, fi);
				cchk(c, ngrid);
				if (nsubshells > 0) {
					m->ShellGreen = dvector(0, nsubshells - 1);
					c = fread(m->ShellGreen, sizeof(double), nsubshells, fi);
					cchk(c, nsubshells);
				}
			}
		}

		fclose(fi);
	} else
		RewriteMeasGreens(td);
}

/*
**
**	RewriteMeasGreens
**
*/
void          RewriteMeasGreens(TOKAMAK * td)
{
	FILE         *fi;
	MEAS         *m;
	size_t        c;
	int           i, nmax, ncoil, nsubshells = 0, ngrid, im;

	nmax = td->PsiGrid->Nsize;
	ncoil = td->NumCoils;
	ngrid = (nmax + 1) * (nmax + 1);

	for (i = 0; i < td->NumShells; i++)
		nsubshells += td->Shells[i]->NumSubShells;

	printf("INFO:	Finding Measure Greens...\n");
	fprintf(LogFile, "INFO:	Finding Measure Greens...\n");
	for (im = 0; im < td->NumMeasures; im++) {
		m = td->Measures[im];
		if (m->FindGreen)
			(*(m->FindGreen)) (td, m);
	}

	/* don't write if not a good filename */
	if ((td->MGname[0] == '\0') || (td->MGname[0] == '*'))
		return;

	fi = fopen(td->MGname, "wb");
	if (!fi)
		nrerror("ERROR:	Could not open file for writing in MeasGreen.");

	/* W R I T E    M E A S   G R E E N S */
	for (im = 0; im < td->NumMeasures; im++) {
		m = td->Measures[im];
		if (m->FindGreen) {
			c = fwrite(m->CoilGreen, sizeof(double), ncoil, fi);
			cchk(c, ncoil);
			c = fwrite_dmatrix(m->PlasmaGreen, 0, nmax, 0, nmax, fi);
			cchk(c, ngrid);
			c = fwrite(m->ShellGreen, sizeof(double), nsubshells, fi);
			cchk(c, nsubshells);
		}
	}

	fclose(fi);

	td->MGreenStatus = MGreenOK;
}

/*
**
**	free_MeasGreens
**
*/
void          free_MeasGreens(TOKAMAK * td)
{
	int           i, im, nmax, ncoil, nsubshells = 0;
	MEAS         *m;

	nmax = td->PsiGrid->Nsize;
	ncoil = td->NumCoils;

	for (i = 0; i < td->NumShells; i++)
		nsubshells += td->Shells[i]->NumSubShells;

	for (im = 0; im < td->NumMeasures; im++) {
		m = td->Measures[im];
		if (m->CoilGreen)
			free_dvector(m->CoilGreen, 0, ncoil - 1);
		if (m->ShellGreen)
			free_dvector(m->ShellGreen, 0, nsubshells - 1);
		if (m->PlasmaGreen)
			free_dmatrix(m->PlasmaGreen, 0, nmax, 0, nmax);
		m->CoilGreen = NULL;
		m->ShellGreen = NULL;
		m->PlasmaGreen = NULL;
	}
}
