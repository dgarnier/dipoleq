/*
** TokaMac v2.0
**
** LoadShellGreens.c
**
** Reads the Shell Green's functions from disk.
** Also, it rewrites the file (erasing if it already exists),
** by calling RewriteShellGreens.
**
** File:		LoadShellGreens.c
** Date:		August 4, 1993
**
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
#include "LoadShellGreens.h"

#define cchk(c,n)	if ((c) != ((size_t) n)) nrerror("ERROR: Could not read Shell Greens.")

extern FILE  *LogFile;

/*
**	Find_SubShellGreens
**
*/
void          Find_SubShellGreens(TOKAMAK * td, SUBSHELL * ss)
{
	int           ic, isc, ix, iz, nmax, ncoil;
	int           sym;
	double       *X, *Z;
	double      **pg;
	double       *cg;
	COIL         *coil;
	SUBCOIL      *subcoil;

	X = td->PsiGrid->X;
	Z = td->PsiGrid->Z;
	nmax = td->PsiGrid->Nsize;
	ncoil = td->NumCoils;
	sym = td->PsiGrid->Symmetric;

	pg = dmatrix(0, nmax, 0, nmax);
	cg = dvector(0, ncoil - 1);

	/* P L A S M A   G R E E N */
	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			pg[ix][iz] = -Green(ss->X, ss->Z, X[ix], Z[iz]);

	/* C O I L   S E T S  */
	for (ic = 0; ic < ncoil; ic++) {
		coil = td->Coils[ic];
		cg[ic] = 0.0;
		for (isc = 0; isc < coil->NumSubCoils; isc++) {
			subcoil = coil->SubCoils[isc];
			cg[ic] += -subcoil->CurrentFraction * Green(ss->X, ss->Z, subcoil->X, subcoil->Z);
			if (sym == UpDownSymmetric)
				cg[ic] += -subcoil->CurrentFraction * Green(ss->X, ss->Z, subcoil->X, -(subcoil->Z));
		}
	}

	ss->PlasmaGreen = pg;
	ss->CoilGreen = cg;
}

/*
**
**	LoadShellGreens
**
*/
void          LoadShellGreens(TOKAMAK * td)
{
	FILE         *fi;
	SHELL        *m;
	SUBSHELL     *ss;
	size_t        c;
	int           nmax, ncoil, ngrid, im, is;

	nmax = td->PsiGrid->Nsize;
	ncoil = td->NumCoils;
	ngrid = (nmax + 1) * (nmax + 1);

	if (td->SGreenStatus == ShellGreenOK) {	/* read ShellGreens from file */
		fi = fopen(td->SGname, "rb");
		if (!fi)
			nrerror("ERROR:	Could not open ShellGreen file for reading.");

		/* R E A D    S H E L L   G R E E N S */
		for (im = 0; im < td->NumShells; im++) {
			m = td->Shells[im];
			if (m->Enabled) {
				for (is = 0; is < m->NumSubShells; is++) {
					ss = m->SubShells[is];
					ss->CoilGreen = dvector(0, ncoil - 1);
					c = fread(ss->CoilGreen, sizeof(double), (size_t) ncoil, fi);
					cchk(c, ncoil);
					ss->PlasmaGreen = dmatrix(0, nmax, 0, nmax);
					c = fread_dmatrix(ss->PlasmaGreen, 0, nmax, 0, nmax, fi);
					cchk(c, ngrid);
				}
			}
		}

		fclose(fi);
	} else
		RewriteShellGreens(td);
}

/*
**
**	RewriteShellGreens
**
*/
void          RewriteShellGreens(TOKAMAK * td)
{
	FILE         *fi;
	SHELL        *m;
	SUBSHELL     *ss;
	size_t        c;
	int           nmax, ncoil, nshell, ngrid, im, is;

	nmax = td->PsiGrid->Nsize;
	ncoil = td->NumCoils;
	nshell = td->NumShells;
	ngrid = (nmax + 1) * (nmax + 1);

	printf("INFO:	Finding Shell Coupling Greens...\n");
	fprintf(LogFile, "INFO:	Finding Shell Coupling Greens...\n");
	for (im = 0; im < nshell; im++) {
		m = td->Shells[im];
		if (m->Enabled) {
			printf("		[%s: ]\n", m->Name);
			fprintf(LogFile, "		[%s]\n", m->Name);
			for (is = 0; is < m->NumSubShells; is++) {
				ss = m->SubShells[is];
				Find_SubShellGreens(td, ss);
				printf("		[%s]\n", ss->Name);
				fprintf(LogFile, "		[%s]\n", ss->Name);
			}
		}
	}

	fi = fopen(td->SGname, "wb");
	if (!fi)
		nrerror("ERROR:	Could not open file for writing in ShellGreen.");

	/* W R I T E    M E A S   G R E E N S */
	for (im = 0; im < td->NumShells; im++) {
		m = td->Shells[im];
		if (m->Enabled) {
			for (is = 0; is < m->NumSubShells; is++) {
				ss = m->SubShells[is];
				c = fwrite(ss->CoilGreen, sizeof(double), (size_t) ncoil, fi);
				cchk(c, ncoil);
				c = fwrite_dmatrix(ss->PlasmaGreen, 0, nmax, 0, nmax, fi);
				cchk(c, ngrid);
			}
		}
	}

	fclose(fi);

	td->SGreenStatus = ShellGreenOK;
}

/*
**
**	free_ShellGreens
**
*/
void          free_ShellGreens(TOKAMAK * td)
{
	int           im, is, nmax, nshell, ncoil;
	SHELL        *m;
	SUBSHELL     *ss;

	nmax = td->PsiGrid->Nsize;
	nshell = td->NumShells;
	ncoil = td->NumCoils;

	for (im = 0; im < nshell; im++) {
		m = td->Shells[im];
		for (is = 0; is < m->NumSubShells; is++) {
			ss = m->SubShells[is];
			if (ss->CoilGreen)
				free_dvector(ss->CoilGreen, 0, ncoil - 1);
			if (ss->PlasmaGreen)
				free_dmatrix(ss->PlasmaGreen, 0, nmax, 0, nmax);
			ss->CoilGreen = NULL;
			ss->PlasmaGreen = NULL;
		}
	}
}
