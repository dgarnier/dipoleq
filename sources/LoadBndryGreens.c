/*
** TokaMac v2.0
**
** LoadBndryGreens.c
**
** Reads the LH Green's functions & CoilGreens from disk.
** These are used to compute the PsiBoundary.
** Also, it rewrites the file (erasing if it already exists),
** by calling RewriteLHGreens.
**
** File:		LoadBndryGreens.c
** Date:		March 21, 1993
**
** Revised:
**
**		August 5, 1993		Added perfectly conducting shells
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdlib.h>
#include "VAX_Alloc.h"
#include <stdio.h>
#include "nrutil.h"
#include "coil.h"
#include "shell.h"
#include "green.h"
#include "tokgreen.h"
#include "tokamak.h"
#include "PsiBoundary.h"
#include "LoadBndryGreens.h"
#include "multitask.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

#define cchk(c,n)	if ((c) != (n + 1)) nrerror("ERROR: Could not read Boundary Greens.")

extern FILE  *LogFile;


#ifdef __cplusplus
extern "C" {
#endif

void          GetBoundaryGreen(LHVEC * lhv, PSIGRID * pg, double x, double z);
void          GetLHGreen(LHVEC * lhv, PSIGRID * pg, int ix, int iz);
void          FindLHGreen(LHARY * lha, PSIGRID * pg);
void          GetCoilGreen(PSIGRID * pg, COIL * c);
void          GetSubShellGreen(PSIGRID * pg, SUBSHELL * sc);
void          GetShellGreen(PSIGRID * pg, SHELL * s);
void          MakeCoilGreenSymmetric(PSIGRID * pg, COIL * c);
void          MakeSubShellGreenSymmetric(PSIGRID * pg, SUBSHELL * c);
void          MakeShellGreenSymmetric(PSIGRID * pg, SHELL * s);



#ifdef __cplusplus
}
#endif


/*
**
**	GetBoundaryGreen
**
** 	For a given (x,z), fill the lhvec by looping around the
** 	computational domain.
**
**
*/
void          GetBoundaryGreen(LHVEC * lhv, PSIGRID * pg, double x, double z)
{
	int           nmax;
	int           ixp, izp;
	double       *Xp, *Zp;

	nmax = pg->Nsize;
	Xp = pg->X;
	Zp = pg->Z;

	/*Top*/
	for (ixp = 0; ixp <= nmax; ixp++)
		lhv->Top[ixp] = Green(x, z, Xp[ixp], Zp[nmax]);

	/*Bottom*/
	for (ixp = 0; ixp <= nmax; ixp++)
		lhv->Bot[ixp] = Green(x, z, Xp[ixp], Zp[0]);

	/*Outside*/
	for (izp = 0; izp <= nmax; izp++)
		lhv->Out[izp] = Green(x, z, Xp[nmax], Zp[izp]);

	/*Inside*/
	for (izp = 0; izp <= nmax; izp++)
		lhv->In[izp] = Green(x, z, Xp[0], Zp[izp]);
}

/*
**
**	GetLHGreen
**
*/
void          GetLHGreen(LHVEC * lhv, PSIGRID * pg, int ix, int iz)
{
	int           nmax;
	int           ixp, izp;
	double        dx, dz, temp;
	double       *Xp, *Zp;

	nmax = pg->Nsize;
	dx = pg->dx;
	dz = pg->dz;
	Xp = pg->X;
	Zp = pg->Z;

	GetBoundaryGreen(lhv, pg, Xp[ix], Zp[iz]);

	/*Top*/
	temp = dx / TWOPI;
	for (ixp = 0; ixp <= nmax; ixp++)
		lhv->Top[ixp] *= temp / Xp[ixp];

	/*Bottom*/
	temp = dx / TWOPI;
	for (ixp = 0; ixp <= nmax; ixp++)
		lhv->Bot[ixp] *= temp / Xp[ixp];

	/*Outside*/
	temp = dz / (TWOPI * Xp[nmax]);
	for (izp = 0; izp <= nmax; izp++)
		lhv->Out[izp] *= temp;

	/*Inside*/
	temp = dz / (TWOPI * Xp[0]);
	for (izp = 0; izp <= nmax; izp++)
		lhv->In[izp] *= temp;
}

/*
**
**	FindLHGreen
**
** 	This is an interesting routine because of the book-keeping involved.
** 	We recognise that for a rectangular computational domain, the Green's
** 	functions are ALWAYS up/down symmetric...even if the plasma and coils
** 	are not.  Thus, we compute the bottom and the bottom half of the
** 	insides and the outsides of the computational domains..
**
**
*/
void          FindLHGreen(LHARY * lha, PSIGRID * pg)
{
	int           ix, iz, nmax;

	nmax = pg->Nsize;

	printf("INFO:	Finding LH Greens....\n");
	fprintf(LogFile, "INFO:	Finding LH Greens....\n");

	for (ix = 0; ix <= nmax; ix++)
		GetLHGreen(lha->Bot[ix], pg, ix, 0);
	printf("		[Bot]\n");
	fprintf(LogFile, "		[Bot]\n");

	for (iz = 0; iz <= nmax / 2; iz++)
		GetLHGreen(lha->In[iz], pg, 0, iz);
	printf("		[In]\n");
	fprintf(LogFile, "		[In]\n");

	for (iz = 0; iz <= nmax / 2; iz++)
		GetLHGreen(lha->Out[iz], pg, nmax, iz);
	printf("		[Out]\n");
	fprintf(LogFile, "		[Out]\n");
}

/*
**
**	GetCoilGreen
**
** 	Fills all of the Greens functions for the external coils...
**
** 	** This gets the actual Green's Function multiplied by TWOPI.
**
**
*/
void          GetCoilGreen(PSIGRID * pg, COIL * c)
{
	int           isc, ix, iz, nmax, nsc;
	COILGREEN    *cgt, *cg;
	SUBCOIL      *sc;

	nmax = pg->Nsize;
	nsc = c->NumSubCoils;

	cg = c->CoilGreen;

	/* initialize cg to 0.0 */
	for (ix = 0; ix <= nmax; ix++)
		cg->Top[ix] = cg->Bot[ix] = 0.0;
	for (iz = 0; iz <= nmax; iz++)
		cg->In[iz] = cg->Out[iz] = 0.0;

	/* create a temporary cgt */
	cgt = new_CoilGreen(nmax);

	for (isc = 0; isc < nsc; isc++) {
		sc = c->SubCoils[isc];
		GetBoundaryGreen((LHVEC *) cgt, pg, sc->X, sc->Z);
		/* Top */
		for (ix = 0; ix <= nmax; ix++)
			cg->Top[ix] += cgt->Top[ix] * sc->CurrentFraction;
		/* Bot */
		for (ix = 0; ix <= nmax; ix++)
			cg->Bot[ix] += cgt->Bot[ix] * sc->CurrentFraction;
		/* Inside */
		for (iz = 0; iz <= nmax; iz++)
			cg->In[iz] += cgt->In[iz] * sc->CurrentFraction;
		/* Outside */
		for (iz = 0; iz <= nmax; iz++)
			cg->Out[iz] += cgt->Out[iz] * sc->CurrentFraction;
	}

	/* destroy the temporary cgt */
	free_CoilGreen(cgt, nmax);
}

/*
**
**	GetSubShellGreen
**
** 	Fills all of the Greens functions for the external coils...
**
** 	** This gets the actual Green's Function multiplied by TWOPI.
**
**
*/
void          GetSubShellGreen(PSIGRID * pg, SUBSHELL * sc)
{
	int           ix, iz, nmax;
	COILGREEN    *cgt;
	SHELLGREEN   *cg;

	nmax = pg->Nsize;

	cg = sc->ShellGreen;

	/* initialize cg to 0.0 */
	for (ix = 0; ix <= nmax; ix++)
		cg->Top[ix] = cg->Bot[ix] = 0.0;
	for (iz = 0; iz <= nmax; iz++)
		cg->In[iz] = cg->Out[iz] = 0.0;

	/* create a temporary cgt */
	cgt = new_CoilGreen(nmax);

	GetBoundaryGreen((LHVEC *) cgt, pg, sc->X, sc->Z);
	/* Top */
	for (ix = 0; ix <= nmax; ix++)
		cg->Top[ix] += cgt->Top[ix];
	/* Bot */
	for (ix = 0; ix <= nmax; ix++)
		cg->Bot[ix] += cgt->Bot[ix];
	/* Inside */
	for (iz = 0; iz <= nmax; iz++)
		cg->In[iz] += cgt->In[iz];
	/* Outside */
	for (iz = 0; iz <= nmax; iz++)
		cg->Out[iz] += cgt->Out[iz];

	/* destroy the temporary cgt */
	free_CoilGreen(cgt, nmax);
}

/*
**
**	GetShellGreen
**
**
*/
void          GetShellGreen(PSIGRID * pg, SHELL * s)
{
	int           i;
	SUBSHELL     *subshell;

	for (i = 0; i < s->NumSubShells; i++) {
		subshell = s->SubShells[i];
		/*
		printf("		[%s]\n", subshell->Name);
		fprintf(LogFile, "		[%s]\n", subshell->Name);
		*/
		GetSubShellGreen(pg, subshell);
	}
}

/*
**
**	MakeCoilGreenSymmetric
**
** 	Make the Green's Function for the Coil Set symmetric.
** 	Effectively adds a mirror-image of the subcoils.
**
*/
void          MakeCoilGreenSymmetric(PSIGRID * pg, COIL * c)
{
	int           ix, iz, nmax;
	COILGREEN    *cg;

	nmax = pg->Nsize;

	cg = c->CoilGreen;

	/*Top and Bottom*/
	for (ix = 0; ix <= nmax; ix++) {
		cg->Top[ix] += cg->Bot[ix];
		cg->Bot[ix] = cg->Top[ix];
	}

	/*In and Out*/
	for (iz = 0; iz <= nmax / 2; iz++) {
		cg->In[iz] += cg->In[nmax - iz];
		cg->In[nmax - iz] = cg->In[iz];
		cg->Out[iz] += cg->Out[nmax - iz];
		cg->Out[nmax - iz] = cg->Out[iz];
	}
}

/*
**
**	MakeSubShellGreenSymmetric
**
** 	Make the Green's Function for a subshell symmetric.
** 	Effectively adds a mirror-image of the subshell.
**
*/
void          MakeSubShellGreenSymmetric(PSIGRID * pg, SUBSHELL * c)
{
	int           ix, iz, nmax;
	SHELLGREEN   *cg;

	nmax = pg->Nsize;

	cg = c->ShellGreen;

	/*Top and Bottom*/
	for (ix = 0; ix <= nmax; ix++) {
		cg->Top[ix] += cg->Bot[ix];
		cg->Bot[ix] = cg->Top[ix];
	}

	/*In and Out*/
	for (iz = 0; iz <= nmax / 2; iz++) {
		cg->In[iz] += cg->In[nmax - iz];
		cg->In[nmax - iz] = cg->In[iz];
		cg->Out[iz] += cg->Out[nmax - iz];
		cg->Out[nmax - iz] = cg->Out[iz];
	}
}

/*
**
**	MakeShellGreenSymmetric
**
**
*/
void          MakeShellGreenSymmetric(PSIGRID * pg, SHELL * s)
{
	int           i;
	SUBSHELL     *subshell;

	for (i = 0; i < s->NumSubShells; i++) {
		subshell = s->SubShells[i];
		MakeSubShellGreenSymmetric(pg, subshell);
	}
}

/*
**
**	LoadBndryGreens
**
*/
void          LoadBndryGreens(TOKAMAK * td)
{
	FILE         *fi;
	LHARY        *lhpg;
	LHVEC        *lhv;
	COIL         *aCoil;
	SHELL        *shell;
	SUBSHELL     *subshell;
	size_t        c;
        unsigned int  nmax;
	int           i, iss;

	nmax = td->PsiGrid->Nsize;

	/* A L L O C A T E   S P A C E */
	lhpg = new_LHary(nmax);
	td->LHPlasmaGreen = lhpg;

	for (i = 0; i < td->NumCoils; i++) {
		aCoil = td->Coils[i];
		if (!(aCoil->CoilGreen))
			aCoil->CoilGreen = new_CoilGreen(nmax);
	}

	for (i = 0; i < td->NumShells; i++) {
		shell = td->Shells[i];
		for (iss = 0; iss < shell->NumSubShells; iss++) {
			subshell = shell->SubShells[iss];
			if (!(subshell->ShellGreen))
				subshell->ShellGreen = new_ShellGreen(nmax);
		}
	}

	if (td->LHGreenStatus == LHGreenOK) {	/* read BndryGreens from file */
		fi = fopen(td->LHname, "rb");
		if (!fi)
			nrerror("ERROR:	Could not open BndryGreens file for reading.");
MULTI;

		/* R E A D    L H G R E E N S */
		for (i = nmax / 2; i >= 0; i--) {
			lhv = lhpg->In[i];
			c = fread(lhv->Top, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
			c = fread(lhv->Bot, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
			c = fread(lhv->In, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
			c = fread(lhv->Out, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
		}
MULTI;

		for (i = nmax / 2; i >= 0; i--) {
			lhv = lhpg->Out[i];
			c = fread(lhv->Top, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
			c = fread(lhv->Bot, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
			c = fread(lhv->In, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
			c = fread(lhv->Out, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
		}
MULTI;
		for (i = nmax; i >= 0; i--) {
			lhv = lhpg->Bot[i];
			c = fread(lhv->Top, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
			c = fread(lhv->Bot, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
			c = fread(lhv->In, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
			c = fread(lhv->Out, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
		}
MULTI;

		/* R E A D    C O I L G R E E N S */
		for (i = 0; i < td->NumCoils; i++) {
			aCoil = td->Coils[i];
			c = fread(aCoil->CoilGreen->Top, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
			c = fread(aCoil->CoilGreen->Bot, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
			c = fread(aCoil->CoilGreen->In, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
			c = fread(aCoil->CoilGreen->Out, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
		}

MULTI;
		/* R E A D    S H E L L G R E E N S */
		for (i = 0; i < td->NumShells; i++) {
			shell = td->Shells[i];
			for (iss = 0; iss < shell->NumSubShells; iss++) {
				subshell = shell->SubShells[iss];
				c = fread(subshell->ShellGreen->Top, sizeof(double), nmax + 1, fi);
				cchk(c, nmax);
				c = fread(subshell->ShellGreen->Bot, sizeof(double), nmax + 1, fi);
				cchk(c, nmax);
				c = fread(subshell->ShellGreen->In, sizeof(double), nmax + 1, fi);
				cchk(c, nmax);
				c = fread(subshell->ShellGreen->Out, sizeof(double), nmax + 1, fi);
				cchk(c, nmax);
			}
		}
MULTI;

		fclose(fi);
	} else
		RewriteBndryGreens(td);
}

/*
**
**	RewriteBndryGreens
**
*/
void          RewriteBndryGreens(TOKAMAK * td)
{
	FILE         *fi;
	LHARY        *lhpg;
	LHVEC        *lhv;
	COIL         *aCoil;
	SHELL        *shell;
	SUBSHELL     *subshell;
	size_t        c;
        unsigned int nmax;
	int          i, iss;

	nmax = td->PsiGrid->Nsize;
	lhpg = td->LHPlasmaGreen;

	FindLHGreen(lhpg, td->PsiGrid);

	printf("INFO:	Finding Coil Boundary Greens...\n");
	fprintf(LogFile, "INFO:	Finding Coil Boundary Greens...\n");
	for (i = 0; i < td->NumCoils; i++) {
		aCoil = td->Coils[i];
		printf("		[%s]\n", aCoil->Name);
		fprintf(LogFile, "		[%s]\n", aCoil->Name);
		GetCoilGreen(td->PsiGrid, aCoil);
		if (td->PsiGrid->Symmetric)
			MakeCoilGreenSymmetric(td->PsiGrid, aCoil);
	}
    MULTI;

	if (td->NumShells > 0) {
		printf("INFO:	Finding Shell Boundary Greens...\n");
		fprintf(LogFile, "INFO:	Finding Shell Boundary Greens...\n");
		for (i = 0; i < td->NumShells; i++) {
			shell = td->Shells[i];
			printf("		[%s]\n", shell->Name);
			fprintf(LogFile, "		[%s]\n", shell->Name);
			GetShellGreen(td->PsiGrid, shell);
			if (td->PsiGrid->Symmetric)
				MakeShellGreenSymmetric(td->PsiGrid, shell);
		}
	}
	if ((td->LHname[0] == '\0') || (td->LHname[0] == '*')) {
		// don't write the greens file
		return;
	}

	fi = fopen(td->LHname, "wb");
	if (!fi)
		nrerror("ERROR:	Could not open BndryGreen file for writing.");

	/* W R I T E    L H G R E E N S */
	for (i = nmax / 2; i >= 0; i--) {
		lhv = lhpg->In[i];
		c = fwrite(lhv->Top, sizeof(double), nmax + 1, fi);
		cchk(c, nmax);
		c = fwrite(lhv->Bot, sizeof(double), nmax + 1, fi);
		cchk(c, nmax);
		c = fwrite(lhv->In, sizeof(double), nmax + 1, fi);
		cchk(c, nmax);
		c = fwrite(lhv->Out, sizeof(double), nmax + 1, fi);
		cchk(c, nmax);
	}
    MULTI;

	for (i = nmax / 2; i >= 0; i--) {
		lhv = lhpg->Out[i];
		c = fwrite(lhv->Top, sizeof(double), nmax + 1, fi);
		cchk(c, nmax);
		c = fwrite(lhv->Bot, sizeof(double), nmax + 1, fi);
		cchk(c, nmax);
		c = fwrite(lhv->In, sizeof(double), nmax + 1, fi);
		cchk(c, nmax);
		c = fwrite(lhv->Out, sizeof(double), nmax + 1, fi);
		cchk(c, nmax);
	}
    MULTI;
	for (i = nmax; i >= 0; i--) {
		lhv = lhpg->Bot[i];
		c = fwrite(lhv->Top, sizeof(double), nmax + 1, fi);
		cchk(c, nmax);
		c = fwrite(lhv->Bot, sizeof(double), nmax + 1, fi);
		cchk(c, nmax);
		c = fwrite(lhv->In, sizeof(double), nmax + 1, fi);
		cchk(c, nmax);
		c = fwrite(lhv->Out, sizeof(double), nmax + 1, fi);
		cchk(c, nmax);
	}
    MULTI;

	/*  W R I T E    C O I L    G R E E N S */
	for (i = 0; i < td->NumCoils; i++) {
		aCoil = td->Coils[i];
		c = fwrite(aCoil->CoilGreen->Top, sizeof(double), nmax + 1, fi);
		cchk(c, nmax);
		c = fwrite(aCoil->CoilGreen->Bot, sizeof(double), nmax + 1, fi);
		cchk(c, nmax);
		c = fwrite(aCoil->CoilGreen->In, sizeof(double), nmax + 1, fi);
		cchk(c, nmax);
		c = fwrite(aCoil->CoilGreen->Out, sizeof(double), nmax + 1, fi);
		cchk(c, nmax);
	}

    MULTI;
	/*  W R I T E    S H E L L    G R E E N S */
	for (i = 0; i < td->NumShells; i++) {
		shell = td->Shells[i];
		for (iss = 0; iss < shell->NumSubShells; iss++) {
			subshell = shell->SubShells[iss];
			c = fwrite(subshell->ShellGreen->Top, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
			c = fwrite(subshell->ShellGreen->Bot, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
			c = fwrite(subshell->ShellGreen->In, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
			c = fwrite(subshell->ShellGreen->Out, sizeof(double), nmax + 1, fi);
			cchk(c, nmax);
		}
	}
    MULTI;

	fclose(fi);
	td->LHGreenStatus = LHGreenOK;
}

/*
**
**	free_BndryGreens
**
*/
void          free_BndryGreens(TOKAMAK * td)
{
	int           i, iss, nmax;
	COIL         *aCoil;
	SHELL        *shell;
	SUBSHELL     *subshell;

	nmax = td->PsiGrid->Nsize;

	if (td->LHPlasmaGreen)
		free_LHary(td->LHPlasmaGreen, nmax);
	td->LHPlasmaGreen = NULL;

	for (i = 0; i < td->NumCoils; i++) {
		aCoil = td->Coils[i];
		if (aCoil->CoilGreen)
			free_CoilGreen(aCoil->CoilGreen, nmax);
		td->Coils[i]->CoilGreen = NULL;
	}
    MULTI;

	for (i = 0; i < td->NumShells; i++) {
		shell = td->Shells[i];
		for (iss = 0; iss < shell->NumSubShells; iss++) {
			subshell = shell->SubShells[iss];
			if (subshell->ShellGreen) {
				free_ShellGreen(subshell->ShellGreen, nmax);
				subshell->ShellGreen = NULL;
			}
		}
	}
}
