/*
** TokaMac v2.0
**
** Restart.c
**
** Routines to read and write the restart file.
** Only the PsiGrid->Current and PsiGrid->Psi arrays
** are writtten.
**
** File:		Restart.c
** Date:		April 2, 1993
**
** Modifications:
**
**	Oct. 3, 1993		Saved solution as well as PSIGRID.
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdio.h>
#include "nrutil.h"
#include "stdio_dmatrix.h"
#include "psigrid.h"
#include "tokamak.h"
#include "Restart.h"
#include "InitJ.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

#define cchk(c,n)	if ((c) != (n + 1)) nrerror("ERROR: Could not read/write RestartFile.")

extern FILE  *LogFile;

void          ReadRestart(char *fn, TOKAMAK * td)
{
	FILE         *fi;
	size_t        count1, count2;
	int           nmax, ix, iz, i;
	double        sum = 0.0, **J;
	PSIGRID      *pg;
	PLASMA       *pl;

	pl = td->Plasma;
	pg = td->PsiGrid;
	nmax = pg->Nsize;

	printf("INFO:	Restarting current from %s.\n", fn);
	fprintf(LogFile, "INFO:	Restarting current from %s.\n", fn);

	fi = fopen(fn, "rb");
        if (!fi) {
            nrinfo("WARNING:	Could not open restart file.");
            InitJ(td->PsiGrid, td->Plasma);
            return;
        }

	/* R E A D    P S I G R I D */

	count1 = fread_dmatrix(pg->Current, 0, nmax, 0, nmax, fi);

	count2 = fread_dmatrix(pg->Psi, 0, nmax, 0, nmax, fi);

	if ((count1 != (nmax + 1) * (nmax + 1)) || (count2 != (nmax + 1) * (nmax + 1))) {
		fclose(fi);
		nrerror("ERROR: Could not read Restart file.");
	}
	J = pg->Current;
	for (ix = 0; ix <= nmax; ix++)
		for (iz = 0; iz <= nmax; iz++)
			sum += J[ix][iz];

	sum *= pg->dx * pg->dz;
	sum /= MU0;

	/* R E A D   P L A S M A   U N K N O W N S */

	if (td->RestartUnkns == RestartOK) {
		count1 = fread(pl->G2p, sizeof(double), MaxPolyTerms + 1, fi);
		cchk(count1, MaxPolyTerms);
		count1 = fread(pl->H, sizeof(double), MaxPolyTerms + 1, fi);
		cchk(count1, MaxPolyTerms);
		count1 = fread(pl->Pp, sizeof(double), MaxPolyTerms + 1, fi);
		cchk(count1, MaxPolyTerms);
		count1 = fread(pl->Rot, sizeof(double), MaxPolyTerms + 1, fi);
		cchk(count1, MaxPolyTerms);
		count1 = fread(pl->Siso, sizeof(double), MaxPolyTerms + 1, fi);
		cchk(count1, MaxPolyTerms);
		count1 = fread(pl->Spar, sizeof(double), MaxPolyTerms + 1, fi);
		cchk(count1, MaxPolyTerms);
		count1 = fread(pl->Sper, sizeof(double), MaxPolyTerms + 1, fi);
		cchk(count1, MaxPolyTerms);

		/* R E A D   C O I L   C U R R E N T S */

		for (i = 0; i < td->NumCoils; i++) {
			count1 = fread(&(td->Coils[i]->CoilCurrent), sizeof(double), 1, fi);
			cchk(count1, 0);
		}
	}
	fclose(fi);

	if (td->RestartUnkns) {
		printf("		The previously saved unknowns were restored.\n");
		fprintf(LogFile, "		The previously saved unknowns were restored.\n");
	} else {
		printf("		The previously saved unknowns were NOT restored.\n");
		fprintf(LogFile, "		The previously saved unknowns were NOT restored.\n");
	}

	printf("		[Ip = %g (A)]\n", sum);
	fprintf(LogFile, "		[Ip = %g (A)]\n", sum);
}

void          WriteRestart(char *fn, TOKAMAK * td)
{
	FILE         *fi;
	size_t        count1, count2;
	int           nmax, i;
	PSIGRID      *pg;
	PLASMA       *pl;

	pl = td->Plasma;
	pg = td->PsiGrid;
	nmax = pg->Nsize;

	fi = fopen(fn, "wb");
	if (!fi)
		nrerror("ERROR:	Could not write restart file.");

	/* S A V E   P S I G R I D  */

	count1 = fwrite_dmatrix(pg->Current, 0, nmax, 0, nmax, fi);

	count2 = fwrite_dmatrix(pg->Psi, 0, nmax, 0, nmax, fi);

	if ((count1 != (nmax + 1) * (nmax + 1)) || (count2 != (nmax + 1) * (nmax + 1))) {
		fclose(fi);
		nrerror("ERROR: Could not write Restart file.");
	}
	/* S A V E   P L A S M A    U N K N O W N S */

	count1 = fwrite(pl->G2p, sizeof(double), MaxPolyTerms + 1, fi);
	cchk(count1, MaxPolyTerms);
	count1 = fwrite(pl->H, sizeof(double), MaxPolyTerms + 1, fi);
	cchk(count1, MaxPolyTerms);
	count1 = fwrite(pl->Pp, sizeof(double), MaxPolyTerms + 1, fi);
	cchk(count1, MaxPolyTerms);
	count1 = fwrite(pl->Rot, sizeof(double), MaxPolyTerms + 1, fi);
	cchk(count1, MaxPolyTerms);
	count1 = fwrite(pl->Siso, sizeof(double), MaxPolyTerms + 1, fi);
	cchk(count1, MaxPolyTerms);
	count1 = fwrite(pl->Spar, sizeof(double), MaxPolyTerms + 1, fi);
	cchk(count1, MaxPolyTerms);
	count1 = fwrite(pl->Sper, sizeof(double), MaxPolyTerms + 1, fi);
	cchk(count1, MaxPolyTerms);

	/* S A V E   C O I L   C U R R E N T S */

	for (i = 0; i < td->NumCoils; i++) {
		count1 = fwrite(&(td->Coils[i]->CoilCurrent), sizeof(double), 1, fi);
		cchk(count1, 0);
	}

	fclose(fi);

	printf("INFO:	Restart file written to %s.\n", fn);
	fprintf(LogFile, "INFO:	Restart file written to %s.\n", fn);

}

#undef cchk
