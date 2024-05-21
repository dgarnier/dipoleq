
#include <stdio.h>
#include <string.h>
#include "nrutil.h"
#include "psigrid.h"
#include "plasma.h"
#include "tokamak.h"
#include "FileOutput.h"

#define TRUE	1
#define FALSE	0
#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307


void NimrodOutput(TOKAMAK *td)
{
	int ipts, jpts, i, j;

	double **P, **Psi, *R, *Z, **dPr, **dPz, **J ;

	PSIGRID *pg;
	PLASMA 	*pl;

	char	fname[256];
	FILE	*fi;
	int		outset = 0;

	pl = td->Plasma;
	pg = td->PsiGrid;

	strncat(fname, td->Oname, 18);	/* take 1st 18 characters */
	strcat(fname, "_nimrod.txt");

	ipts = jpts = pg->Nsize+1;
	R = pg->X;
	Z = pg->Z;
	Psi = pg->Psi;
	dPr = pl->GradPsiX;
	dPz = pl->GradPsiZ;
	P = pl->Piso;
	J = pg->Current;

	fi = fopen(fname, "w");
    if (!fi)
            nrerror("ERROR:	Could not open file for writing.");

		fprintf(fi,"%12d\t\tNumber of Z points\n",ipts);
		fprintf(fi,"%12d\t\tNumber of R points\n",jpts);

	fprintf(fi,"         R      \t     Z      \t  Psi   \t  dPsi/dx    \t  dPsi/dz    \t      P      \t      J\n");

	for (i=0;i<ipts;i++) {
		for (j=0;j<jpts;j++) {
			fprintf(fi,"%21.15e\t%21.15e\t%21.15e\t%21.15e\t%21.15e\t%21.15e\t%21.15e\n",R[i],Z[j],
			    Psi[i][j],
			    dPr[i][j],  dPz[i][j],
			    P[i][j], J[i][j]/MU0);
		}
	}
	fclose(fi);
}
