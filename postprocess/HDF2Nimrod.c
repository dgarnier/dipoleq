
#include <stdio.h>
#include <string.h>
#ifdef __MWERKS__
#include <console.h>
#include <SIOUX.h>
#endif
#ifdef __SC__
#include <console.h>
#endif
#ifdef VAXC
#include <climsgdef.h>
#include <descrip.h>
#endif
#include "nrutil.h"
#include "psigrid.h"
#include "plasma.h"
#include "HDFInput.h"
#include "interpolate.h"
#include "contour.h"

#define TRUE	1
#define FALSE	0
#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307


/* static int		gCount;
static double   *gXvec;
static double   *gZvec; */

void	CountContourStep(double x, double z, double, int flag);
void	RecordContourStep(double x, double z, double, int flag);
int 	GetContour(PSIGRID *pg, double Psi, double **X, double **Z);

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307


int main(int argc, char **argv)
{
	int ipts, jpts, i, j;

	double **P, **Psi, *R, *Z, **dPr, **dPz, **J ;

	PSIGRID 	*pg;
	PLASMA 		*pl;

	char         fn[256] = "dipoleq.hdf";
        char         on[256];
        char		*s;
	FILE         *fi;
        int		outset = 0;

#ifdef __MWERKS__
	unsigned char wtitle[256] = "\pDipole HDF to Text";
#endif

	/* F I L E   I N P U T */

#ifdef __SC__
	argc = ccommand(&argv);
#endif
#ifdef __MWERKS__
	SIOUXSettings.asktosaveonclose = 0;
    SIOUXSettings.autocloseonquit = 0;  /* close this baby when you quit! */
    SIOUXSettings.showstatusline = TRUE;
	SIOUXSetTitle(wtitle);
	argc = ccommand(&argv);
#endif

	/*
	** VAX VMS does not have a standard command line interface
	**
	** For UNIX and other standard C systems, the C command-line
	** interface should be fine.
	**
	** hdf2nimrod input output
	**
	** ignores qualifiers..
	*/
#ifndef VAXC
	if (argc < 2) {
            printf("HDF2Nimrod usage:\n");
            printf("	hdf2nimrod dipoleq_output[.hdf] [nimrod_input[.dat]]\n");
        } else {
            for (i = 1; i < argc; i++)
                if (argv[i][0] != '-') {
                    strncpy(fn,argv[i],256);
                    continue;
                }
            for ( ; i < argc; i++)
                if (argv[i][0] != '-') {
                    strncpy(on,argv[i],256);
                    if (strchr(on,'.') == NULL) strcat(on,".dat");
                    outset = 1;
                    continue;
                }
            if (!outset) {
                strncpy(on,fn,256);
                if ((s=strchr(on,'.')) != NULL) *s = '\0';
                strcat(on,"_nimrod.dat");
            }
            if (strchr(fn,'.') == NULL) strcat(fn,".hdf");
        }
#else
	/*
	**	For VAX VMS systems, we use command line information (CLI)
	**	routines.  A typical input line might look like
	**
	**	$ TokaMac /in=TokIn.dat /log=tlog.out
	**
	**	The qualifiers are optional.
	*/

	int           cli_status;
	char          cli_value[255];
	short         cli_val_len;
	char          vax_fname[63];
	char          vax_lgname[63];

	struct dsc$descriptor_s cli_value_str =
	{255, 14, 1, cli_value};
	$DESCRIPTOR(cli_infile, "IN");
	$DESCRIPTOR(cli_logfile, "LOG");

	cli_status = cli$present(&cli_infile);	/* is /INFILE present ? */
	if (cli_status == CLI$_PRESENT) {
		cli_status = cli$get_value(&cli_infile, &cli_value_str,
&cli_val_len);
		strncpy(vax_fname, cli_value, cli_val_len);
		fn = vax_fname;
	}
	cli_status = cli$present(&cli_logfile);	/* is /LOG present ? */
	if (cli_status == CLI$_PRESENT) {
		cli_status = cli$get_value(&cli_infile, &cli_value_str,
&cli_val_len);
		strncpy(vax_lgname, cli_value, cli_val_len);
		lgn = vax_lgname;
	}
#endif
	printf ("Creates a text file from the 2D data in HDF\n\n");


	pg = HDFPsiGridIn(fn);
	pl = HDFPlasmaIn(pg ,fn);
//	HDFFluxFuncsIn(fn, &npts, &PsiX, &Psi, &P, &G, &Pp, &G2p,
//					&q, &dVdpsi, &Vol, &Shear, &Well, &Jave, &B2ave, &Beta);


	ipts = jpts = pg->Nsize+1;
	R = pg->X;
	Z = pg->Z;
	Psi = pg->Psi;
	dPr = pl->GradPsiX;
	dPz = pl->GradPsiZ;
	P = pl->Piso;
	J = pg->Current;

        fi = fopen(on, "w");
        if (!fi)
            nrerror("ERROR:	Could not open file for writing.");

	    fprintf(fi,"NIMROD input.  From Input File: %s\n",fn);

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
	printf("We are out of here!\n");
	return 0;
}
