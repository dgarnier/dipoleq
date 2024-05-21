/*
 *  Hello World for the CodeWarrior
 *  © 1997-1998 Metrowerks Corp.
 *
 *  Questions and comments to:
 *       <mailto:support@metrowerks.com>
 *       <http://www.metrowerks.com/>
 */

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


static int		gCount;
static double   *gXvec;
static double   *gZvec;

void	CountContourStep(double x, double z, double, int flag);
void	RecordContourStep(double x, double z, double, int flag);
int 	GetContour(PSIGRID *pg, double Psi, double **X, double **Z);
void ComputeCoilCurrentCentroid(char *fn, FILE *out);

// this routine finds the current centroid from the _Conductors.out file by
// adding up all subcoil currents... ignores coils that are not enabled.
void ComputeCoilCurrentCentroid(char *fn, FILE *out)
{
    char	con[256];
    char	*s;
    FILE	*in;
    float	x,z,cur;
    double	x0,z0,current;
    int		n,enabled;

    strncpy(con,fn,256);
    if ((s=strchr(con,'.')) != NULL) *s = '\0';
    strcat(con,"_Conductors.out");
    in = fopen(con,"r");

    if (in == NULL) {
        fprintf(out,"Couldn't open conductor file: %s\n",con);
        return;
    }
    con[0] = '\0';
    x0 = z0 = current = 0.0;
    enabled = 0;
    while ((n=fscanf(in,"%20s %f %f %f",con,&x,&z,&cur)) != EOF) {
 //       printf("n s x z c %d %s %f %f %f\n",n,con,x,z,cur);
        if (n == 1) {
            if (strcmp(con,"On/Off") == 0) {
                n = fscanf(in,"= %d",&enabled);
//                printf("Scanned On/Off = %d\n",enabled);
            }
        }
        if ((n == 4) && (enabled != 0)) {
            x0 += x*cur;
            z0 += z*cur;
            current += cur;
        }
    }
    fclose(in);
    x0 /= current;
    z0 /= current;
    fprintf(out,"Current Centroid R, Z, and total coil current\n");
    fprintf(out,"%10.4f %10.4f %10.4e\n",x0,z0,current);
    return;
}

int main(int argc, char **argv)
{
	int npts;
	double *PsiX, *Psi,  *P, *G,  *Pp, *G2p, *q,  *dVdpsi,  *Vol, *Shear,
			 *Well, *Jave, *B2ave, *Beta;
	double *XC, *ZC, Px, Pz;

	PSIGRID 	*pg;
	PLASMA 		*pl;

	char         fn[256] = "LDXTest.HDF";
        char	     on[256];
	char         *lgn = "DipLog.out";
        char		*s;
	int           i, j, ncpts;
        int		outset = 0;
	FILE         *fi;

#ifdef __MWERKS__
	unsigned char wtitle[256] = "\pPostProcessor v0.9";
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
	** tokamak -f infile -l logfile
	**
	** where all qualifiers are optional.
	*/
#ifndef VAXC
	if (argc < 2) {
            printf("contourHDF usage:\n");
            printf("	countourhdf dipoleq_output[.hdf] [countour[.dat]]\n");
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
                strcat(on,"_contour.dat");
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
		cli_status = cli$get_value(&cli_infile, &cli_value_str, &cli_val_len);
		strncpy(vax_fname, cli_value, cli_val_len);
		fn = vax_fname;
	}
	cli_status = cli$present(&cli_logfile);	/* is /LOG present ? */
	if (cli_status == CLI$_PRESENT) {
		cli_status = cli$get_value(&cli_infile, &cli_value_str, &cli_val_len);
		strncpy(vax_lgname, cli_value, cli_val_len);
		lgn = vax_lgname;
	}
#endif
	printf ("Handy Dipole post processor\n\n");

        fi = fopen(fn,"r");
        if (fi == NULL) nrerror("ERROR: Could not open HDF file.");
        fclose(fi);

	pg = HDFPsiGridIn(fn);
	pl = HDFPlasmaIn(pg ,fn);
	HDFFluxFuncsIn(fn, &npts, &PsiX, &Psi, &P, &G, &Pp, &G2p,
					&q, &dVdpsi, &Vol, &Shear, &Well, &Jave, &B2ave, &Beta);


        fi = fopen(on, "w");
        if (!fi)
			nrerror("ERROR:	Could not open file for writing.");

	fprintf(fi,"Dipole Contour Tracer.  From Input File: %s\n",fn);
        ComputeCoilCurrentCentroid(fn, fi);
        fprintf(fi,"Number of contours:\n%d\n",npts);
	fprintf(fi,"Cont#\tPsiNorm\tPsi\tPress\tdPdPsi\tV\t# points\n");
	fprintf(fi,"     R      \t     Z      \t  dPsi/dx   \t  dPsi/dz   \n");


	for (i=0;i<npts;i++) {

		printf("Doing Contour Trace %s %d\n",fn,i);

		ncpts = GetContour(pg,Psi[i],&XC,&ZC);

//		fprintf(fi,"%12g\t\tNormalized Psi of Contour\n", PsiX[i]);
//		fprintf(fi,"%12g\t\tPsi (webers) \n",Psi[i]);
//		fprintf(fi,"%12g\t\tp(psi) (pascals) \n", P[i]);
//		fprintf(fi,"%12g\t\tdp/dpsi (pascals/weber)\n",Pp[i]);
//		fprintf(fi,"%12g\t\tV (m^3/weber)\n",dVdpsi[i]);
//		fprintf(fi,"%12d\t\tNumber of contour points\n",ncpts);
//		fprintf(fi,"     R      \t     Z      \t  dPsi/dx   \t  dPsi/dz   \n");
		fprintf(fi,"%5d\t%5f\t%12g\t%12g\t%12g\t%7g\t%5d\n",
			i, PsiX[i], Psi[i], P[i], Pp[i], dVdpsi[i], ncpts);
		for (j=0;j<ncpts;j++) {
			Px = interpolate(pg, pl->GradPsiX, XC[j], ZC[j]);
			Pz = interpolate(pg, pl->GradPsiZ, XC[j], ZC[j]);
			fprintf(fi,"%12g\t%12g\t%12g\t%12g\n",XC[j],ZC[j],Px,Pz);
		}
		free_dvector(XC,0,ncpts);
		free_dvector(ZC,0,ncpts);
	}
	fclose(fi);
	printf("We are out of here!\n");
	return 0;
}


void		CountContourStep(double x, double z, double psi, int flag)
{
#pragma unused(x,z,psi)
	switch (flag) {
	  case CONTOUR_START:
	  	  gCount = 0;
		  break;
	  case CONTOUR_TRACE:
	  case CONTOUR_STOP:
	  	  gCount++;
		  break;
	}
}

void		RecordContourStep(double x, double z, double psi, int flag)
{
#pragma unused(psi)
	switch (flag) {
	  case CONTOUR_START:
//	  	  printf("Starting Contour ... ");
	  	  gCount = 0;
		  break;
	  case CONTOUR_STOP:
//	  	  printf("Finished: %d\n",gCount+1);
	  case CONTOUR_TRACE:
	  	  gCount++;
		  break;
	}
	gXvec[gCount] = x;
	gZvec[gCount] = z;
}

int 	GetContour(PSIGRID *pg, double Psi, double **X, double **Z)
{
		int nmax, npts;

		nmax = pg->Nsize;

		contour(pg->X, pg->Z, pg->Psi, 0, nmax, 0, nmax, Psi, CONTOUR_ONLY_CLOSED,
			CONTOUR_NO_MIDPOINT, CountContourStep);

		npts = gCount;

		*X = gXvec = dvector(0,npts);
		*Z = gZvec = dvector(0,npts);

		contour(pg->X, pg->Z, pg->Psi, 0, nmax, 0, nmax, Psi, CONTOUR_ONLY_CLOSED,
			CONTOUR_NO_MIDPOINT, RecordContourStep);

		return npts;
}
