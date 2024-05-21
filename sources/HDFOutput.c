/*
** Dipole 0.9
**
** HDFOutput.c
**
** Routines to write binary HDF files for graphical
** data analysis.
**
** File:		HDFOutput.c
** Date:		April 2, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdlib.h>
#include "VAX_Alloc.h"
#include <stdio.h>
#include <string.h>
#include "nrutil.h"
#include <unistd.h>
#include "multitask.h"

/* The following machine definitions are needed for the HDF library */
#ifdef THINK_C
#define 	MAC 	1
#endif
#ifdef __SC__
#define 	MAC 	1
#endif
#ifdef VAXC
#define 	VMS 	1
#endif

#ifdef __cplusplus
extern "C"
{
#endif

/* #include "hdf.h" */
#include "mfhdf.h"

#ifdef __cplusplus
}
#endif

typedef float float32;

#include "psigrid.h"
#include "plasma.h"
#include "HDFOutput.h"

#define SDCHK(x)	if ((x) == FAIL) { \
						 fprintf(stdout,"when calling %s\n",#x); \
						 HEprint(stdout,0); \
                         nrerror("ERROR: SD-HDF error in HDFOutput.");  \
                     }

#define COMPRESS_TYPE 	COMP_CODE_DEFLATE
#define COMPRESS_LEVEL 	6
#define COMPRESS(x) 	SDCHK(SDsetcompress (x, COMP_CODE_DEFLATE, &c_info))

static comp_info c_info = {.deflate.level = 6};	/* Compression structure */


#define PI	3.14159265358979323
#define MU0	1.25663706e-06
#define TWOPI	6.283185307

float32      *AryToFloat32(double **ad, int nmax, double multiplier);
float32      *VecToFloat32(double *vd, int nmax, double multiplier);

#define	TO_FLOAT32(d,a,npts) for(i=0;i<npts;i++) a[i] = (float32)d[i]

/*
**	AryToFloat32
**
**
*/
float32      *AryToFloat32(double **ad, int nmax, double multiplier)
{
	int           ix, iz;
	float32      *a, *a0;

	a0 = a = (float32 *) malloc((unsigned) ((nmax + 1) * (nmax + 1)) * sizeof(float32));
	if (!a)
		nrerror("ERROR: Allocation failure in AryToFloat32.");

	for (iz = 0; iz <= nmax; iz++)	/* slowest */
		for (ix = 0; ix <= nmax; ix++)	/* fastest */
			*a++ = (float32) (ad[ix][iz] * multiplier);

	return a0;
}

/*
**	VecToFloat32
**
**
*/
float32      *VecToFloat32(double *vd, int nmax, double multiplier)
{
	int           ix;
	float32      *v, *v0;

	v0 = v = (float32 *) malloc((unsigned) (nmax + 1) * sizeof(float32));
	if (!v)
		nrerror("ERROR: Allocation failure in VecToFloat32.");

	for (ix = 0; ix <= nmax; ix++)
		*v++ = (float32)( vd[ix] * multiplier );

	return v0;
}

/*
**	HDFPsiGrid
**
**
*/

#define DEBUG_HDF 0

#if DEBUG_HDF
char folder[FILENAME_MAX];
char curFolder[FILENAME_MAX];
char oldFolder[FILENAME_MAX];

#define TELL_FOLDER(x) getcwd( folder, FILENAME_MAX ); \
        printf("%s, the current Folder is: %s\n", x, folder)
#else
#define TELL_FOLDER(x) {}
#endif

void          HDFPsiGrid(PSIGRID * pg, char *Oname)
{
	int           nmax;
	int32         dims[2], starts[2];
	float32      *a, *v, maxv, minv;
	char          fname[FILENAME_MAX] = "";
//	char		curFolder[FILENAME_MAX];
	int32		sd_id, status, sds_id, dim_id;

	nmax = pg->Nsize;

	dims[0] = nmax + 1;
	dims[1] = nmax + 1;
	starts[0]=0;
	starts[1]=0;

	/* this fixes a stupid bug in HDF 4.1r3 */
//	getcwd(curFolder,FILENAME_MAX);

//	strcpy(fname,curFolder);
	strncpy(fname,Oname,18);
	strcat(fname,".hdf");

	TELL_FOLDER("Before HDF start");

	/* S E T U P    H D F */
	MULTI;

	SDCHK(sd_id = SDstart(fname,DFACC_CREATE));  /* overwrite old file */

	/* Do Current first and dimensions */

	TELL_FOLDER("After HDF start");


	MULTI;

    TELL_FOLDER("Write 0");

	sds_id = SDcreate(sd_id,CUR_NAME,DFNT_FLOAT32,2,dims);
	SDCHK(sds_id);

    TELL_FOLDER("Write 1");

	COMPRESS(sds_id);
    TELL_FOLDER("Compress 1");

	status = SDsetdatastrs(sds_id,CUR_NAME,"A/m^2","E10.3", "cartesian");
	SDCHK(status);
    TELL_FOLDER("Write 2");

	a = AryToFloat32(pg->Current, nmax, 1.0 / MU0);
	status = SDwritedata(sds_id,starts,NULL,dims,(VOIDP)a);
	free(a);
	SDCHK(status);
    TELL_FOLDER("Write 3");

	MULTI;
   TELL_FOLDER("Write 3a");

	/* Define X dimension */
	dim_id = SDgetdimid(sds_id, 1);
	SDCHK(dim_id);
   TELL_FOLDER("Write 3b");

	status = SDsetdimname(dim_id,DIMX_NAME);
	SDCHK(status);
   TELL_FOLDER("Write 3c");

	status = SDsetdimval_comp(dim_id,SD_DIMVAL_BW_COMP);
	SDCHK(status);
    TELL_FOLDER("Write 4");

	v = VecToFloat32(pg->X, nmax, 1.0);
	status = SDsetdimscale(dim_id, dims[0], DFNT_FLOAT32, (VOIDP) v);
	free(v);
	SDCHK(status);

	status = SDsetdimstrs(dim_id,DIMX_NAME,"m","F7.4");
	SDCHK(status);

	MULTI;
    TELL_FOLDER("Write 5");

	/* now define Z dimension */
	dim_id = SDgetdimid(sds_id, 0);
	SDCHK(dim_id);

	status = SDsetdimname(dim_id,DIMZ_NAME);
	SDCHK(status);

	status = SDsetdimval_comp(dim_id,SD_DIMVAL_BW_COMP);
	SDCHK(status);

	v = VecToFloat32(pg->Z, nmax, 1.0);
	status = SDsetdimscale(dim_id, dims[1],DFNT_FLOAT32, (VOIDP) v);
	free(v);
	SDCHK(status);

	status = SDsetdimstrs(dim_id,DIMZ_NAME,"m","F7.4");
	SDCHK(status);

	status = SDendaccess (sds_id);
	SDCHK(status);
    TELL_FOLDER("Write 6");

	/* PSI */
	MULTI;

	sds_id = SDcreate(sd_id,PSI_NAME,DFNT_FLOAT32,2,dims);
	SDCHK(sds_id);

	COMPRESS(sds_id);

	status = SDsetdatastrs(sds_id,PSI_NAME,"weber","E10.3", "cartesian");
	SDCHK(status);

	a = AryToFloat32(pg->Psi, nmax, 1.0);
	status = SDwritedata(sds_id,starts,NULL,dims,(VOIDP)a);
	free(a);
	SDCHK(status);
    TELL_FOLDER("Write 7");

	/* also store the range from PsiAxis to PsiLim (PsiSep) */
	MULTI;

	maxv = (float32) pg->PsiLim;
	minv = (float32) pg->PsiAxis;

	status = SDsetrange(sds_id,(VOIDP) &maxv, (VOIDP) &minv);

	SDCHK(dim_id = SDgetdimid(sds_id, 1));
	SDCHK(status = SDsetdimname(dim_id,DIMX_NAME));
	SDCHK(dim_id = SDgetdimid(sds_id, 0));
	SDCHK(status = SDsetdimname(dim_id,DIMZ_NAME));
	SDCHK(status = SDendaccess (sds_id));
    TELL_FOLDER("Write 8");

	/* Residuals */
	MULTI;

	sds_id = SDcreate(sd_id,RES_NAME,DFNT_FLOAT32,2,dims);
	SDCHK(sds_id);

	COMPRESS(sds_id);

	status = SDsetdatastrs(sds_id,RES_NAME,"A/m^2","E10.3", "cartesian");
	SDCHK(status);

	a = AryToFloat32(pg->Residual, nmax, 1.0 / MU0);
	status = SDwritedata(sds_id,starts,NULL,dims,(VOIDP)a);
	free(a);
	SDCHK(status);
    TELL_FOLDER("Write 9");

	SDCHK(dim_id = SDgetdimid(sds_id, 1));
	SDCHK(status = SDsetdimname(dim_id,DIMX_NAME));
	SDCHK(dim_id = SDgetdimid(sds_id, 0));
	SDCHK(status = SDsetdimname(dim_id,DIMZ_NAME));
	SDCHK(status = SDendaccess (sds_id));

	TELL_FOLDER("After much HDF psigrid writing");
	status = SDend (sd_id);

	TELL_FOLDER("After closing file");

	MULTI;

		/* this fixes a stupid bug in HDF 4.1r3 */
//	chdir(curFolder);

	TELL_FOLDER("Just reset the folder");

}


/*
**	HDFPlasma
**
**
*/
void          HDFPlasma(PLASMA * pl, PSIGRID * pg, char *Oname)
{
	int32         dims[2], starts[2], nmax;
	float32      *a;
	char          fname[FILENAME_MAX] = "";
//	char		curFolder[FILENAME_MAX];

	int32		sd_id, status, sds_id, dim_id;
	double		**ad;
	int			ix,iz;

	nmax = pg->Nsize;

	dims[0] = nmax + 1;
	dims[1] = nmax + 1;
	starts[0]=0;
	starts[1]=0;

	/* this fixes a stupid bug in HDF 4.1r3 */
//	getcwd(curFolder,FILENAME_MAX);

//	strcpy(fname,curFolder);
//	strncat(fname,Oname,18);
	strncpy(fname,Oname,18);
	strcat(fname,".hdf");

	ad = dmatrix(0, nmax, 0, nmax);	/* for temporary storage */

	TELL_FOLDER("Just about to open");

	/* S E T U P    H D F */
    MULTI;
	SDCHK(sd_id = SDstart(fname,DFACC_WRITE));  /* concat into existing file */

	TELL_FOLDER("Opened 2nd time");

	/* B 2 */
	MULTI;
	sds_id = SDcreate(sd_id,MODB_NAME,DFNT_FLOAT32,2,dims);
	SDCHK(sds_id);

	COMPRESS(sds_id);

	status = SDsetdatastrs(sds_id,MODB_NAME,"T2", "F7.4", "cartesian");
	SDCHK(status);

	a = AryToFloat32(pl->B2, nmax, 1.0);
	status = SDwritedata(sds_id,starts,NULL,dims,(VOIDP)a);
	free(a);
	SDCHK(status);

	SDCHK(dim_id = SDgetdimid(sds_id, 1));
	SDCHK(status = SDsetdimname(dim_id,DIMX_NAME));
	SDCHK(dim_id = SDgetdimid(sds_id, 0));
	SDCHK(status = SDsetdimname(dim_id,DIMZ_NAME));
	SDCHK(status = SDendaccess (sds_id));

		/* Bp_x */
		MULTI;
	sds_id = SDcreate(sd_id,BpX_NAME,DFNT_FLOAT32,2,dims);
	SDCHK(sds_id);

	COMPRESS(sds_id);

	status = SDsetdatastrs(sds_id,BpX_NAME, "T", "F7.4", "cartesian");
	SDCHK(status);

	for (ix = 0; ix <= nmax; ix++)
		for (iz = 0; iz <= nmax; iz++)
			ad[ix][iz] = pl->GradPsiZ[ix][iz] / TWOPI / pg->X[ix];
	a = AryToFloat32(ad, nmax, 1.0);

	status = SDwritedata(sds_id,starts,NULL,dims,(VOIDP)a);
	free(a);
	SDCHK(status);

	SDCHK(dim_id = SDgetdimid(sds_id, 1));
	SDCHK(status = SDsetdimname(dim_id,DIMX_NAME));
	SDCHK(dim_id = SDgetdimid(sds_id, 0));
	SDCHK(status = SDsetdimname(dim_id,DIMZ_NAME));
	SDCHK(status = SDendaccess (sds_id));

	/* Bp_z */
	MULTI;
	sds_id = SDcreate(sd_id,BpZ_NAME,DFNT_FLOAT32,2,dims);
	SDCHK(sds_id);

	COMPRESS(sds_id);

	status = SDsetdatastrs(sds_id,BpZ_NAME, "T", "F7.4", "cartesian");
	SDCHK(status);

	for (ix = 0; ix <= nmax; ix++)
		for (iz = 0; iz <= nmax; iz++)
			ad[ix][iz] = -pl->GradPsiX[ix][iz] / TWOPI / pg->X[ix];
	a = AryToFloat32(ad, nmax, 1.0);
	status = SDwritedata(sds_id,starts,NULL,dims,(VOIDP)a);
	free(a);
	SDCHK(status);

	SDCHK(dim_id = SDgetdimid(sds_id, 1));
	SDCHK(status = SDsetdimname(dim_id,DIMX_NAME));
	SDCHK(dim_id = SDgetdimid(sds_id, 0));
	SDCHK(status = SDsetdimname(dim_id,DIMZ_NAME));
	SDCHK(status = SDendaccess (sds_id));

	/* G */
	MULTI;
	sds_id = SDcreate(sd_id,TFLUX_NAME,DFNT_FLOAT32,2,dims);
	SDCHK(sds_id);

	COMPRESS(sds_id);

	status = SDsetdatastrs(sds_id,TFLUX_NAME," ", "F7.4", "cartesian");
	SDCHK(status);

	a = AryToFloat32(pl->G, nmax, 1.0);
	status = SDwritedata(sds_id,starts,NULL,dims,(VOIDP)a);
	free(a);
	SDCHK(status);

	SDCHK(dim_id = SDgetdimid(sds_id, 1));
	SDCHK(status = SDsetdimname(dim_id,DIMX_NAME));
	SDCHK(dim_id = SDgetdimid(sds_id, 0));
	SDCHK(status = SDsetdimname(dim_id,DIMZ_NAME));
	SDCHK(status = SDendaccess (sds_id));

	/* P R E S S U R E */
	MULTI;
	sds_id = SDcreate(sd_id,PRESS_NAME,DFNT_FLOAT32,2,dims);
	SDCHK(sds_id);

	COMPRESS(sds_id);

	status = SDsetdatastrs(sds_id,PRESS_NAME, "pascal", "E10.3", "cartesian");
	SDCHK(status);

	switch (pl->ModelType) {
      default :
			a = AryToFloat32(pl->Piso, nmax, 1.0 / MU0);
			status = SDwritedata(sds_id,starts,NULL,dims,(VOIDP)a);
			free(a);
			SDCHK(status);
		  	break;
	  case Plasma_IsoFlow:

		  break;
	  case Plasma_AnisoNoFlow:

		  break;
	  case Plasma_AnisoFlow:

		  break;
	}

	SDCHK(dim_id = SDgetdimid(sds_id, 1));
	SDCHK(status = SDsetdimname(dim_id,DIMX_NAME));
	SDCHK(dim_id = SDgetdimid(sds_id, 0));
	SDCHK(status = SDsetdimname(dim_id,DIMZ_NAME));
	SDCHK(status = SDendaccess (sds_id));

		/* B E T A */
		MULTI;
	SDCHK(sds_id = SDcreate(sd_id,BETA_NAME,DFNT_FLOAT32,2,dims));

	COMPRESS(sds_id);

	SDCHK(status = SDsetdatastrs(sds_id,BETA_NAME, "", "E10.3", "cartesian"));

	switch (pl->ModelType) {
      default :
      	    for (ix = 1; ix <= nmax; ix++)
        	    for (iz = 1; iz <= nmax; iz++)
        	       ad[ix][iz] = 2*pl->Piso[ix][iz]/pl->B2[ix][iz];
			a = AryToFloat32(ad, nmax, 1.0);
			status = SDwritedata(sds_id,starts,NULL,dims,(VOIDP)a);
			free(a);
			SDCHK(status);
		  	break;
	  case Plasma_IsoFlow:

		  break;
	  case Plasma_AnisoNoFlow:

		  break;
	  case Plasma_AnisoFlow:

		  break;
	}

	SDCHK(dim_id = SDgetdimid(sds_id, 1));
	SDCHK(status = SDsetdimname(dim_id,DIMX_NAME));
	SDCHK(dim_id = SDgetdimid(sds_id, 0));
	SDCHK(status = SDsetdimname(dim_id,DIMZ_NAME));
	SDCHK(status = SDendaccess (sds_id));

	free_dmatrix(ad, 0, nmax, 0, nmax);

	SDCHK(status = SDend (sd_id));

	TELL_FOLDER("After closing file");

	MULTI;

		/* this fixes a stupid bug in HDF 4.1r3 */
//	chdir(curFolder);

	TELL_FOLDER("Just reset the folder");

}

/*
**	HDFFluxFuncs
**
**
*/
void	HDFFluxFuncs(char *Oname, int npts, double *PsiX,
					double *Psi, double *P, double *G, double *Pp, double *G2p,
					double *q, double *dVdpsi, double *Vol, double *Shear,
					double *Well, double *Jave, double *B2ave, double *Beta)
{
	int32         dim[1], start[1];
	float32      *a;

	int32		sd_id, sds_id, dim_id;
	char          fname[FILENAME_MAX] = "";
//	char		curFolder[FILENAME_MAX];

	int			i;

	a = (float32 *) malloc (npts*sizeof(float32));

	dim[0] = npts;
	start[0]=0;


	/* this fixes a stupid bug in HDF 4.1r3 */
//	getcwd(curFolder,FILENAME_MAX);

//	strcpy(fname,curFolder);
//	strncat(fname,Oname,18);
	strncpy(fname,Oname,18);
	strcat(fname,".hdf");


	TELL_FOLDER("Just about to open");

	/* S E T U P    H D F */
    MULTI;
	SDCHK(sd_id = SDstart(fname,DFACC_WRITE));  /* concat into existing file */

	TELL_FOLDER("Opened 3rd time");

	/* Psi */
	MULTI;
	SDCHK( sds_id = SDcreate(sd_id,PSI_1D,DFNT_FLOAT32,1,dim) );

	/* Define PsiX dimension */
	MULTI;
	SDCHK( dim_id = SDgetdimid(sds_id, 0) );
	SDCHK( SDsetdimname(dim_id,PSIX_NAME) );

	SDCHK( SDsetdimval_comp(dim_id,SD_DIMVAL_BW_COMP)    );
	TO_FLOAT32(PsiX,a,npts);

	SDCHK( SDsetdimscale(dim_id, npts, DFNT_FLOAT32, (VOIDP) a));
	SDCHK( SDsetdimstrs(dim_id,"Psi","normalized","F7.4")     );

	SDCHK( SDsetdatastrs(sds_id,"Psi","weber","E10.3", "poloidal flux") );
	TO_FLOAT32(Psi,a,npts);
	SDCHK( SDwritedata(sds_id,start,NULL,dim,(VOIDP) a) );
	SDCHK( SDendaccess(sds_id) );

	/* Pressure */
	MULTI;
	SDCHK( sds_id = SDcreate(sd_id,PRESS_1D,DFNT_FLOAT32,1,dim) );
	SDCHK( dim_id = SDgetdimid(sds_id, 0) );
	SDCHK( SDsetdimname(dim_id,PSIX_NAME) ); /* use the same psix */
	SDCHK( SDsetdatastrs(sds_id,"Pressure","pascals","E10.3", "poloidal flux") );
	TO_FLOAT32(P,a,npts);
	SDCHK( SDwritedata(sds_id,start,NULL,dim,(VOIDP) a) );
	SDCHK( SDendaccess(sds_id) );

	/* G */
	MULTI;
	SDCHK( sds_id = SDcreate(sd_id,G_1D,DFNT_FLOAT32,1,dim) );
	SDCHK( dim_id = SDgetdimid(sds_id, 0) );
	SDCHK( SDsetdimname(dim_id,PSIX_NAME) ); /* use the same psix */
	SDCHK( SDsetdatastrs(sds_id,"ToroidalFlux","T-m^2","E10.3", "poloidal flux") );
	TO_FLOAT32(G,a,npts);
	SDCHK( SDwritedata(sds_id,start,NULL,dim,(VOIDP) a) );
	SDCHK( SDendaccess(sds_id) );

	/* dP_dPsi */
	MULTI;
	SDCHK( sds_id = SDcreate(sd_id,PP_1D,DFNT_FLOAT32,1,dim) );
	SDCHK( dim_id = SDgetdimid(sds_id, 0) );
	SDCHK( SDsetdimname(dim_id,PSIX_NAME) ); /* use the same psix */
	SDCHK( SDsetdatastrs(sds_id,"dP/dPsi","pascals/weber","E10.3", "poloidal flux") );
	TO_FLOAT32(Pp,a,npts);
	SDCHK( SDwritedata(sds_id,start,NULL,dim,(VOIDP) a) );
	SDCHK( SDendaccess(sds_id) );

	/* g2p */
	MULTI;
	SDCHK( sds_id = SDcreate(sd_id,G2P_1D,DFNT_FLOAT32,1,dim) );
	SDCHK( dim_id = SDgetdimid(sds_id, 0) );
	SDCHK( SDsetdimname(dim_id,PSIX_NAME) ); /* use the same psix */
	SDCHK( SDsetdatastrs(sds_id,"dG^2/dPsi","T-m^2","E10.3", "poloidal flux") );
	TO_FLOAT32(G2p,a,npts);
	SDCHK( SDwritedata(sds_id,start,NULL,dim,(VOIDP) a) );
	SDCHK( SDendaccess(sds_id) );

	/* q */
	MULTI;
	SDCHK( sds_id = SDcreate(sd_id,Q_1D,DFNT_FLOAT32,1,dim) );
	SDCHK( dim_id = SDgetdimid(sds_id, 0) );
	SDCHK( SDsetdimname(dim_id,PSIX_NAME) ); /* use the same psix */
	SDCHK( SDsetdatastrs(sds_id,"q","","E10.3", "poloidal flux") );
	TO_FLOAT32(q,a,npts);
	SDCHK( SDwritedata(sds_id,start,NULL,dim,(VOIDP) a) );
	SDCHK( SDendaccess(sds_id) );

	/* dV/dPsi */
	MULTI;
	SDCHK( sds_id = SDcreate(sd_id,V_1D,DFNT_FLOAT32,1,dim) );
	SDCHK( dim_id = SDgetdimid(sds_id, 0) );
	SDCHK( SDsetdimname(dim_id,PSIX_NAME) ); /* use the same psix */
	SDCHK( SDsetdatastrs(sds_id,"V","m / T","E10.3", "poloidal flux") );
	TO_FLOAT32(dVdpsi,a,npts);
	SDCHK( SDwritedata(sds_id,start,NULL,dim,(VOIDP) a) );
	SDCHK( SDendaccess(sds_id) );

	/* Vol */
	MULTI;
	SDCHK( sds_id = SDcreate(sd_id,VOL_1D,DFNT_FLOAT32,1,dim) );
	SDCHK( dim_id = SDgetdimid(sds_id, 0) );
	SDCHK( SDsetdimname(dim_id,PSIX_NAME) ); /* use the same psix */
	SDCHK( SDsetdatastrs(sds_id,"Integrated Volume","m^3","E10.3", "poloidal flux") );
	TO_FLOAT32(Vol,a,npts);
	SDCHK( SDwritedata(sds_id,start,NULL,dim,(VOIDP) a) );
	SDCHK( SDendaccess(sds_id) );

	/* Shear */
	MULTI;
	SDCHK( sds_id = SDcreate(sd_id,SHEAR_1D,DFNT_FLOAT32,1,dim) );
	SDCHK( dim_id = SDgetdimid(sds_id, 0) );
	SDCHK( SDsetdimname(dim_id,PSIX_NAME) ); /* use the same psix */
	SDCHK( SDsetdatastrs(sds_id,"Shear","","E10.3", "poloidal flux") );
	TO_FLOAT32(Shear,a,npts);
	SDCHK( SDwritedata(sds_id,start,NULL,dim,(VOIDP) a) );
	SDCHK( SDendaccess(sds_id) );

	/* Well */
	MULTI;
	SDCHK( sds_id = SDcreate(sd_id,WELL_1D,DFNT_FLOAT32,1,dim) );
	SDCHK( dim_id = SDgetdimid(sds_id, 0) );
	SDCHK( SDsetdimname(dim_id,PSIX_NAME) ); /* use the same psix */
	SDCHK( SDsetdatastrs(sds_id,"Well","","E10.3", "poloidal flux") );
	TO_FLOAT32(Well,a,npts);
	SDCHK( SDwritedata(sds_id,start,NULL,dim,(VOIDP) a) );
	SDCHK( SDendaccess(sds_id) );

	/* J ave */
	MULTI;
	SDCHK( sds_id = SDcreate(sd_id,J_1D,DFNT_FLOAT32,1,dim) );
	SDCHK( dim_id = SDgetdimid(sds_id, 0) );
	SDCHK( SDsetdimname(dim_id,PSIX_NAME) ); /* use the same psix */
	SDCHK( SDsetdatastrs(sds_id,"<J>","A/m^2","E10.3", "poloidal flux") );
	TO_FLOAT32(Jave,a,npts);
	SDCHK( SDwritedata(sds_id,start,NULL,dim,(VOIDP) a) );
	SDCHK( SDendaccess(sds_id) );

	/* B^2 ave */
	MULTI;
	SDCHK( sds_id = SDcreate(sd_id,B2_1D,DFNT_FLOAT32,1,dim) );
	SDCHK( dim_id = SDgetdimid(sds_id, 0) );
	SDCHK( SDsetdimname(dim_id,PSIX_NAME) ); /* use the same psix */
	SDCHK( SDsetdatastrs(sds_id,"<B^2>","T^2","E10.3", "poloidal flux") );
	TO_FLOAT32(B2ave,a,npts);
	SDCHK( SDwritedata(sds_id,start,NULL,dim,(VOIDP) a) );
	SDCHK( SDendaccess(sds_id) );

	/* Beta ave */
	MULTI;
	SDCHK( sds_id = SDcreate(sd_id,BETA_1D,DFNT_FLOAT32,1,dim) );
	SDCHK( dim_id = SDgetdimid(sds_id, 0) );
	SDCHK( SDsetdimname(dim_id,PSIX_NAME) ); /* use the same psix */
	SDCHK( SDsetdatastrs(sds_id,"<Beta>","","E10.3", "") );
	TO_FLOAT32(Beta,a,npts);
	SDCHK( SDwritedata(sds_id,start,NULL,dim,(VOIDP) a) );
	SDCHK( SDendaccess(sds_id) );


	TELL_FOLDER("After much HDF flux profile writing");
	SDCHK(SDend (sd_id));

	TELL_FOLDER("After closing file");

	free(a);


		/* this fixes a stupid bug in HDF 4.1r3 */
//	chdir(curFolder);

	TELL_FOLDER("After changing folder");


}


/*
**	HDFBoundary
**
*/
void    HDFBoundary(const char *Oname, const char *Vname, double psiLabel,
                     double *X, double *Z, int len)
{
	int32         	dims[2], starts[2];
	float32      	*a, psi;
	char          	fname[FILENAME_MAX] = "";
//	char		curFolder[FILENAME_MAX];

	int32		sd_id, status, sds_id;
	int		i;

	dims[0] = 	len;
	dims[1] = 	2;
	starts[0]=	0;
	starts[1]=	0;

	/* this fixes a stupid bug in HDF 4.1r3 */
//	getcwd(curFolder,FILENAME_MAX);

//	strcpy(fname,curFolder);
//	strncat(fname,Oname,18);
	strncpy(fname,Oname,18);
	strcat(fname,".hdf");

//	ad = dmatrix(0, nmax, 0, nmax);	/* for temporary storage */

	TELL_FOLDER("Just about to open");

	/* S E T U P    H D F */
        MULTI;
	SDCHK(sd_id = SDstart(fname,DFACC_WRITE));  /* concat into existing file */

	TELL_FOLDER("Opened 2nd time");

	/* create the boundary... */
	MULTI;
	sds_id = SDcreate(sd_id,Vname,DFNT_FLOAT32,2,dims);
	SDCHK(sds_id);

	COMPRESS(sds_id);

        psi = (float32) psiLabel;
        SDCHK(SDsetattr(sds_id,"PsiNormalized",DFNT_FLOAT32,1,(VOIDP)&psi));

	status = SDsetdatastrs(sds_id,Vname,"m", "F7.4", "cartesian");
	SDCHK(status);

        a = (float32 *) malloc((len+1)*2*sizeof(float32));
        if (!a)
            nrerror("Failure to allocate temporary storarge in HDFBoundary.");
        for (i=0; i<len; i++) {
            a[2*i]   = (float32) X[i];
            a[2*i+1] = (float32) Z[i];
        }

        status = SDwritedata(sds_id,starts,NULL,dims,(VOIDP)a);
	free(a);
	SDCHK(status);

	SDCHK(status = SDendaccess (sds_id));

	SDCHK(status = SDend (sd_id));

	TELL_FOLDER("After closing file");

	MULTI;

		/* this fixes a stupid bug in HDF 4.1r3 */
//	chdir(curFolder);

	TELL_FOLDER("Just reset the folder");
}
