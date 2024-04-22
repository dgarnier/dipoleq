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

#ifdef __cplusplus
extern "C"
{
#endif

#include "hdf5.h"
#include "hdf5_hl.h"

#ifdef __cplusplus
}
#endif

typedef float float32;

#include "psigrid.h"
#include "plasma.h"
#include "HDFOutput.h"

#define H5CHK(x)	if ((x) < 0) do { \
						H5Eprint2(H5E_DEFAULT, stderr); \
						nrerror("ERROR: HDF5 error in HDFOutput."); \
					} while (0)


#define PI	3.14159265358979323
#define MU0	1.25663706e-06
#define TWOPI	6.283185307
 
/*
**	ScaleArray
**
**
*/
static void ScaleArray(double *a, double **ad, int nmax, double multiplier)
{
	int ix, iz;

	for (iz = 0; iz <= nmax; iz++)	   /* slowest */
		for (ix = 0; ix <= nmax; ix++) /* fastest */
			*a++ = ad[ix][iz] * multiplier;
}

/*
**	ScaleVec
**
**
*/
static void ScaleVec(double *v, double *vd, int nmax, double multiplier)
{
	int ix;

	for (ix = 0; ix <= nmax; ix++)
		*v++ = vd[ix] * multiplier;
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


static void HDFWrite1D(double * data, const char *name, const char *units, hid_t loc, const hsize_t *dims, hid_t dsc1)
{
	hid_t ds_id;

	H5CHK(H5LTmake_dataset_double(loc, name, 1, dims, data));
	H5CHK(H5LTset_attribute_string(loc, name, "UNITS", units));
	H5CHK(ds_id = H5Dopen2(loc, name, H5P_DEFAULT));
	H5CHK(H5DSattach_scale(ds_id, dsc1, 0));
	H5CHK(H5Dclose(ds_id));
}

static void HDFWrite2D(double *data, const char *name, const char *units, hid_t loc, hid_t dsp2d, hid_t dsc1, hid_t dsc2) {
	hid_t ds_id, attr_id;

	H5CHK(ds_id = H5Dcreate2(loc, name, H5T_NATIVE_DOUBLE, dsp2d, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
	H5CHK(H5Dwrite(ds_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data));
	H5CHK(H5DSattach_scale(ds_id, dsc1, 0));
	H5CHK(H5DSattach_scale(ds_id, dsc2, 1));
	H5CHK(H5Dclose(ds_id));

	// use high level interface to set attributes
	H5CHK(H5LTset_attribute_string(loc, name, "UNITS", units));
	// H5CHK(H5LTset_attribute_string(file, Vname, "DIMENSION", "cartesian"));
	// H5CHK(H5LTset_attribute_string(file, Vname, "FORMAT", "F7.4"));
}

void          HDFPsiGrid(PSIGRID * pg, char *Oname)
{
	char   		fname[FILENAME_MAX] = "";
	hid_t  		file, dsp1d, dsp2d, dsc1, dsc2;
	hsize_t		dims[2], nmax;
	double *	a;

	nmax = pg->Nsize;
	dims[0] = dims[1] = nmax + 1;
	a = (double *) malloc(((nmax + 1) * (nmax + 1)) * sizeof(double));

	strncpy(fname, Oname, 18);
	strcat(fname, ".h5");

	/* S E T U P    H D F */

	/* Create a new file using default properties. Overwrite the old file. */
	H5CHK(file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
	
	/* Do dimensions first, aka Dimension Scales */

	/* dataspaces */
	H5CHK(dsp1d = H5Screate_simple(1, dims, NULL));
	H5CHK(dsp2d = H5Screate_simple(2, dims, NULL));

	/* dimension datasets -> datascale datasets*/
	H5CHK(dsc1 = H5Dcreate2(file, DIMX_NAME, H5T_NATIVE_DOUBLE, dsp1d, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
	H5CHK(dsc2 = H5Dcreate2(file, DIMZ_NAME, H5T_NATIVE_DOUBLE, dsp1d, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
	H5CHK(H5Dwrite(dsc1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pg->X));
	H5CHK(H5Dwrite(dsc2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pg->Z));
	H5CHK(H5DSset_scale(dsc1, DIMX_NAME));
	H5CHK(H5DSset_scale(dsc2, DIMZ_NAME));
	H5CHK(H5LTset_attribute_string(dsc1, ".", "UNITS", "m"));
	H5CHK(H5LTset_attribute_string(dsc2, ".", "UNITS", "m"));

	/* Do Current */
	ScaleArray(a, pg->Current, nmax, 1.0 / MU0);
	HDFWrite2D(a, CUR_NAME, "A/m^2", file, dsp2d, dsc1, dsc2);

	/* Do Psi */
	ScaleArray(a, pg->Psi, nmax, 1.0);
	HDFWrite2D(a, PSI_NAME, "Wb", file, dsp2d, dsc1, dsc2);
	H5LTset_attribute_double(file, PSI_NAME, "PsiAxis", &pg->PsiAxis, 1);
	H5LTset_attribute_double(file, PSI_NAME, "PsiLim", &pg->PsiLim, 1);
	H5LTset_attribute_double(file, PSI_NAME, "DelPsi", &pg->DelPsi, 1);

	/* Do Residuals */
	ScaleArray(a, pg->Residual, nmax, 1.0 / MU0);
	HDFWrite2D(a, RES_NAME, "A/m^2", file, dsp2d, dsc1, dsc2);
	

	/* Close the file */
	H5CHK(H5Dclose(dsc1));
	H5CHK(H5Dclose(dsc2));
	H5CHK(H5Sclose(dsp1d));
	H5CHK(H5Sclose(dsp2d));
	H5CHK(H5Fclose(file));

	free(a);
}


/*
**	HDFPlasma
**
**
*/
void          HDFPlasma(PLASMA * pl, PSIGRID * pg, char *Oname)
{
	char fname[FILENAME_MAX] = "";
	hid_t file, dsp1d, dsp2d, dsc1, dsc2;
	hsize_t dims[2], nmax;
	double *a, *ap;

	nmax = pg->Nsize;
	dims[0] = dims[1] = nmax + 1;

	strncpy(fname, Oname, 18);
	strcat(fname, ".h5");
	a = (double *) malloc(((nmax + 1) * (nmax + 1)) * sizeof(double));

	TELL_FOLDER("Just about to open");

	/* S E T U P    H D F */
    MULTI;
	H5CHK(file = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT));

	/* dataspaces */
	H5CHK(dsp1d = H5Screate_simple(1, dims, NULL));
	H5CHK(dsp2d = H5Screate_simple(2, dims, NULL));

	/* read dimension scales */
	H5CHK(dsc1 = H5Dopen2(file, DIMX_NAME, H5P_DEFAULT));
	H5CHK(dsc2 = H5Dopen2(file, DIMZ_NAME, H5P_DEFAULT));


	TELL_FOLDER("Opened 2nd time");	
	
	/* B 2 */
	MULTI;
	ScaleArray(a, pl->B2, nmax, 1.0);
	HDFWrite2D(a, MODB_NAME, "T^2", file, dsp2d, dsc1, dsc2);
	
	/* Bp_x */
	MULTI;
	ap = a;
	for (int ix = 0; ix <= nmax; ix++)
		for (int iz = 0; iz <= nmax; iz++)
			*ap++ = pl->GradPsiZ[ix][iz] / TWOPI / pg->X[ix];

	HDFWrite2D(a, BpX_NAME, "T", file, dsp2d, dsc1, dsc2);
	
	/* Bp_z */
	MULTI;
	ap = a;
	for (int ix = 0; ix <= nmax; ix++)
		for (int iz = 0; iz <= nmax; iz++)
			*ap++ = -pl->GradPsiX[ix][iz] / TWOPI / pg->X[ix];

	HDFWrite2D(a, BpZ_NAME, "T", file, dsp2d, dsc1, dsc2);
	
	/* G */
	MULTI;
	ScaleArray(a, pl->G, nmax, 1.0);
	HDFWrite2D(a, TFLUX_NAME, "Wb", file, dsp2d, dsc1, dsc2);
	
	/* P R E S S U R E */
	MULTI;
	switch (pl->ModelType) {
      default :
			ScaleArray(a, pl->Piso, nmax, 1.0 / MU0);
			HDFWrite2D(a, PRESS_NAME, "Pa", file, dsp2d, dsc1, dsc2);
		  	break;
	  case Plasma_IsoFlow:

		  break;
	  case Plasma_AnisoNoFlow:

		  break;
	  case Plasma_AnisoFlow:

		  break;
	}
	
	/* B E T A */
	MULTI;
	switch (pl->ModelType) {
	  default :
	  		ap = a;
	  	    for (int ix = 1; ix <= nmax; ix++)
			    for (int iz = 1; iz <= nmax; iz++)
			       *ap++ = 2*pl->Piso[ix][iz]/pl->B2[ix][iz];
			HDFWrite2D(a, BETA_NAME, "", file, dsp2d, dsc1, dsc2);
		  	break;
	  case Plasma_IsoFlow:

		  break;
	  case Plasma_AnisoNoFlow:

		  break;
	  case Plasma_AnisoFlow:

		  break;
	}

	/* Close the file */
	H5CHK(H5Dclose(dsc1));
	H5CHK(H5Dclose(dsc2));
	H5CHK(H5Sclose(dsp1d));
	H5CHK(H5Sclose(dsp2d));
	H5CHK(H5Fclose(file));

	free(a);
	
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
	char fname[FILENAME_MAX] = "";
	hid_t file, dsp1d, dsrho;
	hsize_t dims[1];

	dims[0] = npts;

	strncpy(fname, Oname, 18);
	strcat(fname, ".h5");

	TELL_FOLDER("Just about to open");

	/* S E T U P    H D F */
	MULTI;
	H5CHK(file = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT));

	/* dataspaces */
	// H5CHK(dsp1d = H5Screate_simple(1, dims, NULL));

	/* psi_x dimension scales */

	H5CHK(H5LTmake_dataset_double(file, PSIX_NAME, 1, dims, PsiX));
	H5CHK(H5LTset_attribute_string(file, PSIX_NAME, "UNITS", "1"));
	H5CHK(dsrho = H5Dopen2(file, PSIX_NAME, H5P_DEFAULT));
	// H5CHK(dsrho = H5Dcreate2(file, PSIX_NAME, H5T_NATIVE_DOUBLE, dsp1d, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
	// H5CHK(H5Dwrite(dsrho, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, PsiX));
	H5CHK(H5DSset_scale(dsrho, PSIX_NAME));
	//H5CHK(H5LTset_attribute_string(dsrho, ".", "UNITS", "1"));

	TELL_FOLDER("Opened 3rd time");

	/* Psi */
	MULTI;
	HDFWrite1D(Psi, PSI_1D, "Wb", file, dims, dsrho);

	/* Pressure */
	MULTI;
	HDFWrite1D(P, PRESS_1D, "Pa", file, dims, dsrho);
	
	/* G */
	MULTI;
	HDFWrite1D(G, G_1D, "Wb", file, dims, dsrho);

	/* dP_dPsi */
	MULTI;
	HDFWrite1D(Pp, PP_1D, "Pa/Wb", file, dims, dsrho);

	/* g2p */
	MULTI;
	HDFWrite1D(G2p, G2P_1D, "Wb", file, dims, dsrho);

	/* dV/dPsi */
	MULTI;
	HDFWrite1D(dVdpsi, V_1D, "m^3/Wb", file, dims, dsrho);

	/* Vol */
	MULTI;
	HDFWrite1D(Vol, VOL_1D, "m^3", file, dims, dsrho);

	/* Shear */
	MULTI;
	HDFWrite1D(Shear, SHEAR_1D, "1/m^2", file, dims, dsrho);
	
	/* Well */
	MULTI;
	HDFWrite1D(Well, WELL_1D, "Wb", file, dims, dsrho);

	/* J ave */
	MULTI;
	HDFWrite1D(Jave, J_1D, "A/m^2", file, dims, dsrho);

	/* B^2 ave */
	MULTI;
	HDFWrite1D(B2ave, B2_1D, "T^2", file, dims, dsrho);

	/* Beta ave */
	MULTI;
	HDFWrite1D(Beta, BETA_1D, "1", file, dims, dsrho);

	/* Close the file */
	H5CHK(H5Dclose(dsrho));
	H5CHK(H5Fclose(file));

}


/*
**	HDFBoundary
**
*/
void    HDFBoundary(const char *Oname, const char *Vname, double psiLabel, 
                     double *X, double *Z, int len)
{
	char fname[FILENAME_MAX] = "";
	hid_t file, dsp1d;
	hsize_t dims[2];
	double *a;

	strncpy(fname, Oname, 18);
	strcat(fname, ".h5");

	TELL_FOLDER("Just about to open");

	/* S E T U P    H D F */
	MULTI;
	H5CHK(file = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT));

	dims[0] = 	len;
	dims[1] = 	2;

	TELL_FOLDER("Opened 4th or 5th time");

	/* Boundary */
	MULTI;
	a = (double *) malloc(len * 2 * sizeof(double));
	if (!a)
		nrerror("Failure to allocate temporary storarge in HDFBoundary.");
	for (int i = 0; i < len; i++)
	{
		a[i*2]   = X[i];
		a[i*2+1] = Z[i];
	}

	H5CHK(H5LTmake_dataset_double(file, Vname, 2, dims, a)); /* actually wants c-contigous memory representation */
	H5CHK(H5LTset_attribute_double(file, Vname, "PsiNormalized", &psiLabel, 1));
	H5CHK(H5LTset_attribute_string(file, Vname, "UNITS", "m"));
	H5CHK(H5LTset_attribute_string(file, Vname, "DIMENSION", "cartesian"));
	H5CHK(H5LTset_attribute_string(file, Vname, "FORMAT", "F7.4"));

	/* Close the file */
	H5CHK(H5Fclose(file));

	free(a);

	TELL_FOLDER("After closing file");	
	
	MULTI;

	TELL_FOLDER("Just reset the folder");
}



