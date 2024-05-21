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
#include "limiter.h"
#include "HDFOutput.h"

#define H5CHK(x)	if ((x) < 0) do { \
						H5Eprint2(H5E_DEFAULT, stderr); \
						nrerror("ERROR: HDF5 error in HDFOutput."); \
					} while (0)


#define PI	3.14159265358979323
#define MU0	1.25663706e-06
#define TWOPI	6.283185307

#define ITER_IJ(ix, iz, nmin, nmax) \
	for (int iz = nmin; iz <= nmax; iz++)	  /* slowest */ \
		for (int ix = nmin; ix <= nmax; ix++) /* fastest */ \


/*
**	ScaleArray
**
**
*/
static void ScaleArray(double *a, double **ad, int nmax, double multiplier)
{
	ITER_IJ(ix, iz, 0, nmax) *a++ = ad[ix][iz] * multiplier;
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

static void HDFWrite0D(double * data, const char *name, const char *units, hid_t loc)
{
	H5CHK(H5LTmake_dataset_double(loc, name, 0, NULL, data));
	H5CHK(H5LTset_attribute_string(loc, name, "UNITS", units));
}

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
	hid_t  		file, dsp1d, dsp2d, dsc1, dsc2, g0d, g1d, g2d, gbd;
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
	H5CHK(H5LTset_attribute_string(file, ".", "TITLE", "Dipole 0.9 Equilibrium Data"));
	H5CHK(H5LTset_attribute_string(file, ".", "VERSION", "0.9"));
	H5CHK(H5LTset_attribute_string(file, ".", "ONAME", Oname));

	/* Create 3 groups, 0d, 1d, and 2d */
	H5CHK(g0d = H5Gcreate2(file, SCALAR_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
	H5CHK(g1d = H5Gcreate2(file, FLUX_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
	H5CHK(g2d = H5Gcreate2(file, GRID_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
	H5CHK(gbd = H5Gcreate2(file, BOUND_GROUP, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
	H5CHK(H5Gclose(g1d)); /* close these now */
	H5CHK(H5Gclose(gbd));

	/* 0D first, (magnetic axis) */
	HDFWrite0D(&pg->XMagAxis, RMAGX_0D, "m", g0d);
	HDFWrite0D(&pg->ZMagAxis, ZMAGX_0D, "m", g0d);

	/* Do dimensions first, aka Dimension Scales */

	/* dataspaces */
	H5CHK(dsp1d = H5Screate_simple(1, dims, NULL));
	H5CHK(dsp2d = H5Screate_simple(2, dims, NULL));

	/* dimension datasets -> datascale datasets*/
	H5CHK(dsc1 = H5Dcreate2(g2d, DIMX_NAME, H5T_NATIVE_DOUBLE, dsp1d, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
	H5CHK(dsc2 = H5Dcreate2(g2d, DIMZ_NAME, H5T_NATIVE_DOUBLE, dsp1d, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
	H5CHK(H5Dwrite(dsc1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pg->X));
	H5CHK(H5Dwrite(dsc2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pg->Z));
	H5CHK(H5DSset_scale(dsc1, DIMX_NAME));
	H5CHK(H5DSset_scale(dsc2, DIMZ_NAME));
	H5CHK(H5LTset_attribute_string(dsc1, ".", "UNITS", "m"));
	H5CHK(H5LTset_attribute_string(dsc2, ".", "UNITS", "m"));

	/* Do Current */
	ScaleArray(a, pg->Current, nmax, 1.0 / MU0);
	HDFWrite2D(a, CUR_NAME, "A/m^2", g2d, dsp2d, dsc1, dsc2);

	/* Do Psi */
	ScaleArray(a, pg->Psi, nmax, 1.0);
	HDFWrite2D(a, PSI_NAME, "Wb", g2d, dsp2d, dsc1, dsc2);
	/* attributes can't have units, but they are in Wb */
	H5LTset_attribute_double(g2d, PSI_NAME, PSIFCFS_0D, &pg->PsiAxis, 1);
	H5LTset_attribute_double(g2d, PSI_NAME, PSILCFS_0D, &pg->PsiLim, 1);
	HDFWrite0D(&pg->PsiMagAxis, PSIAXIS_0D, "Wb", g0d);
	HDFWrite0D(&pg->PsiAxis, PSIFCFS_0D, "Wb", g0d);
	HDFWrite0D(&pg->PsiLim, PSILCFS_0D, "Wb", g0d);




	/* Do Residuals */
	ScaleArray(a, pg->Residual, nmax, 1.0 / MU0);
	HDFWrite2D(a, RES_NAME, "A/m^2", g2d, dsp2d, dsc1, dsc2);


	/* Close the file */
	H5CHK(H5Dclose(dsc1));
	H5CHK(H5Dclose(dsc2));
	H5CHK(H5Sclose(dsp1d));
	H5CHK(H5Sclose(dsp2d));
	H5CHK(H5Gclose(g0d));
	H5CHK(H5Gclose(g2d));
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
	hid_t file, dsp2d, dsc1, dsc2, g0d, g2d;
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
	H5CHK(g0d = H5Gopen2(file, SCALAR_GROUP, H5P_DEFAULT));
	H5CHK(g2d = H5Gopen2(file, GRID_GROUP, H5P_DEFAULT));

	/* dataspaces */
	H5CHK(dsp2d = H5Screate_simple(2, dims, NULL));

	/* read dimension scales */
	H5CHK(dsc1 = H5Dopen2(g2d, DIMX_NAME, H5P_DEFAULT));
	H5CHK(dsc2 = H5Dopen2(g2d, DIMZ_NAME, H5P_DEFAULT));


	TELL_FOLDER("Opened 2nd time");
	/* 0D values */
	HDFWrite0D(&pl->R0, R0_0D, "m", g0d);
	HDFWrite0D(&pl->Z0, Z0_0D, "m", g0d);
	HDFWrite0D(&pl->B0, BT_0D, "T", g0d);
	HDFWrite0D(&pl->Ip, IP_0D, "A", g0d);

	/* B 2 */
	MULTI;
	ScaleArray(a, pl->B2, nmax, 1.0);
	HDFWrite2D(a, MODB_NAME, "T^2", g2d, dsp2d, dsc1, dsc2);

	/* Bp_x */
	MULTI;
	ap = a;
	ITER_IJ(ix, iz, 0, nmax) *ap++ = pl->GradPsiZ[ix][iz] / TWOPI / pg->X[ix];
	HDFWrite2D(a, BpX_NAME, "T", g2d, dsp2d, dsc1, dsc2);

	/* Bp_z */
	MULTI;
	ap = a;
	ITER_IJ(ix, iz, 0, nmax) *ap++ = -pl->GradPsiX[ix][iz] / TWOPI / pg->X[ix];

	HDFWrite2D(a, BpZ_NAME, "T", g2d, dsp2d, dsc1, dsc2);

	/* G */
	MULTI;
	ScaleArray(a, pl->G, nmax, 1.0);
	HDFWrite2D(a, TFLUX_NAME, "1", g2d, dsp2d, dsc1, dsc2);

		/* P R E S S U R E */
	MULTI;
	switch (pl->ModelType) {
      default :
			ScaleArray(a, pl->Piso, nmax, 1.0 / MU0);
			HDFWrite2D(a, PRESS_NAME, "Pa", g2d, dsp2d, dsc1, dsc2);
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
			ITER_IJ(ix, iz, 0, nmax) *ap++ = 2*pl->Piso[ix][iz]/pl->B2[ix][iz];
			HDFWrite2D(a, BETA_NAME, "", g2d, dsp2d, dsc1, dsc2);
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
	H5CHK(H5Sclose(dsp2d));
	H5CHK(H5Gclose(g2d));
	H5CHK(H5Gclose(g0d));
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
					double *Well, double *Jave, double *B2ave, double *Beta,
					double *BetaMax, double *XBetaMax, double *ZBetaMax, double *BBetaMax,
					double *BMax, double *XBMax, double *ZBMax)
{
	char fname[FILENAME_MAX] = "";
	hid_t file, dsp1d, dsrho, g1d;
	hsize_t dims[1];

	dims[0] = npts;

	strncpy(fname, Oname, 18);
	strcat(fname, ".h5");

	TELL_FOLDER("Just about to open");

	/* S E T U P    H D F */
	MULTI;
	H5CHK(file = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT));
	H5CHK(g1d = H5Gopen2(file, FLUX_GROUP, H5P_DEFAULT));

	/* dataspaces */
	// H5CHK(dsp1d = H5Screate_simple(1, dims, NULL));

	/* psi_x dimension scales */

	H5CHK(H5LTmake_dataset_double(g1d, PSIX_NAME, 1, dims, PsiX));
	H5CHK(H5LTset_attribute_string(g1d, PSIX_NAME, "UNITS", "1"));
	H5CHK(dsrho = H5Dopen2(g1d, PSIX_NAME, H5P_DEFAULT));
	H5CHK(H5DSset_scale(dsrho, PSIX_NAME));

	TELL_FOLDER("Opened 3rd time");

	/* Psi */
	MULTI;
	HDFWrite1D(Psi, PSI_1D, "Wb", g1d, dims, dsrho);

	/* Pressure */
	MULTI;
	HDFWrite1D(P, PRESS_1D, "Pa", g1d, dims, dsrho);

	/* G */
	MULTI;
	HDFWrite1D(G, G_1D, "1", g1d, dims, dsrho);

	/* dP_dPsi */
	MULTI;
	HDFWrite1D(Pp, PP_1D, "Pa/Wb", g1d, dims, dsrho);

	/* g2p */
	MULTI;
	HDFWrite1D(G2p, G2P_1D, "1/Wb", g1d, dims, dsrho);

	/* dV/dPsi */
	MULTI;
	HDFWrite1D(dVdpsi, V_1D, "m^3/Wb", g1d, dims, dsrho);

	/* Vol */
	MULTI;
	HDFWrite1D(Vol, VOL_1D, "m^3", g1d, dims, dsrho);

	/* q */
	MULTI;
	HDFWrite1D(q, Q_1D, "", g1d, dims, dsrho);

	/* Shear */
	MULTI;
	HDFWrite1D(Shear, SHEAR_1D, "1/m^2", g1d, dims, dsrho);

	/* Well */
	MULTI;
	HDFWrite1D(Well, WELL_1D, "Wb", g1d, dims, dsrho);

	/* J ave */
	MULTI;
	HDFWrite1D(Jave, J_1D, "A/m^2", g1d, dims, dsrho);

	/* B^2 ave */
	MULTI;
	HDFWrite1D(B2ave, B2_1D, "T^2", g1d, dims, dsrho);

	/* Beta ave */
	MULTI;
	HDFWrite1D(Beta, BETA_1D, "", g1d, dims, dsrho);

	/* Beta max */
	MULTI;
	HDFWrite1D(BetaMax,    BETAMAX_1D, "", g1d, dims, dsrho);
	HDFWrite1D(XBetaMax, X_BETAMAX_1D, "m", g1d, dims, dsrho);
	HDFWrite1D(ZBetaMax, Z_BETAMAX_1D, "m", g1d, dims, dsrho);
	HDFWrite1D(BBetaMax, B_BETAMAX_1D, "T", g1d, dims, dsrho);

	/* B max */
	MULTI;
	HDFWrite1D(BMax,    BMAX_1D, "T", g1d, dims, dsrho);
	HDFWrite1D(XBMax, X_BMAX_1D, "m", g1d, dims, dsrho);
	HDFWrite1D(ZBMax, Z_BMAX_1D, "m", g1d, dims, dsrho);

	/* Close the file */
	H5CHK(H5Dclose(dsrho));
	H5CHK(H5Gclose(g1d));
	H5CHK(H5Fclose(file));

}

/*
** HDFLimiters
**
*/

void	HDFLimiters(const char *Oname, LIMITER **lims, int nlims)
{
	char fname[FILENAME_MAX] = "";
	hid_t file, dsp1d, gbd;
	hsize_t dims[3];
	int nsegp = 0, nsegn = 0;
	double *a, *ap;

	strncpy(fname, Oname, 18);
	strcat(fname, ".h5");

	TELL_FOLDER("Just about to open");

	/* S E T U P    H D F */
	MULTI;
	H5CHK(file = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT));
	H5CHK(gbd = H5Gopen2(file, BOUND_GROUP, H5P_DEFAULT));

	TELL_FOLDER("Opened 6th or time");

	for (int i = 0; i < nlims; i++)
	{
		if (lims[i]->Enabled > 0)
			nsegp++;
		if (lims[i]->Enabled < 0)
			nsegn++;
	}
	if (nsegp > 0) {
		/* Segments */
		dims[0] = nsegp;
		dims[1] = dims[2] = 2;
		a = ap = (double *) malloc(nsegp * 4 * sizeof(double));
		for (int i = 0; i < nlims; i++)
			if (lims[i]->Enabled > 0) {
				*ap++ = lims[i]->X1;
				*ap++ = lims[i]->Z1;
				*ap++ = lims[i]->X2;
				*ap++ = lims[i]->Z2;
			}
		H5CHK(H5LTmake_dataset_double(gbd, OLIM_NAME, 3, dims, a));
		H5CHK(H5LTset_attribute_string(gbd, OLIM_NAME, "UNITS", "m"));
		H5CHK(H5LTset_attribute_string(gbd, OLIM_NAME, "DIMENSION", "cylindrical"));
		H5CHK(H5LTset_attribute_string(gbd, OLIM_NAME, "FORMAT", "F7.4"));
		free(a);
	}
	if (nsegn > 0) {
		/* Segments */
		dims[0] = nsegn;
		dims[1] = dims[2] = 2;
		a = ap = (double *) malloc(nsegn * 4 * sizeof(double));
		for (int i = 0; i < nlims; i++)
			if (lims[i]->Enabled < 0) {
				*ap++ = lims[i]->X1;
				*ap++ = lims[i]->Z1;
				*ap++ = lims[i]->X2;
				*ap++ = lims[i]->Z2;
			}
		H5CHK(H5LTmake_dataset_double(gbd, ILIM_NAME, 3, dims, a));
		H5CHK(H5LTset_attribute_string(gbd, ILIM_NAME, "UNITS", "m"));
		H5CHK(H5LTset_attribute_string(gbd, ILIM_NAME, "DIMENSION", "cylindrical"));
		H5CHK(H5LTset_attribute_string(gbd, ILIM_NAME, "FORMAT", "F7.4"));
		free(a);
	}

	/* Close the file */
	H5CHK(H5Gclose(gbd));
	H5CHK(H5Fclose(file));

	TELL_FOLDER("After closing file");
}

/*
**	HDFBoundary
**
*/
void    HDFBoundary(const char *Oname, const char *Vname, double psiLabel,
                     double *X, double *Z, int len)
{
	char fname[FILENAME_MAX] = "";
	hid_t file, dsp1d, gid;
	hsize_t dims[2];
	double *a;

	strncpy(fname, Oname, 18);
	strcat(fname, ".h5");

	TELL_FOLDER("Just about to open");

	/* S E T U P    H D F */
	MULTI;
	H5CHK(file = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT));
	H5CHK(gid = H5Gopen2(file, BOUND_GROUP, H5P_DEFAULT));

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

	H5CHK(H5LTmake_dataset_double(gid, Vname, 2, dims, a)); /* actually wants c-contigous memory representation */
	H5CHK(H5LTset_attribute_double(gid, Vname, "PsiNormalized", &psiLabel, 1));
	H5CHK(H5LTset_attribute_string(gid, Vname, "UNITS", "m"));
	H5CHK(H5LTset_attribute_string(gid, Vname, "DIMENSION", "cartesian"));
	H5CHK(H5LTset_attribute_string(gid, Vname, "FORMAT", "F7.4"));

	/* Close the file */
	H5CHK(H5Gclose(gid));
	H5CHK(H5Fclose(file));

	free(a);

	TELL_FOLDER("After closing file");

	MULTI;

	TELL_FOLDER("Just reset the folder");
}
