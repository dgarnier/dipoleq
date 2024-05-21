/*
** Dipole 0.9
**
** HDFInput.c
**
** Routines to write binary HDF files for graphical
** data analysis.
**
** File:		HDFInput.c
** Date:		Jan 6, 1999
**
** Routine list:
**
**
**
** (c) D. Garnier -- Columbia University
*/

#include <stdlib.h>
#include "VAX_Alloc.h"
#include <stdio.h>
#include <string.h>
#include "nrutil.h"

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

#include "mfhdf.h"

#include "psigrid.h"
#include "plasma.h"
#include "HDFInput.h"

#define SDCHK(x)	if ((x) == FAIL) { \
						 fprintf(stdout,"when calling %s\n",#x); \
						 HEprint(stdout,0); \
                         nrerror("ERROR: SD-HDF error in HDFOutput.");  \
                     }

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307



float32		*AryFloat32(int nmax);
float32		*VecFloat32(int nmax);
void		Float32ToAry(double **da, float32 *sa, int nmax, double multiplier);
void		Float32ToVec(double *dv, float32 *sv, int nmax, double multiplier);
double   	*Float32ToVecP(float32 *sv, int nmax, double multiplier);


/*
**	AryFloat32
**
**
*/
float32      *AryFloat32(int nmax)
{
	float32      *a0;

	a0 =  (float32 *) malloc((unsigned) ((nmax + 1) * (nmax + 1)) * sizeof(float32));
	if (!a0)
		nrerror("ERROR: Allocation failure in AryFloat32.");

	return a0;
}

float32		*VecFloat32(int nmax)
{
	float32		*v;

	v = (float32 *) malloc((unsigned) (nmax + 1) * sizeof(float32));
	if (!v)
		nrerror("ERROR: Allocation failure in VecToFloat32.");

	return v;
}

void		Float32ToAry(double **da, float32 *sa, int nmax, double multiplier)
{
	int ix,iz;
	for (iz = 0; iz <= nmax; iz++)	/* slowest */
		for (ix = 0; ix <= nmax; ix++)	/* fastest */
			da[ix][iz] = ((double) *(sa++)) * multiplier;
}


/*
**	Float32ToVec
**
**
*/

void		Float32ToVec(double *dv, float32 *sv, int nmax, double multiplier)
{
	int           ix;
	for (ix = 0; ix <= nmax; ix++)
		dv[ix] = (double) *sv++ * multiplier;
}

double   *Float32ToVecP(float32 *sv, int nmax, double multiplier)
{
	double *dv;
	int           ix;

	dv = dvector(0, nmax);
	for (ix = 0; ix <= nmax; ix++)
		dv[ix] = (double) *sv++ * multiplier;

	return dv;

}

/*
**	HDFPsiGrid
**
**
*/

PSIGRID	*HDFPsiGridIn(char *Iname)
{
	PSIGRID		*pg;
	int         nmax;
	int32       starts[2], rank = 2;
	float32     *data, *v, maxv, minv;
	int32		sd_id, status, sds_id, sds_index, dim_id;

	int32 		dim_sizes[MAX_VAR_DIMS];
    int32 		data_type, n_attrs;
    char        name[MAX_NC_NAME];


	pg = new_PsiGrid();

	starts[0]=0;
	starts[1]=0;


	/* S E T U P    H D F */


	SDCHK(sd_id = SDstart(Iname,DFACC_READ));

	SDCHK(sds_index = SDnametoindex(sd_id, PSI_NAME));

	SDCHK(sds_id = SDselect(sd_id,sds_index));

	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes,
                           &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);

    pg->Nsize = nmax = dim_sizes[0]-1;

    printf("Setting Grid Size to %d\n",nmax);

	/* initialize PsiGrid */
	init_PsiGrid(pg);

	/* read in Psi */
	data = AryFloat32(nmax);
	status = SDreaddata(sds_id,starts,NULL,dim_sizes, (VOIDP) data);
	Float32ToAry(pg->Psi,data,nmax,1.0);
	SDCHK(status);


		/* also store the range from PsiAxis to PsiLim (PsiSep) */

	status = SDgetrange(sds_id,(VOIDP) &maxv, (VOIDP) &minv);

	pg->PsiLim  = (double) maxv;
	pg->PsiAxis = (double) minv;
	pg->DelPsi = pg->PsiLim - pg->PsiAxis;


	/* Define X dimension */
	dim_id = SDgetdimid(sds_id, 1);
	SDCHK(dim_id);

	v = VecFloat32(nmax);
	status = SDgetdimscale(dim_id, (VOIDP) v);
	Float32ToVec(pg->X,v,nmax,1.0);
	SDCHK(status);

	pg->Xmin = pg->X[0];
	pg->Xmax = pg->X[nmax];
	pg->dx = (pg->Xmax - pg->Xmin)/nmax;


	/* now define Z dimension */
	dim_id = SDgetdimid(sds_id, 0);
	SDCHK(dim_id);
	status = SDgetdimscale(dim_id, (VOIDP) v);
	Float32ToVec(pg->Z,v,nmax,1.0);
	SDCHK(status);
	free(v);

	pg->Zmin = pg->Z[0];
	pg->Zmax = pg->Z[nmax];
	pg->dz = (pg->Zmax - pg->Zmin)/nmax;

	status = SDendaccess (sds_id);
	SDCHK(status);

	/* Current */

	sds_index = SDnametoindex(sd_id, CUR_NAME);
	SDCHK(sds_index);

	sds_id = SDselect(sd_id,sds_index);
	SDCHK(sds_id);

	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes,
                           &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);

	status = SDreaddata(sds_id,starts,NULL,dim_sizes, (VOIDP) data);
	Float32ToAry(pg->Current,data,nmax,MU0);
	SDCHK(status);

	status = SDendaccess (sds_id);
	SDCHK(status);



	/* Residuals */


	sds_index = SDnametoindex(sd_id, RES_NAME);
	SDCHK(sds_index);

	sds_id = SDselect(sd_id,sds_index);
	SDCHK(sds_id);

	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes,
                           &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);

	status = SDreaddata(sds_id,starts,NULL,dim_sizes, (VOIDP) data);
	Float32ToAry(pg->Residual,data,nmax,MU0);
	SDCHK(status);

	status = SDendaccess (sds_id);
	SDCHK(status);

	status = SDend (sd_id);
	free(data);

	return pg;

}

/*
**	HDFPlasmaIn
**
**
*/
PLASMA *     HDFPlasmaIn( PSIGRID * pg, char *Oname)
{
	PLASMA		*pl;
	int32         dims[2], starts[2], nmax, rank = 2;
	float32      *data;
	int32		sd_id, status, sds_id, sds_index;
	int			ix,iz;

	int32 		dim_sizes[MAX_VAR_DIMS];
    int32 		data_type, n_attrs;
    char        name[MAX_NC_NAME];


	nmax = pg->Nsize;


	pl = new_Plasma();
	pl->Nsize = nmax;
	pl->ModelType = Plasma_Std;
//	init_Plasma(pl);

	pl->B2 = dmatrix(0, nmax, 0, nmax);
	pl->GradPsiX = dmatrix(0, nmax, 0, nmax);
	pl->GradPsiZ = dmatrix(0, nmax, 0, nmax);
	pl->GradPsi2 = dmatrix(0, nmax, 0, nmax);

	pl->Bt = dmatrix(0, nmax, 0, nmax);
	pl->G = dmatrix(0, nmax, 0, nmax);
    pl->Piso = dmatrix(0, nmax, 0, nmax);


	dims[0] = nmax + 1;
	dims[1] = nmax + 1;
	starts[0]=0;
	starts[1]=0;

	/* S E T U P    H D F */

	sd_id = SDstart(Oname,DFACC_READ);  /* concat into existing file */
	SDCHK(sd_id);


	/* B 2 */

	sds_index = SDnametoindex(sd_id, MODB_NAME);
	SDCHK(sds_index);

	sds_id = SDselect(sd_id,sds_index);
	SDCHK(sds_id);

	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes,
                           &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);

	data = AryFloat32(nmax);
	status = SDreaddata(sds_id,starts,NULL,dim_sizes, (VOIDP) data);
	Float32ToAry(pl->B2,data,nmax,1.0);
	free(data);
	SDCHK(status);

	status = SDendaccess (sds_id);
	SDCHK(status);


	/* Bp_x */

	sds_index = SDnametoindex(sd_id, BpX_NAME);
	SDCHK(sds_index);

	sds_id = SDselect(sd_id,sds_index);
	SDCHK(sds_id);

	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes,
                           &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);

	data = AryFloat32(nmax);
	status = SDreaddata(sds_id,starts,NULL,dim_sizes, (VOIDP) data);
	Float32ToAry(pl->GradPsiZ,data,nmax,TWOPI);
	free(data);
	SDCHK(status);

	for (ix = 0; ix <= nmax; ix++)
		for (iz = 0; iz <= nmax; iz++)
			pl->GradPsiZ[ix][iz] = pl->GradPsiZ[ix][iz] * pg->X[ix];

	status = SDendaccess (sds_id);
	SDCHK(status);

	/* Bp_z */

	sds_index = SDnametoindex(sd_id, BpZ_NAME);
	SDCHK(sds_index);

	sds_id = SDselect(sd_id,sds_index);
	SDCHK(sds_id);

	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes,
                           &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);

	data = AryFloat32(nmax);
	status = SDreaddata(sds_id,starts,NULL,dim_sizes, (VOIDP) data);
	Float32ToAry(pl->GradPsiX,data,nmax,-TWOPI);
	free(data);
	SDCHK(status);

	for (ix = 0; ix <= nmax; ix++)
		for (iz = 0; iz <= nmax; iz++)
			pl->GradPsiX[ix][iz] = pl->GradPsiX[ix][iz] * pg->X[ix];

	status = SDendaccess (sds_id);
	SDCHK(status);


	/* G */

	sds_index = SDnametoindex(sd_id, TFLUX_NAME);
	SDCHK(sds_index);

	sds_id = SDselect(sd_id,sds_index);
	SDCHK(sds_id);

	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes,
                           &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);

	data = AryFloat32(nmax);
	status = SDreaddata(sds_id,starts,NULL,dim_sizes, (VOIDP) data);
	Float32ToAry(pl->G,data,nmax,1.0);
	free(data);
	SDCHK(status);

	status = SDendaccess (sds_id);
	SDCHK(status);


	/* P R E S S U R E */
	sds_index = SDnametoindex(sd_id, PRESS_NAME);
	SDCHK(sds_index);

	sds_id = SDselect(sd_id,sds_index);
	SDCHK(sds_id);

	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes,
                           &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);

	data = AryFloat32(nmax);
	status = SDreaddata(sds_id,starts,NULL,dim_sizes, (VOIDP) data);
	Float32ToAry(pl->Piso,data,nmax,1.0);
	free(data);
	SDCHK(status);

	status = SDendaccess (sds_id);
	SDCHK(status);

	/* end */
	status = SDend (sd_id);
	SDCHK(status);

	return pl;

}



/*
**	HDFFluxFuncs
**
**
*/
void	HDFFluxFuncsIn(char *Oname, int *nptsIn, double **PsiX,
					double **Psi, double **P, double **G, double **Pp, double **G2p,
					double **q, double **dVdpsi, double **Vol, double **Shear,
					double **Well, double **Jave, double **B2ave, double **Beta)
{
	int32         dim[1], start[1];
	float32      *a;
	int32		sd_id, sds_id, dim_id, sds_index;
	int			npts;

	int32 		dim_sizes[MAX_VAR_DIMS];
    int32 		data_type, n_attrs, rank;
    char        name[MAX_NC_NAME];




	/* S E T U P    H D F */

	SDCHK( sd_id = SDstart(Oname,DFACC_READ) );  /* concat into existing file */

	/* Psi */

	SDCHK(sds_index = SDnametoindex(sd_id, PSI_1D));
	SDCHK(sds_id = SDselect(sd_id,sds_index));
	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes,
                           &data_type, &n_attrs));
    printf("Found %ld points for 1D data", dim_sizes[0]);
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);

    *nptsIn = npts = dim_sizes[0];
    a = (float32 *) malloc (npts*sizeof(float32));
    dim[0] = npts;
	start[0]=0;

	SDCHK(SDreaddata(sds_id,start,NULL,dim_sizes, (VOIDP) a));
	*Psi = Float32ToVecP(a, npts-1, 1.0);

	/* Define X dimension */
	SDCHK( dim_id = SDgetdimid(sds_id, 0));
	SDCHK( SDgetdimscale(dim_id, (VOIDP) a) );
	*PsiX = Float32ToVecP(a,npts-1,1.0);
	SDCHK( SDendaccess (sds_id) );

	/* Pressure */
	SDCHK(sds_index = SDnametoindex(sd_id, PRESS_1D));
	SDCHK(sds_id = SDselect(sd_id,sds_index));
	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);
   	SDCHK(SDreaddata(sds_id,start,NULL,dim_sizes, (VOIDP) a));
	*P = Float32ToVecP(a, npts-1, 1.0);
	SDCHK( SDendaccess (sds_id) );

	/* G */
	SDCHK(sds_index = SDnametoindex(sd_id, G_1D));
	SDCHK(sds_id = SDselect(sd_id,sds_index));
	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);
   	SDCHK(SDreaddata(sds_id,start,NULL,dim_sizes, (VOIDP) a));
	*G = Float32ToVecP(a, npts-1, 1.0);
	SDCHK( SDendaccess (sds_id) );


	/* dP_dPsi */
	SDCHK(sds_index = SDnametoindex(sd_id, PP_1D));
	SDCHK(sds_id = SDselect(sd_id,sds_index));
	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);
   	SDCHK(SDreaddata(sds_id,start,NULL,dim_sizes, (VOIDP) a));
	*Pp = Float32ToVecP(a, npts-1, 1.0);
	SDCHK( SDendaccess (sds_id) );


	/* g2p */
	SDCHK(sds_index = SDnametoindex(sd_id, G2P_1D));
	SDCHK(sds_id = SDselect(sd_id,sds_index));
	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);
   	SDCHK(SDreaddata(sds_id,start,NULL,dim_sizes, (VOIDP) a));
	*G2p = Float32ToVecP(a, npts-1, 1.0);
	SDCHK( SDendaccess (sds_id) );


	/* q */
	SDCHK(sds_index = SDnametoindex(sd_id, Q_1D));
	SDCHK(sds_id = SDselect(sd_id,sds_index));
	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);
   	SDCHK(SDreaddata(sds_id,start,NULL,dim_sizes, (VOIDP) a));
	*q = Float32ToVecP(a, npts-1, 1.0);
	SDCHK( SDendaccess (sds_id) );


	/* dV/dPsi */
	SDCHK(sds_index = SDnametoindex(sd_id, V_1D));
	SDCHK(sds_id = SDselect(sd_id,sds_index));
	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);
   	SDCHK(SDreaddata(sds_id,start,NULL,dim_sizes, (VOIDP) a));
	*dVdpsi = Float32ToVecP(a, npts-1, 1.0);
	SDCHK( SDendaccess (sds_id) );


	/* Vol */
	SDCHK(sds_index = SDnametoindex(sd_id, VOL_1D));
	SDCHK(sds_id = SDselect(sd_id,sds_index));
	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);
   	SDCHK(SDreaddata(sds_id,start,NULL,dim_sizes, (VOIDP) a));
	*Vol = Float32ToVecP(a, npts-1, 1.0);
	SDCHK( SDendaccess (sds_id) );


	/* Shear */
	SDCHK(sds_index = SDnametoindex(sd_id, SHEAR_1D));
	SDCHK(sds_id = SDselect(sd_id,sds_index));
	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);
   	SDCHK(SDreaddata(sds_id,start,NULL,dim_sizes, (VOIDP) a));
	*Shear = Float32ToVecP(a, npts-1, 1.0);
	SDCHK( SDendaccess (sds_id) );


	/* Well */
	SDCHK(sds_index = SDnametoindex(sd_id, WELL_1D));
	SDCHK(sds_id = SDselect(sd_id,sds_index));
	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);
   	SDCHK(SDreaddata(sds_id,start,NULL,dim_sizes, (VOIDP) a));
	*Well = Float32ToVecP(a, npts-1, 1.0);
	SDCHK( SDendaccess (sds_id) );

	/* B^2 ave */
	SDCHK(sds_index = SDnametoindex(sd_id, B2_1D));
	SDCHK(sds_id = SDselect(sd_id,sds_index));
	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);
   	SDCHK(SDreaddata(sds_id,start,NULL,dim_sizes, (VOIDP) a));
	*B2ave = Float32ToVecP(a, npts-1, 1.0);
	SDCHK( SDendaccess (sds_id) );

	/* Beta ave */
	SDCHK(sds_index = SDnametoindex(sd_id, BETA_1D));
	SDCHK(sds_id = SDselect(sd_id,sds_index));
	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);
   	SDCHK(SDreaddata(sds_id,start,NULL,dim_sizes, (VOIDP) a));
	*Beta = Float32ToVecP(a, npts-1, 1.0);
	SDCHK( SDendaccess (sds_id) );

	/*J ave */
	SDCHK(sds_index = SDnametoindex(sd_id, J_1D));
	SDCHK(sds_id = SDselect(sd_id,sds_index));
	SDCHK(SDgetinfo (sds_id, name, &rank, dim_sizes, &data_type, &n_attrs));
    printf("Data Set %s: rank %ld, data_type %ld, n_attrs %ld\n",name,rank,data_type,n_attrs);
   	SDCHK(SDreaddata(sds_id,start,NULL,dim_sizes, (VOIDP) a));
	*Jave = Float32ToVecP(a, npts-1, 1.0);
	SDCHK( SDendaccess (sds_id) );


	SDCHK( SDend (sd_id) );

	free(a);
}
