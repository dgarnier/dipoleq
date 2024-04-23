/*
** TokaMac v2.0
**
** HDFOutput.h
**
** Routines to write binary HDF files for graphical
** data analysis.
**
** File:		HDFOutput.h
** Date:		April 2, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _HDFOUTPUT_

#define _HDFOUTPUT_		1

#include "psigrid.h"
#include "plasma.h"
#include "limiter.h"

/* Groups */
#define GRID_GROUP	"/Grid"
#define BOUND_GROUP	"/Boundaries"
#define FLUX_GROUP	"/FluxFunctions"
#define SCALAR_GROUP "/Scalars"

/* DataSets */
#define CUR_NAME 	"Current"
#define PSI_NAME 	"Psi"
#define RES_NAME 	"Residuals"
#define DIMX_NAME 	"R"
#define DIMZ_NAME 	"Z"
#define PSIX_NAME 	"PsiX"
#define MODB_NAME 	"Mod B2"
#define BpX_NAME  	"Bp_R"
#define BpZ_NAME  	"Bp_Z"
#define TFLUX_NAME 	"ToroidalFlux"
#define PRESS_NAME 	"Pressure"
#define BETA_NAME 	"Beta"
#define LCFS_NAME	"LCFS"
#define FCFS_NAME	"FCFS"
#define PSI_1D 		"Psi_1D"
#define PRESS_1D 	"pres"
#define G_1D 		"fpol"
#define PP_1D 		"pprime"
#define G2P_1D 		"dG2dPsi_1D"
#define FFP_1D 		"ffprime"
#define Q_1D 		"qpsi"
#define V_1D 		"dVdPsi_1D"
#define VOL_1D 		"Vol_1D"
#define SHEAR_1D 	"Shear_1D"
#define WELL_1D 	"Well_1D"
#define B2_1D 		"B2ave_1D"
#define BETA_1D 	"Beta_1D"
#define J_1D 		"Jave_1D"
#define IP_0D 		"cpasma" // EFIT name for plasma current
#define BT_0D 		"bcentr" // EFIT name for toroidal field
#define R0_0D 		"rcentr" // EFIT name for reference radius
#define PSIAXIS_0D 	"simagx" // EFIT name for psi axis
#define PSILIM_0D   "sibdry" // EFIT name for psi axis
#define RMAGX_0D 	"rmagx" // EFIT name for R axis
#define ZMAGX_0D 	"zmagx" // EFIT name for Z axis
#define OLIM_NAME   "lim"   // EFIT name for limiter
#define ILIM_NAME   "ilim"  // inner limiter
  

#ifdef __cplusplus
extern "C" {
#endif

void	HDFPsiGrid(PSIGRID *,char *);
void 	HDFLimiters(const char *Oname, LIMITER **lims, int nlims);
void 	HDFBoundary(const char *Oname, const char *Vname, double, double *, double *, int);
void    HDFPPsi(PSIGRID *, char *);
void	HDFPlasma(PLASMA *,PSIGRID *,char *);
void	HDFFluxFuncs(char *Oname, int npts, double *PsiX, 
					double *Psi, double *P, double *G, double *Pp, double *G2p, 
					double *q, double *dVdpsi, double *Vol, double *Shear, 
					double *Well, double *Jave, double *B2ave, double *Beta);

#ifdef __cplusplus
}
#endif
#endif

