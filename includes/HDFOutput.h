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
#define PSIX_NAME 	"PsiNorm"
#define MODB_NAME 	"B2"
#define BpX_NAME  	"Bp_R"
#define BpZ_NAME  	"Bp_Z"
#define TFLUX_NAME 	"ToroidalFlux"
#define PRESS_NAME 	"Pressure"
#define BETA_NAME 	"Beta"
#define LCFS_NAME	"LCFS"
#define FCFS_NAME	"FCFS"
#define PSI_1D 		"psi"
#define PRESS_1D 	"ppsi"
#define G_1D 		"Gpsi"
#define PP_1D 		"pprime"
#define G2P_1D 		"G2prime"
#define FFP_1D 		"ffprime"
#define Q_1D 		"qpsi"
#define V_1D 		"Vprime"
#define VOL_1D 		"Vpsi"
#define SHEAR_1D 	"Shear"
#define WELL_1D 	"Well"
#define B2_1D 		"B2Ave"
#define BETA_1D 	"BetaAve"
#define J_1D 		"JAve"
#define BETAMAX_1D  "BetaMax"
#define X_BETAMAX_1D "RBetaMax"
#define Z_BETAMAX_1D "ZBetaMax"
#define B_BETAMAX_1D "BBetaMax"
#define BMAX_1D 	"BMax"
#define X_BMAX_1D 	"RBMax"
#define Z_BMAX_1D 	"ZBMax"

#define IP_0D 		"Ip"      // name for plasma current
#define BT_0D 		"B0"      // EFIT name for toroidal field
#define R0_0D 		"R0"      // EFIT name for reference radius
#define Z0_0D 		"Z0"      // EFIT name for reference radius
#define R0Z0_0D 	"R0Z0"    // Scale factor for G -> F conversion
#define PSIAXIS_0D 	"PsiMagX" // psi on magnetic axis
#define PSIFCFS_0D 	"PsiFCFS" // psi at boundary near or at magnetic axis
#define PSILCFS_0D  "PsiLCFS" // psi at boundary at outer limiter or separatrix
#define RMAGX_0D 	"RMagX"   // R axis
#define ZMAGX_0D 	"ZMagX"   // Z axis
#define OLIM_NAME   "olim"    // EFIT name for limiter
#define ILIM_NAME   "ilim"    // inner limiter


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
					double *Well, double *Jave, double *B2ave, double *Beta,
					double *BetaMax, double *XBetaMax, double *ZBetaMax, double *BBetaMax,
					double *BMax, double *XBMax, double *ZBMax);

#ifdef __cplusplus
}
#endif
#endif
