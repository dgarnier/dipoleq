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

#define CUR_NAME 	"Current"
#define PSI_NAME 	"Psi"
#define RES_NAME 	"Residuals"
#define DIMX_NAME 	"x"
#define DIMZ_NAME 	"z"
#define PSIX_NAME 	"PsiX"
#define MODB_NAME 	"Mod B2"
#define BpX_NAME  	"Bp_x"
#define BpZ_NAME  	"Bp_z"
#define TFLUX_NAME 	"ToroidalFlux"
#define PRESS_NAME 	"Pressure"
#define BETA_NAME 	"Beta"
#define LCFS_NAME	"LCFSBoundary"
#define FCFS_NAME	"FCFSBoundary"
#define PSI_1D 		"Psi_1D"
#define PRESS_1D 	"Pressure_1D"
#define G_1D 		"G_1D"
#define PP_1D 		"dP/dPsi_1D"
#define G2P_1D 		"dG2/dPsi_1D"
#define Q_1D 		"q_1D"
#define V_1D 		"dV/dPsi_1D"
#define VOL_1D 		"Vol_1D"
#define SHEAR_1D 	"Shear_1D"
#define WELL_1D 	"Well_1D"
#define B2_1D 		"B2ave_1D"
#define BETA_1D 	"Beta_1D"
#define J_1D 		"Jave_1D"

void	HDFPsiGrid(PSIGRID *,char *);
void	HDFBoundary(const char *Oname, const char *Vname, double, double *, double *, int);
void    HDFPPsi(PSIGRID *, char *);
void	HDFPlasma(PLASMA *,PSIGRID *,char *);
void	HDFFluxFuncs(char *Oname, int npts, double *PsiX, 
					double *Psi, double *P, double *G, double *Pp, double *G2p, 
					double *q, double *dVdpsi, double *Vol, double *Shear, 
					double *Well, double *Jave, double *B2ave, double *Beta);


#endif

