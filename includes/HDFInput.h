/*
** Dipole v0.9
**
** HDFInput.h
**
** Routines to read binary HDF files for post-processing
** data analysis.
**
** File:		HDFInput.h
** Date:		Jan 6, 1998
**
** Routine list:
**
**
**
** (c) D. Garnier -- Columbia University
*/

#ifndef _HDFINPUT_

#define _HDFINPUT_		1

#include "psigrid.h"
#include "plasma.h"
#include "HDFOutput.h"

#ifdef __cplusplus
extern "C" {
#endif

PSIGRID 	*HDFPsiGridIn(char *);
PLASMA		*HDFPlasmaIn(PSIGRID *,char *);
void		HDFFluxFuncsIn(char *Oname, int *nptsIn, double **PsiX,
					double **Psi, double **P, double **G, double **Pp, double **G2p,
					double **q, double **dVdpsi, double **Vol, double **Shear,
					double **Well, double **Jave, double **B2ave, double **Beta);

#ifdef __cplusplus
}
#endif

#endif
