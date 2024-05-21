/*
** Dipole v0.9
**
** J_DipoleStd.h
**
**
**
** File:		J_DipoleStd.h
** Date:		December 11, 1998
**
** Revisions:
**
**
**
**
** (c) D. Garnier -- Columbia University
*/

#ifndef _J_DIPOLESTD_

#define _J_DIPOLESTD_		1

#include "tokamak.h"

#ifdef __cplusplus
extern "C" {
#endif

void 		J_DipoleStd(TOKAMAK *,double **);
double 		J_DipoleStd_Loc(TOKAMAK * , int , int);
double		Pp_DipoleStd_Loc(TOKAMAK *, double Psi);
double		P_DipoleStd_Loc(TOKAMAK *, double Psi);
double		G2p_DipoleStd_Loc(TOKAMAK *, double Psi);
double		G2_DipoleStd_Loc(TOKAMAK *, double Psi);

#ifdef __cplusplus
}
#endif
#endif
