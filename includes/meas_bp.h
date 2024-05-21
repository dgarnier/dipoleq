/*
** TokaMac v2.0
**
** meas_bp.h
**
** Header file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_bp_Green(TOKAMAK *td, MEAS *m)
**
**				meas_bp_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_bp_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_bp.h
** Date:		March 24, 1993
**
**	USAGE:
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _MEAS_BP_

#define	_MEAS_BP_		1

#include "measurement.h"
#include "tokamak.h"

#ifdef __cplusplus
extern "C" {
#endif

double  g_bp(double x, double z, double xc, double zc);
void	meas_bp_Green(TOKAMAK *, MEAS *);
void	meas_bp_Fit(TOKAMAK *, MEAS *);
void	meas_bp_Now(TOKAMAK *, MEAS *);
void	meas_bp_L(TOKAMAK *, MEAS *, double *);

#ifdef __cplusplus
}
#endif
#endif
