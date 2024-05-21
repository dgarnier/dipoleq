/*
** TokaMac v2.0
**
** meas_ppsix.h
**
** Header file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_ppsix_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_ppsix_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_ppsix.h
** Date:		October 26, 1993
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _MEAS_PPSIX_

#define	_MEAS_PPSIX_		1

#include "measurement.h"
#include "tokamak.h"

void	meas_ppsix_Fit(TOKAMAK *, MEAS *);
void	meas_ppsix_Now(TOKAMAK *, MEAS *);
void	meas_ppsix_L(TOKAMAK *, MEAS *, double *);

#endif
