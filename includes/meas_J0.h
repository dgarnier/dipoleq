/*
** TokaMac v2.0
**
** meas_J0.h
**
** Header file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_J0_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_J0_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_J0.h
** Date:		February 21, 1996
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _MEAS_J0_

#define	_MEAS_J0_		1

#include "measurement.h"
#include "tokamak.h"

void	meas_J0_Fit(TOKAMAK *, MEAS *);
void	meas_J0_Now(TOKAMAK *, MEAS *);
void	meas_J0_L(TOKAMAK *, MEAS *, double *);

#endif
