/*
** TokaMac v2.0
**
** meas_circle.h
**
** Header file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_circle_Green(TOKAMAK *td, MEAS *m)
**
**				meas_circle_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_circle_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_circle.h
** Date:		March 24, 1993
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _MEAS_CIRCLE_

#define	_MEAS_CIRCLE_		1

#include "measurement.h"
#include "tokamak.h"


void	meas_circle_Green(TOKAMAK *, MEAS *);
void	meas_circle_Fit(TOKAMAK *, MEAS *);
void	meas_circle_Now(TOKAMAK *, MEAS *);
void	meas_circle_L(TOKAMAK *, MEAS *, double *);

#endif
