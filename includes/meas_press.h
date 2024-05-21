/*
** TokaMac v2.0
**
** meas_press.h
**
** Header file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_press_Green(TOKAMAK *td, MEAS *m)
**
**				meas_press_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_press_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_press.h
** Date:		March 24, 1993
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _MEAS_PRESS_

#define	_MEAS_PRESS_		1

#include "measurement.h"
#include "tokamak.h"

void	meas_press_Fit(TOKAMAK *, MEAS *);
void	meas_press_Now(TOKAMAK *, MEAS *);
void	meas_press_L(TOKAMAK *, MEAS *, double *);

#endif
