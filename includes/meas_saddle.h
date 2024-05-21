/*
** TokaMac v2.0
**
** meas_saddle.h
**
** Header file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_saddle_Green(TOKAMAK *td, MEAS *m)
**
**				meas_saddle_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_saddle_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_saddle.h
** Date:		March 24, 1993
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _MEAS_SADDLE_

#define	_MEAS_SADDLE_		1

#include "measurement.h"
#include "tokamak.h"


void	meas_saddle_Green(TOKAMAK *, MEAS *);
void	meas_saddle_Fit(TOKAMAK *, MEAS *);
void	meas_saddle_Now(TOKAMAK *, MEAS *);
void	meas_saddle_L(TOKAMAK *, MEAS *, double *);

#endif
