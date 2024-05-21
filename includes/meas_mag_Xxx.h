/*
** TokaMac v2.0
**
** meas_mag_Fit.h
**
** Source file for local poloidal field measurements.
**
** Every magnetic measurement has a common subroutine:
**
**				meas_mag_Fit(TOKAMAK *td, MEAS *m)
**
**
** File:		meas_mag_Fit.h
** Date:		March 24, 1993
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _MEAS_MAG_XXX_

#define _MEAS_MAG_XXX_		1

#include "measurement.h"
#include "tokamak.h"

void meas_mag_Fit(TOKAMAK *, MEAS *);
void meas_mag_Now(TOKAMAK *, MEAS *);
void meas_mag_L(TOKAMAK * , MEAS * , double *);


#endif
