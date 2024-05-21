/*
** TokaMac v2.0
**
** meas_coilcur.h
**
** Header file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_coilcur_Green(TOKAMAK *td, MEAS *m)
**
**				meas_coilcur_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_coilcur_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_coilcur.h
** Date:		March 24, 1993
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _MEAS_COILCUR_

#define	_MEAS_COILCUR_		1

#include "measurement.h"
#include "tokamak.h"

void	meas_coilcur_Fit(TOKAMAK *, MEAS *);
void	meas_coilcur_Now(TOKAMAK *, MEAS *);
void	meas_coilcur_L(TOKAMAK *, MEAS *, double *);

#endif
