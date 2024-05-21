/*
** TokaMac v2.0
**
** meas_pnorm.h
**
** Header file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_pnorm_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_pnorm_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_pnorm.h
** Date:		November 15, 1993
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _MEAS_PNORM_

#define	_MEAS_PNORM_		1

#include "measurement.h"
#include "tokamak.h"

void	meas_pnorm_Fit(TOKAMAK *, MEAS *);
void	meas_pnorm_Now(TOKAMAK *, MEAS *);
void	meas_pnorm_L(TOKAMAK *, MEAS *, double *);

#endif
