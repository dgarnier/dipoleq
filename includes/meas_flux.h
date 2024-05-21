/*
** TokaMac v2.0
**
** meas_flux.h
**
** Header file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_flux_Green(TOKAMAK *td, MEAS *m)
**
**				meas_flux_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_flux_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_flux.h
** Date:		March 24, 1993
**
**	USAGE:
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _MEAS_FLUX_

#define	_MEAS_FLUX_		1

#include "measurement.h"
#include "tokamak.h"


void	meas_flux_Green(TOKAMAK *, MEAS *);
void	meas_flux_Fit(TOKAMAK *, MEAS *);
void	meas_flux_Now(TOKAMAK *, MEAS *);
void	meas_flux_L(TOKAMAK *, MEAS *, double *);

#endif
