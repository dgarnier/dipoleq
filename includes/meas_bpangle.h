/*
** TokaMac v2.0
**
** meas_bpangle.h
**
** Header file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_bpangle_Green(TOKAMAK *td, MEAS *m)
**
**				meas_bpangle_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_bpangle_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_bpangle.h
** Date:		November 15, 1993
**
**	USAGE:
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _MEAS_BPANGLE_

#define	_MEAS_BPANGLE_		1

#include "measurement.h"
#include "tokamak.h"


void	meas_bpangle_Green(TOKAMAK *, MEAS *);
void	meas_bpangle_Fit(TOKAMAK *, MEAS *);
void	meas_bpangle_Now(TOKAMAK *, MEAS *);
void	meas_bpangle_L(TOKAMAK *, MEAS *, double *);

#endif
