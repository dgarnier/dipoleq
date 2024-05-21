/*
** TokaMac v2.0
**
** meas_plasmacur.h
**
** Header file for local poloidal field measurements.
**
** Every measurement must define the following subroutines:
**
**				meas_plasmacur_Green(TOKAMAK *td, MEAS *m)
**
**				meas_plasmacur_Fit(TOKAMAK *td, MEAS *m)
**
**				meas_plasmacur_L(TOKAMAK *td, MEAS *m, double *L)
**
**
** File:		meas_plasmacur.h
** Date:		March 24, 1993
**
**	USAGE:
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _MEAS_PLASMACUR_

#define	_MEAS_PLASMACUR_		1

#include "measurement.h"
#include "tokamak.h"

void	meas_plasmacur_Fit(TOKAMAK *, MEAS *);
void	meas_plasmacur_Now(TOKAMAK *, MEAS *);
void	meas_plasmacur_L(TOKAMAK *, MEAS *, double *);

#endif
