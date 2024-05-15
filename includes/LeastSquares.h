/*
**	TokaMac v2.0
**
**	LeastSquares.h
**
**	This file  contains the least-squares fitting routines which
**	adjust the plasma parameters in order to fit the
**	experimental measurements.
**
**
**	File:		LeastSquares.h
**	Date:		March 20, 1993
**
**	Modifcations:
**
**		Oct. 18, 1993		Created Copy & Rewrite utilities
**
**
**	(c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _LEASTSQUARES_

#define _LEASTSQUARES_		1

#include "tokamak.h"

#ifdef __cplusplus
extern "C" {
#endif

void CopyUnknowns(TOKAMAK *, double *);
void RewriteUnknowns(TOKAMAK *, double *);
void LeastSquares(TOKAMAK *,int );

#ifdef __cplusplus
}
#endif

#endif
