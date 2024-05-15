/*
** TokaMac v2.0
**
** PsiBoundary.c
**
** Computes the value of psi around the computational domain
** for an axisymmetric tokamak.
**
** The plasma current and a coil set are computed separately.
**
** NOTE: You must run PlasmaBoundary BEFORE running CoilBoundary.
**
**
** File:		PsiBoundary.c
** Date:		March 20, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _PSIBOUNDARY_

#define _PSIBOUNDARY_ 			1

#include "tokgreen.h"
#include "coil.h"
#include "psigrid.h"
#include "tokamak.h"

#if defined(__cplusplus)
extern "C" {
#endif

void PsiPlasmaBoundary(LHARY *,PSIGRID *);
void PsiCoilBoundary(PSIGRID *, COIL *);
void PsiBoundary(TOKAMAK *);

#if defined(__cplusplus)
}
#endif

#endif
