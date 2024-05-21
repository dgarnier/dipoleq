/*
** TokaMac v2.0
**
** J_IsoNoFlow.h
**
**
**
** File:		J_IsoNoFlow.h
** Date:		March 22, 1993
**
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _J_ISONOFLOW_

#define _J_ISONOFLOW_		1

#include "tokamak.h"

#ifdef __cplusplus
extern "C" {
#endif

void J_IsoNoFlow(TOKAMAK *,double **, double *, double *);
double J_IsoNoFlow_Loc(TOKAMAK * , int , int , double *, double *);

#ifdef __cplusplus
}
#endif

#endif
