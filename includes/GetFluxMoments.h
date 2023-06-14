/*
** TokaMac v2.0
**
** GetFluxMoments.h
**
**
**
** File:		GetFluxMoments.h
** Date:		April 10, 1993
**
** Routine list:
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _GETFLUXMOMENTS_

#define _GETFLUXMOMENTS_	1

#ifdef __cplusplus
extern "C" {
#endif

#include "psigrid.h"

void Trace_Count(double , double , double , int );
void FindTheta(void);
void GetFluxContour(PSIGRID *, double , double **, double **, int *);
void GetFluxMoments(PSIGRID * , double , double *, double *, int );

#ifdef __cplusplus
}
#endif

#endif
