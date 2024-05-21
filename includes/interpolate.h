/*
** TokaMac v2.0
**
** interpolate.h
**
**
** File:		interpolate.h
** Date:		March 12, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _INTERPOLATE_

#define _INTERPOLATE_	1

#include "psigrid.h"

#ifdef __cplusplus
extern "C" {
#endif

double  interpolate(PSIGRID * , double **, double , double );
double  interpolate_int(PSIGRID * , int **, double , double );

#ifdef __cplusplus
}
#endif

#endif
