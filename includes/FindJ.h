/*
** TokaMac v2.0
**
** FindJ.h
**
**
**
** File:		FindJ.h
** Date:		March 25, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/


#ifndef _FINDJ_

#define	_FINDJ_		1

#include "tokamak.h"

#ifdef __cplusplus
extern "C" {
#endif

void ZeroJ(TOKAMAK *);
void FindJ(TOKAMAK *);
double  FindJ_Loc(TOKAMAK * , int , int );

#ifdef __cplusplus
}
#endif

#endif
