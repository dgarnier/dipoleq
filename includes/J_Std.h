/*
** TokaMac v2.0
**
** J_Std.h
**
**
**
** File:		J_Std.h
** Date:		March 22, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _J_STD_

#define _J_STD_		1

#ifdef __cplusplus
extern "C" {
#endif

#include "tokamak.h"

void J_Std(TOKAMAK *,double **, double , double );
double J_Std_Loc(TOKAMAK * , int , int , double , double );
#ifdef __cplusplus
}    /* extern "C" */
#endif

#endif
