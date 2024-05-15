/*
** TokaMac v2.0
**
** InitJ.h
**
** This file contains routines which initializes the plasma
** current using info within the Plasma data structure.
**
**
** File:		includes:InitJ.h
** Date:		March 20, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _INITJ_

#define _INITJ_ 	1

#include "psigrid.h"
#include "plasma.h"

#ifdef __cplusplus
extern "C" {
#endif

void InitJ(PSIGRID * , PLASMA * );

#ifdef __cplusplus
}
#endif

#endif
