/*
** TokaMac v2.0
**
** AddCoilJ.h
**
** This file contains routines which add any coil currents
** located within the computational domain.
**
**
** File:		includes:AddCoilJ.h
** Date:		March 20, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _ADDCOILJ_

#define _ADDCOILJ_ 	1

#include "psigrid.h"
#include "coil.h"
#include "tokamak.h"

#ifdef __cplusplus
extern "C" {
#endif

void AddSubCoilJ(PSIGRID * , COIL * );
void AddCoilJ(TOKAMAK *);

#ifdef __cplusplus
}
#endif

#endif
