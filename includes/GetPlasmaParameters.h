/*
** TokaMac v2.0
**
** GetPlasmaParameters.h
**
**
**
** File:		GetPlasmaParameters.h
** Date:		March 31, 1993
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _GETPLASMAPARAMETERS_

#define _GETPLASMAPARAMETERS_	1

#include "tokamak.h"

#ifdef __cplusplus
extern "C" {
#endif

void	GetGradPsi(TOKAMAK *);
void	GetPlasmaParameters(TOKAMAK *);

#ifdef __cplusplus
}
#endif
#endif
