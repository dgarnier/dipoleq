/*
** TokaMac v2.0
**
** FindMeasFit.h
**
** Finds the computed values for each measurement.
** This assumes that FindJ has been previously run in order
** to fill the Plasma and PsiGrid arrays.
**
** File:		FindMeasFit.h
** Date:		March 24, 1993
**
**	USAGE:
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _FINDMEASFIT_

#define _FINDMEASFIT_		1

#include "tokamak.h"

#ifdef __cplusplus
extern "C" {
#endif

void FindMeasFit(TOKAMAK *td);
void FindMeasNow(TOKAMAK *td);

#ifdef __cplusplus
}
#endif

#endif
