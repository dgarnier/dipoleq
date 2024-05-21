/*
** TokaMac v2.0
**
** FileOutput.h
**
** This file defines subroutines to write
** various output quantities.
**
** File:		FileOutput.h
** Date:		March 20, 1993
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _FILEOUTPUT_

#define _FILEOUTPUT_ 1

#include "tokamak.h"

#ifdef __cplusplus
extern "C" {
#endif

void InValuesOutput(TOKAMAK *);
void ConductorsOutput(TOKAMAK * );
void PsiGridOutput(TOKAMAK *);
void PlasmaOutput(TOKAMAK *);
void MeasOutput(TOKAMAK *);
void FluxProfileOutput(TOKAMAK *);
void BndMomentsOutput(TOKAMAK *);
void DCONOutput(TOKAMAK * td);
void EQGRUMOutput(TOKAMAK *);
void GS2Output(TOKAMAK * td);

#if NIMRODOUTPUT

void NimrodOutput(TOKAMAK *td);

#endif

#ifdef __cplusplus
}
#endif

#endif
