/*
** TokaMac v2.0
**
** Restart.h
**
** Routines to read and write the restart file.
** Only the PsiGrid->Current and PsiGrid->Psi arrays
** are writtten.
**
** File:		Restart.h
** Date:		April 2, 1993
**
** Modifications:
**
**	Oct. 3, 1993		Restart solution as well as PsiGrid
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _RESTART_

#define _RESTART_	1

#include "tokamak.h"

#ifdef __cplusplus
extern "C" {
#endif

void WriteRestart(char *, TOKAMAK *);
void ReadRestart(char *, TOKAMAK *);

#ifdef __cplusplus
}   // extern "C"
#endif
#endif
