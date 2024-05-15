/*
** TokaMac v2.0
**
** AddShellJ.h
**
** This file contains routines which add any shell currents
** located within the computational domain.
**
**
** File:		includes:AddShellJ.h
** Date:		August 3, 1993
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _ADDSHELLJ_

#define _ADDSHELLJ_ 	1

#include "psigrid.h"
#include "shell.h"
#include "tokamak.h"

#ifdef __cplusplus
extern "C" {
#endif

void AddSubShellJ(PSIGRID * , SHELL * );
void AddShellJ(TOKAMAK *);

#ifdef __cplusplus
}
#endif

#endif
