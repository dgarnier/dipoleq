/*
** TokaMac v2.0
**
** LoadShellGreens.h
**
**
** File:		LoadShellGreens.h
** Date:		August 4, 1993
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _LOADSHELLGREENS_

#define _LOADSHELLGREENS_ 	1

#ifdef __cplusplus
extern "C" {
#endif

void 	LoadShellGreens(TOKAMAK * );
void    free_ShellGreens(TOKAMAK * );
void    RewriteShellGreens(TOKAMAK * );

#ifdef __cplusplus
}
#endif

#endif
