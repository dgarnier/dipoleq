/*
** TokaMac v2.0
**
** VAX_Alloc.h
**
**	This file makes some standard definitions in order to use the VMS
**	dynamic memory allocators.  These calls use the virtual memory checkers.
**
**
** File:		VAX_Alloc.h
** Date:		May 22, 1993
**
**
**
** (c) M. Mauel -- Columbia University
*/

#ifndef _VAX_ALLOC_

#define _VAX_ALLOC_	1

#ifdef VAXC

#define malloc	VAXC$MALLOC_OPT
#define calloc	VAXC$CALLOC_OPT
#define free	VAXC$FREE_OPT
#define cfree	VAXC$CFREE_OPT
#define realloc	VAXC$REALLOC_OPT

#endif

#ifdef IDL
#include <stdio.h>
#include "export.h"
#define malloc(x) IDL_MemAlloc(x,NULL,IDL_MSG_RET)
#define free(x)   IDL_MemFree(x,NULL,IDL_MSG_RET)
#endif



#endif
