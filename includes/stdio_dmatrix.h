/*
**	stdio_dmatrix.h
**
**	A utility to read and write dmatrix created
**	by nrutil.c
**
**
** File:		stdio_dmatrix.h
** Date:		March 24, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _STDIO_DMATRIX_

#define _STDIO_DMATRIX_		1

#include <stdio.h>

size_t fread_dmatrix(double **, long , long , long , long, FILE *);
size_t fwrite_dmatrix(double **, long , long , long , long, FILE *);

#endif
