/*
** TokaMac v2.0
**
** L U Decomposition from Numerical Recipies, 2nd Edition.
**
** File:		include:ludcmp.h
** Date:		August 4, 1993
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _LUDCMP_

#define _LUDCMP_ 1

#ifdef __cplusplus

extern "C" {

#endif /* __cplusplus */


void ludcmp(double **, int, int *, double *);
void lubksb(double **, int, int *, double *);


#ifdef __cplusplus

}

#endif /* __cplusplus */


#endif
