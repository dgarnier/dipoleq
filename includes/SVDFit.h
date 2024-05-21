/*
**	TokaMac v2.0
**
**	SVDFit.h
**
**	Routines from Numerical Recipes
**	Used to compute the LeastSquares best fit
**  by applying Singular Value Decomposition.
**
**
**	File:		SVDFit.h
**	Date:		March 21, 1993
**
**	Routine list:
**
**
**
**	(c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _SVDFIT_

#define _SVDFIT_		1

#ifdef __cplusplus
extern "C" {
#endif

void svdcmp(double **,int ,int ,double *,double **);
void svbksb(double **, double *, double**, int ,int ,double *, double *);
void svdfit(double **,double **, double *,double **,double *,double *,int ,int ,double *);
void svdvar(double **, int , double *, double **);

#ifdef __cplusplus
}   // extern "C"
#endif

#endif
