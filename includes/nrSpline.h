/*
**	nrSpline.h
**
**	Cubic Spline Interpolation from Numerical Recipies.
**
**	Modified for double precision math by M. E. Mauel
**	Dept. of Applied Physics
**	May 4, 1993
**
**
*/

#ifndef _NRSPLINE_

#define _NRSPLINE_

#ifdef __cplusplus

extern "C" {

#endif /* __cplusplus */

void spline(double *, double *, int , double , double , double *);
void splint(double *, double *, double *, int , double , double *);
void splint_dervs(double *, double *, double *, int , double , double *, double *);

void splie2(double *, double *, double **, int , int , double **);
void splie2_alt(double *, double *, double **, int , int , double **);

void splin2(double *, double *, double **, double **, int , int ,
	double , double , double *);
void splin2_dervs(double *, double *, double **, double **, double **,
 				int , int , double , double , double *, double *, double *);

#ifdef __cplusplus

}

#endif /* __cplusplus */


#endif
