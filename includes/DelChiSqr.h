/*
**	DelChiSqr.h
**
**	This routine returns the value of Delta-ChiSqr
**	corresponding to a confidence level of "p" and
**	a degree of freedom, "nu".
**
**	From Numerical Recipies in C, 2nd Ed., p. 697.
**
**	M. E. Mauel
**	Columbia University
**
**	October 18, 1993
**
*/

#ifndef _DELCHISQR_

#define _DELCHISQR_	1

#ifdef __cplusplus
extern "C" {
#endif

double  gammq(double , double );
double	DelChiSqr(double , int );

#ifdef __cplusplus
}
#endif
#endif
