/*
**	nrRomberg.h
**
**	Romberg integration from Numerical Recipies.
**
**	Main routine:
**
**	integral = qromb(function, a, b);
**
**	where
**			function is defined as
**				double function(double x);
**
**			a, b are the limits of integration
**
**	Modified for double precision math by M. E. Mauel
**	Dept. of Applied Physics
**	May 4, 1993
**
*/

#ifndef _NRROMBERG_

#define _NRROMBERG_ 1

#ifdef __cplusplus
extern "C" {
#endif

double qromb(double (*)(double), double , double );
double qsimp(double (*)(double), double , double );

#ifdef __cplusplus
}
#endif

#endif
