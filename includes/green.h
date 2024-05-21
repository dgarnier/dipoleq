/*
 * Green.h
 *
 * Interface file for Green.c
 *
 * Usage:
 *
 *	g = Green(x, z, xc, zc);   where (xc,zc) is the location of a coil
 *  GetdGreen(G,dGx,dGz, x,z, xc,zc);  where dGx == d(G)/dx
 *
 * (c) L. Bai & M. Mauel, Dept. Applied Physics, Columbia University
 * June 19, 1992
 */

#ifndef _GREEN_

#define _GREEN_ 	1

#ifdef __cplusplus
extern "C" {
#endif

double Green(double, double, double, double);
void GetdGreen(double *,double *,double *,double,double,double,double);

#ifdef __cplusplus
}
#endif


#endif
