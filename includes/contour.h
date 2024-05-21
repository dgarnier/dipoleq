/***
 * Public Interfaces for Contour.c
 *
 * � M. E. Mauel -- Dept. of Applied Physics, Columbia University
 *               -- March 20, 1992
 ***/

#ifndef _CONTOUR_

#define _CONTOUR_				1

#ifdef __cplusplus
extern "C" {
#endif


#define CONTOUR_START			1
#define CONTOUR_STOP			2
#define CONTOUR_TRACE			3
#define CONTOUR_DONT_DRAW		4
#define CONTOUR_DRAW			5

#define CONTOUR_ALL				1
#define CONTOUR_ONLY_CLOSED 	2
#define CONTOUR_ONLY_OPEN		3

#define CONTOUR_NO_MIDPOINT		0
#define CONTOUR_MIDPOINT		1

 void contour(	double *,			/* x */
 				double *,			/* y */
 				double **,			/* z */
 				int,				/* ixmin */
 				int,				/* ixmax */
 				int,				/* iymin */
 				int,				/* iymax */
 				double,				/* height */
 				int,				/* CONTOUR_WHICH */
 				int,				/* CONTOUR_MIDPOINT */
 				void (*)(double,double,double,int)		/* a procedure */
 				);

/* DTG ... 2/2/99.  Added this procedure to do quicker single contours for modelling */

 void contour_point(	double *,			/* x */
 						double *,			/* y */
 						double **,			/* z */
 						int,				/* ixmin */
 						int,				/* ixmax */
 						int,				/* iymin */
 						int,				/* iymax */
 						double,				/* x to start contour */
 						double,             /* y to start contour */
 						int,				/* CONTOUR_MIDPOINT */
 						void (*)(double,double,double,int)		/* a procedure */
 					);
#ifdef __cplusplus
}    /* extern "C" */
#endif
#endif
