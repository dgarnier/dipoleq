/*
** TokaMac v2.0
**
** Limiter header file.
** This files defines the basic data structure used
** to store tokamak limiters.  A limiter consists of
** line along which the outermost flux surface can
** be located.
**
** File:		include:limiter.h
** Date:		February 12, 1993
**
** Example:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _LIMITER_

#define _LIMITER_ 1

#define	Limiter_Off 	0
#define Limiter_On		1
#define Limiter_Inner	-1

typedef struct limiter {
	int					Enabled;	/* if non-zero, then use this limiter */
    char				Name[32];
	double 				X1,Z1;		/* The coordinates of the endpoints of */
	double				X2,Z2;		/* the limiter line.	*/
	double				PsiMin;		/* the minimum value of Psi along limiter */
	double				Xmin,Zmin;	/* the coordinates of the min of Psi along limiter */
	} LIMITER;

#ifdef __cplusplus
extern "C" {
#endif

LIMITER *new_Limiter();
void free_Limiter(LIMITER *);

#ifdef __cplusplus
}
#endif

#endif
