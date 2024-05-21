/*
** TokaMac v2.0
**
** Separatrix header file.
** This files defines the basic data structure used
** to store tokamak separatrix.  A separatrix is
** basically a region within which to look for a
** saddle point of Psi.
**
** File:		include:separatrix.h
** Date:		February 12, 1993
**
** Example:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
**
** 12/7/98 DTG ... made some changes to allow separatrices
** to be more useful when computing for dipoles.
*/

#ifndef _SEPARATRIX_

#define _SEPARATRIX_ 1

#define	Separatrix_Off 		0
#define Separatrix_On		1

#define YesSeparatrix		1
#define NoSeparatrix		0

typedef struct separatrix {
    char			Name[32];
	int				Enabled;
	int				IsSeparatrix;	/* if non-zero, then there is a separartix */
	double 			X1,Z1;
	double			X2,Z2;
	double			XC,ZC;      /* DTG 12/7/98
	                              the position of the "center of the plasma from this sep */
	double			PsiSep;		/* the value of Psi at the separatrix */
	double			Xs,Zs;		/* the location of the separatrix */
	} SEPARATRIX;

#ifdef __cplusplus
extern "C" {
#endif

SEPARATRIX *new_Separatrix(void);
void free_Separatrix(SEPARATRIX *);

int	IsValidSeparatrix(SEPARATRIX *, double, double);

int IsPtDivertor(SEPARATRIX *, double , double , double , double );

#ifdef __cplusplus
}
#endif
#endif
