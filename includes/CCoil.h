/*
** TokaMac v2.0
**
** Coil header file.
** This files defines the basic data structure used
** to store tokamak coil sets.  A coil set consists of several
** coils (or winding packs).
**
** File:		include:coil.h
** Date:		February 9, 1993
**
** Example:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _COIL_

#define _COIL_ 1

#define	Coil_Off 	0
#define Coil_On		1


/*
** COILGREEN
**
** This stores the LH Green's Function for a coil set.
** It contains four arrays for the top, bottom, inside, and
** outside of the computational grid.
**
*/

typedef struct coilgreen {
	double		*Top;
	double		*Bot;
	double		*In;
	double		*Out;
} COILGREEN;


/*
**
** SUBCOIL
**
** This is a particular winding of a coil set.
**
*/

typedef struct subcoil {
    char		Name[32];				/* the subcoil name */
	double		X,Z;					/* the location of the subcoil */
	double		CurrentFraction;		/* the fraction of the total current */
} SUBCOIL;


/*
**
** COIL
**
** A coil set is a PF coil system with a particular measured current.
**
** [Note: the subcoil array has dimensions [0..NumSubCoils-1]. ]
**
*/

typedef struct coil {
	int				Enabled;		/* if non-zero, then use this coil */
	double			CoilCurrent;	/* the current flowing through the coil set */
    char			Name[32];		/* the name of this coil set */
	int				NumSubCoils;	/* the number of subcoils */
	SUBCOIL			**SubCoils;		/* an array of SUBCOIL pointers */
	COILGREEN		*CoilGreen;		/* the Green Funcs-top, bot, inside, outside */
} COIL;


/*
**
** Function Prototypes
**
*/

#ifdef __cplusplus
extern "C" {
#endif

COILGREEN *new_CoilGreen(int );
void free_CoilGreen(COILGREEN *,int );

SUBCOIL *new_SubCoil(void);

COIL *new_Coil(int );
void free_Coil(COIL *,int );
void add_SubCoil(COIL *,SUBCOIL *);

#ifdef __cplusplus
}
#endif

#endif
