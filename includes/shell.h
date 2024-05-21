/*
** TokaMac v2.0
**
** Shell header file.
** This files defines the basic data structure used
** to store tokamak shell sets.  A shell set consists of several
** subshells representing the shell as a series of filaments.
**
** File:		include:shell.h
** Date:		August 3, 1993
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _SHELL_

#define _SHELL_ 		1

#define	Shell_Off 		0
#define Shell_On		1


/*
** SHELLGREEN
**
** This stores the LH Green's Function for a subshell.
** It contains four arrays for the top, bottom, inside, and
** outside of the computational grid.
**
*/

typedef struct shellgreen {
	double		*Top;
	double		*Bot;
	double		*In;
	double		*Out;
} SHELLGREEN;


/*
**
** SUBSHELL
**
** This is a particular winding of a subshell (filament).
**
*/

typedef struct subshell {
    char		Name[32];			/* the subshell name */
	double		X,Z;				/* the location of the subshell */
	double		Radius;				/* the (cross-sectional) radius of the subshell */
	double		Current;			/* the current within this filament */
	SHELLGREEN	*ShellGreen;		/* the Green Funcs-top, bot, inside, outside */
	double		**PlasmaGreen;		/* the Green's Funcs between plasma and subcoil */
	double		*CoilGreen;			/* the Green's Funcs between the coil sets and subcoil */
} SUBSHELL;


/*
**
** SHELL
**
** A shell is a collection of subshells (or filaments) which act as a
** perfectly conductor with ZERO net current.
**
** [Note: the subshell array has dimensions [0..NumSubShells-1]. ]
**
*/

typedef struct shell {
	int				Enabled;		/* if non-zero, then use this shell */
    char			Name[32];		/* the name of this shell set */
	int				NumSubShells;	/* the number of subshells */
	SUBSHELL		**SubShells;	/* an array of SUBSHELL pointers */
} SHELL;


/*
**
** Function Prototypes
**
*/

#ifdef __cplusplus
extern "C" {
#endif

SHELLGREEN *new_ShellGreen(int );
void free_ShellGreen(SHELLGREEN *,int );

SUBSHELL *new_SubShell(void);
void free_SubShell(SUBSHELL *,int , int);
double	Self_Inductance(SUBSHELL *);

SHELL *new_Shell(int );
void free_Shell(SHELL *,int , int);
void add_SubShell(SHELL *,SUBSHELL *);

#ifdef __cplusplus
}
#endif


#endif
