/*
** TokaMac v2.0
**
** shell.c
**
** This file define basic global creation and destruction
** subroutines for a PERFECTLY conducting shell or conductor.
** These shells do NOT contain net current.
**
** File:		shell.c
** Date:		August 3, 1993
**
**
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdlib.h>
#include "VAX_Alloc.h"
#include <string.h>
#include <math.h>
#include "nrutil.h"
#include "shell.h"

#ifdef __cplusplus
extern "C" {
#endif

SUBSHELL      *copy_SubShell(SUBSHELL * sc);

#ifdef __cplusplus
}
#endif




/*
**
**	new_SubShell
**
*/

SUBSHELL     *new_SubShell(void)
{
	SUBSHELL     *sc;

	sc = (SUBSHELL *) malloc((unsigned) sizeof(SUBSHELL));
	if (!sc)
		nrerror("ERROR: Allocation error in new_SubShell");

	sc->X = 1.0;
	sc->Z = 0.0;
	sc->Radius = 0.001;
	sc->Current = 0.0;
	sc->ShellGreen = NULL;
	sc->PlasmaGreen = NULL;
	sc->CoilGreen = NULL;
	strcpy(sc->Name, " - ");

	return sc;
}

/*
**
**	free_SubShell
**
*/

void          free_SubShell(SUBSHELL * c, int nmax, int ncoils)
{

	if (c->ShellGreen)
		free_ShellGreen(c->ShellGreen, nmax);
	if (c->PlasmaGreen)
		free_dmatrix(c->PlasmaGreen, 0, nmax, 0, nmax);
	if (c->CoilGreen)
		free_dvector(c->CoilGreen, 0, ncoils - 1);

	free(c);
}

/*
**	Danger:  This doesn't really make a copy of the
**	the SubShell.  Instead, it copies the pointers
**	so that the Green's functions do not need to be reallocated.
**
*/
SUBSHELL      *copy_SubShell(SUBSHELL * sc)
{
	SUBSHELL      *new_sc;

	new_sc = new_SubShell();
	strcpy(new_sc->Name, sc->Name);
	new_sc->X = sc->X;
	new_sc->Z = sc->Z;
	new_sc->Radius = sc->Radius;
	new_sc->Current = sc->Current;
	new_sc->ShellGreen = sc->ShellGreen;	/* copy pointers only */
	new_sc->PlasmaGreen = sc->PlasmaGreen;	/* copy pointers only */
	new_sc->CoilGreen = sc->CoilGreen;		/* copy pointers only */

	return new_sc;
}

/*
**
**	add_SubShell
**
*/

void          add_SubShell(SHELL * c, SUBSHELL * sc)
{
	SUBSHELL    **newAry;
	int           is, scnum;

	scnum = c->NumSubShells;

	/*   realloc(..) appears to fail with VAXC !!
	newAry = (SUBSHELL **) realloc(c->SubShells, (unsigned) (c->NumSubShells + 1) * sizeof(SUBSHELL *));
	if (!newAry)
		nrerror("allocation error in add_SubShell");
	*/

	newAry = (SUBSHELL **) malloc((size_t) (scnum + 1) * sizeof(SUBSHELL *));
	if (newAry == NULL)
		nrerror("allocation error in add_SubShell");

	for (is = 0; is < scnum; is++) {
		newAry[is] = copy_SubShell(c->SubShells[is]);
		free(c->SubShells[is]);		/* do not free Green's functions */
	}
	free(c->SubShells);

	c->SubShells = newAry;
	c->SubShells[c->NumSubShells] = sc;
	c->NumSubShells += 1;
}

/*
**
** 	new_Shell
**
*/

SHELL        *new_Shell(int NumSubShells)
{
	SHELL        *c;
	int           isc;

	c = (SHELL *) malloc((unsigned) sizeof(SHELL));
	if (!c)
		nrerror("ERROR: Allocation error in new_Shell.");

	c->NumSubShells = NumSubShells;
	c->Enabled = Shell_Off;
	strcpy(c->Name, " - ");

	c->SubShells = (SUBSHELL **) malloc((unsigned) NumSubShells * sizeof(SUBSHELL *));
	if (!c->SubShells)
		nrerror("ERROR: Allocation error in new_Shell.");
	for (isc = 0; isc < NumSubShells; isc++)
		c->SubShells[isc] = new_SubShell();

	return c;
}

void          free_Shell(SHELL * c, int nmax, int ncoils)
{
	int           isc;

	for (isc = 0; isc < c->NumSubShells; isc++)
		if (c->SubShells[isc])
			free_SubShell(c->SubShells[isc], nmax, ncoils);
	free(c->SubShells);
	free(c);
}

SHELLGREEN   *new_ShellGreen(int nmax)
{
	SHELLGREEN   *cg;
	double       *gvec;

	/* alloc space for new lhvec */
	cg = (SHELLGREEN *) malloc((unsigned) sizeof(SHELLGREEN));
	if (!cg)
		nrerror("ERROR: Allocation failure in new_ShellGreen.");

	gvec = dvector(0, nmax);
	cg->Top = gvec;

	gvec = dvector(0, nmax);
	cg->Bot = gvec;

	gvec = dvector(0, nmax);
	cg->In = gvec;

	gvec = dvector(0, nmax);
	cg->Out = gvec;

	return cg;
}

void          free_ShellGreen(SHELLGREEN * cg, int nmax)
{
	free_dvector(cg->Top, 0, nmax);
	free_dvector(cg->Bot, 0, nmax);
	free_dvector(cg->In, 0, nmax);
	free_dvector(cg->Out, 0, nmax);

	free(cg);
}

/*
**	Self_Inductance
**
*/
double        Self_Inductance(SUBSHELL * s)
{
	return s->X * (log(8.0 * s->X / s->Radius) - 1.75);
}
