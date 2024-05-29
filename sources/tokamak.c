/*
** TokaMac v2.0
**
** tokamak.c
**
** This file contains routines which define the tokamak
** data structure.
**
** In addition, an event structure is defined which is used to
** pass action requests to the computational "engine".
**
** A user interface must be defined separately in order to
** pass data and instructions to the code.
**
** Since this is a data analysis code, we will use MKS units
** for all physical quantities.  (Some dimensionless variables
** are used as noted.)
**
** File:		tokamak.c
** Date:		February 27, 1993
**
** Revisions:
**
**		August 3, 1993		Added perfectly conducting shells
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdlib.h>
#include "VAX_Alloc.h"
#include <string.h>
#include <time.h>
#include "nrutil.h"
#include "coil.h"
#include "shell.h"
#include "green.h"
#include "limiter.h"
#include "separatrix.h"
#include "measurement.h"
#include "plasma.h"
#include "psigrid.h"
#include "tokgreen.h"
#include "tokamak.h"

TOKAMAK      *new_Tokamak()
{
	TOKAMAK      *td;
//	PSIGRID      *pg;
//	PLASMA       *pl;

	td = (TOKAMAK *) malloc((unsigned) sizeof(TOKAMAK));
	if (!td)
		nrerror("ERROR: Allocation error in new_Tokamak.");

	td->MaxIterFixed = 0;
	td->MaxIterFree = 0;
	td->IterFixed = 0;
	td->IterFree = 0;
	td->LHGreenStatus = LHGreenNotOK;
	td->MGreenStatus = MGreenNotOK;
	td->SGreenStatus = ShellGreenNotOK;
	td->SInductStatus = ShellInductNotOK;
	td->RestartStatus = RestartNotOK;
	td->RestartUnkns = RestartNotOK;
	td->NumEqualEq = 0;
	td->Confidence = 0.683;
	td->NumMCarloEq = 0;
	td->NumMCarloData = 0;
	td->MaxIterMCarlo = 14;

	strcpy(td->Name, " - ");
	strcpy(td->Info, " - ");
	strcpy(td->Iname, "TokIn.dat");
	strcpy(td->Oname, "TokOut");
	strcpy(td->LHname, "TokBndryGreens.bin");
	strcpy(td->MGname, "TokMeasGreens.bin");
	strcpy(td->SGname, "TokShellGreens.bin");
	strcpy(td->SMname, "TokShellInduct.bin");
	strcpy(td->RSname, "TokRestart.bin");
	strcpy(td->Start, " - ");
	strcpy(td->Stop, " - ");

	td->NumCoils = 0;
	td->NumShells = 0;
	td->NumSeps = 0;
	td->NumLimiters = 0;
	td->NumMeasures = 0;
	td->NumUnkns = 0;

	td->Plasma = new_Plasma();
	td->PsiGrid = new_PsiGrid();
	td->Coils = NULL;
	td->Shells = NULL;
	td->Limiters = NULL;
	td->Seps = NULL;
	td->Measures = NULL;
	td->LHPlasmaGreen = NULL;

	td->Covar = NULL;
	td->UnknVectors = NULL;
	td->SValues = NULL;

	return td;
}

/*
**
**	init_Tokamak
**
**	This subroutine allocates memory for the tokamak
**	data structure.
**
*/

void          init_Tokamak(TOKAMAK * td)
{
	int           i;

	init_Plasma(td->Plasma);
	init_PsiGrid(td->PsiGrid);

	td->Coils = (COIL **) malloc((unsigned) td->NumCoils * sizeof(COIL *));
	if (!td->Coils)
		nrerror("ERROR: Allocation error in init_Tokamak.");
	for (i = 0; i < td->NumCoils; i++)
		td->Coils[i] = new_Coil(0);

	td->Shells = (SHELL **) malloc((unsigned) td->NumShells * sizeof(SHELL *));
	if ((td->NumShells > 0) && !td->Shells)
		nrerror("ERROR: Allocation error in init_Tokamak.");
	for (i = 0; i < td->NumShells; i++)
		td->Shells[i] = new_Shell(0);

	td->Limiters = (LIMITER **) malloc((unsigned) td->NumLimiters * sizeof(LIMITER *));
	if (!td->Limiters)
		nrerror("ERROR: Allocation error in init_Tokamak.");
	for (i = 0; i < td->NumLimiters; i++)
		td->Limiters[i] = NULL;

	td->Seps = (SEPARATRIX **) malloc((unsigned) td->NumSeps * sizeof(SEPARATRIX *));
	if ((td->NumSeps > 0) && !td->Seps)
		nrerror("ERROR: Allocation error in init_Tokamak.");
	for (i = 0; i < td->NumSeps; i++)
		td->Seps[i] = NULL;

	td->Measures = (MEAS **) malloc((unsigned) td->NumMeasures * sizeof(MEAS *));
	if (!td->Measures)
		nrerror("ERROR: Allocation error in init_Tokamak.");
	for (i = 0; i < td->NumMeasures; i++)
		td->Measures[i] = NULL;
}

void          free_Tokamak(TOKAMAK * td, int full)
{
	int           i, num_subshells = 0;
	int           nmax = td->PsiGrid->Nsize;
	int           ncoils = td->NumCoils;
	int           nshells = td->NumShells;

	if (td->PsiGrid)
		free_PsiGrid(td->PsiGrid);

	if (td->Plasma)
		free_Plasma(td->Plasma);

	if (td->Coils)
		for (i = 0; i < ncoils; i++)
			if (td->Coils[i])
				free_Coil(td->Coils[i], nmax);

	if (td->Shells)
		for (i = 0; i < nshells; i++)
			if (td->Shells[i]) {
				num_subshells += td->Shells[i]->NumSubShells;
				free_Shell(td->Shells[i], nmax, ncoils);
			}
	if (td->Limiters)
		for (i = 0; i < td->NumLimiters; i++)
			if (td->Limiters[i])
				free_Limiter(td->Limiters[i]);

	if (td->Seps)
		for (i = 0; i < td->NumSeps; i++)
			if (td->Seps[i])
				free_Separatrix(td->Seps[i]);

	if (td->Measures)
		for (i = 0; i < td->NumMeasures; i++)
			if (td->Measures[i])
				free_Measure(td->Measures[i], nmax, ncoils, num_subshells);

	if (td->LHPlasmaGreen)
		free_LHary(td->LHPlasmaGreen, nmax);

	if (td->Covar)
		free_dmatrix(td->Covar, 1, td->NumUnkns, 1, td->NumUnkns);

	if (td->UnknVectors)
		free_dmatrix(td->UnknVectors, 1, td->NumUnkns, 1, td->NumUnkns);

	if (td->SValues)
		free_dvector(td->SValues, 1, td->NumUnkns);

	if (full && td)
		free(td);
}

void          SetStartTime(TOKAMAK * td)
{
	time_t        now;
	struct tm    *date;

	now = time(NULL);
	date = localtime(&now);

	strftime(td->Start, TIMEBUF, "%A, %b %d, %Y %I:%M:%S %p", date);
}

void          SetStopTime(TOKAMAK * td)
{
	time_t        now;
	struct tm    *date;

	now = time(NULL);
	date = localtime(&now);

	strftime(td->Stop, TIMEBUF, "%A, %b %d, %Y %I:%M:%S %p", date);
}
