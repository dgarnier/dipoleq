/*
**
** DipolEq v 1.0
** -- stolen completely from
**
** TokaMac v2.0
**
** TokaMac.c
**
** Main routine
**
** File:		TokaMac.c
** Date:		April 2, 1993
**
** Revisions:
**
**		August 5, 1993		Added perfectly conducting shells
**		October 18, 1993	Enable creation of "equal" equilibria
**
**
** Flow chart:
**
**	1.	Read input file.
**	2.	Initialize current either from file or parabolic profile.
**	3.	Cycle through iterations of (a) free-boundary equilibria, and
**		(b) leastsquares minimizations.
**	4.	After free-boundary equilibrium converges or after the given
**		number of iterations, stop cycle.
**	5.	Write out several output files.
**	6.	End.
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdlib.h>
#include "VAX_Alloc.h"
#include <string.h>
#include <stdio.h>
#ifdef __MWERKS__
#include <console.h>
#endif
#ifdef __SC__
#include <console.h>
#endif
#ifdef VAXC
#include <climsgdef.h>
#include <descrip.h>
#endif
#include "nrutil.h"
#include "psigrid.h"
#include "measurement.h"
#include "coil.h"
#include "tokamak.h"
#include "FileInput.h"
#include "FileOutput.h"
#include "FindJ.h"
#include "InitJ.h"
#include "AddCoilJ.h"
#include "AddShellJ.h"
#include "LoadBndryGreens.h"
#include "LoadMeasGreens.h"
#include "PlasmaBoundary.h"
#include "PsiBoundary.h"
#include "LeastSquares.h"
#include "GetPlasmaParameters.h"
#include "Restart.h"
#include "dUnkn.h"
#include "DelChiSqr.h"
#include "Find_ShellCurrent.h"
#include "FindMeasFit.h"

#define TRUE	1
#define FALSE	0

void          IterateSolution(TOKAMAK *, int *);
void          DoFreeBoundary(TOKAMAK *, int);
void          DoFixedBoundary(TOKAMAK *);

void          MakeEqualEq(TOKAMAK *);
void          MakeEqualEq2(TOKAMAK *);
#define MAKEEQUALEQ		MakeEqualEq

void          MakeMCarloEq(TOKAMAK *);
void          MakeMCarloData(TOKAMAK * );

#ifdef __cplusplus
#define EXTERN extern "C"
#else 
#define EXTERN extern
#endif

EXTERN FILE  *LogFile;
EXTERN long   GLOBAL_RAN3_SEED;	/* initilize to some negative number */

/*
**
**	M A I N   P R O G R A M
**
*/

main(int argc, char **argv)
{
	TOKAMAK      *td;
	char         *fn = "DipIn.dat";
	char         *lgn = "DipLog.out";
	int           i, IsFirst = TRUE;

	/* F I L E   I N P U T */

#ifdef __SC__
	argc = ccommand(&argv);
#endif
#ifdef __MWERKS__
	argc = ccommand(&argv);
#endif

	/*
	** VAX VMS does not have a standard command line interface
	**
	** For UNIX and other standard C systems, the C command-line
	** interface should be fine.
	**
	** tokamak -f infile -l logfile
	**
	** where all qualifiers are optional.
	*/
#ifndef VAXC
	if (argc > 2)
		for (i = 1; i < argc; i = i + 2) {
			if ((strcmp(argv[i], "-f") == 0) && (argv[i + 1]))
				fn = argv[i + 1];
			if ((strcmp(argv[i], "-l") == 0) && (argv[i + 1]))
				lgn = argv[i + 1];
		}
#else
	/*
	**	For VAX VMS systems, we use command line information (CLI)
	**	routines.  A typical input line might look like
	**
	**	$ TokaMac /in=TokIn.dat /log=tlog.out
	**
	**	The qualifiers are optional.
	*/

	int           cli_status;
	char          cli_value[255];
	short         cli_val_len;
	char          vax_fname[63];
	char          vax_lgname[63];

	struct dsc$descriptor_s cli_value_str =
	{255, 14, 1, cli_value};
	$DESCRIPTOR(cli_infile, "IN");
	$DESCRIPTOR(cli_logfile, "LOG");

	cli_status = cli$present(&cli_infile);	/* is /INFILE present ? */
	if (cli_status == CLI$_PRESENT) {
		cli_status = cli$get_value(&cli_infile, &cli_value_str, &cli_val_len);
		strncpy(vax_fname, cli_value, cli_val_len);
		fn = vax_fname;
	}
	cli_status = cli$present(&cli_logfile);	/* is /LOG present ? */
	if (cli_status == CLI$_PRESENT) {
		cli_status = cli$get_value(&cli_infile, &cli_value_str, &cli_val_len);
		strncpy(vax_lgname, cli_value, cli_val_len);
		lgn = vax_lgname;
	}
#endif

	LogFile = fopen(lgn, "w");
	if (!LogFile)
		nrerror("ERROR:	Could not open logfile for writting.");

	printf("\nTokaMac Version 2.0b3\n");
	fprintf(LogFile, "\nTokaMac Version 2.0b3\n");
	printf("TOKAMAC:	Starting new equilibrium fit using input from %s.\n", fn);
	fprintf(LogFile, "TOKAMAC:	Starting new equilibrium fit using input from %s.\n", fn);

	td = FileInput(fn);
	SetStartTime(td);
	InValuesOutput(td);

	printf("		Start time: %s\n", td->Start);
	fprintf(LogFile, "		Start time: %s\n", td->Start);

	/*  I N I T I A L I Z E   C U R R E N T */

	if (td->RestartStatus == RestartOK)
		ReadRestart(td->RSname, td);
	else
		InitJ(td->PsiGrid, td->Plasma);

	/*  F I N D    S O L U T I O N */

	IterateSolution(td, &IsFirst);

	/*  S A V E    R E S U L T S */

	WriteRestart(td->RSname, td);

	GetPlasmaParameters(td);
	SetStopTime(td);

	printf("INFO:	Writing output files.\n", td->Stop);
	fprintf(LogFile, "INFO:	Writing output files.\n", td->Stop);
	PsiGridOutput(td);
	ConductorsOutput(td);
	PlasmaOutput(td);
	FluxProfileOutput(td);
	MeasOutput(td);
	BndMomentsOutput(td);
	EQGRUMOutput(td);
	DCONOutput(td);

	if (td->NumEqualEq > 0)
		MAKEEQUALEQ(td);

	if (td->NumMCarloEq > 0)
		MakeMCarloEq(td);

	if (td->NumMCarloData > 0)
		MakeMCarloData(td);

	printf("\n		TokaMac Finished at %s.\n", td->Stop);
	fprintf(LogFile, "\n		TokaMac Finished at %s.\n", td->Stop);
	free_Tokamak(td, TRUE);

	fclose(LogFile);
}

/*
**	F R E E   B O U N D A R Y   S O L U T I O N
*/
void          DoFreeBoundary(TOKAMAK * td, int IsFirstTime)
{
	int           is;

	if (td->NumShells > 0)
		Find_ShellCurrent(td);

	LoadBndryGreens(td);
	PsiBoundary(td);
	free_BndryGreens(td);

	AddCoilJ(td);

	AddShellJ(td);

	GoPDE(td->PsiGrid);

	PlasmaBoundary(td);

	LoadMeasGreens(td);
	LeastSquares(td, IsFirstTime);
	free_MeasGreens(td);

	FindJ(td);
}

/*
**	F I X E D   B O U N D A R Y   S O L U T I O N
*/
void          DoFixedBoundary(TOKAMAK * td)
{
	int           is;

	AddCoilJ(td);

	AddShellJ(td);

	GoPDE(td->PsiGrid);

	PlasmaBoundary(td);

	LoadMeasGreens(td);
	LeastSquares(td, FALSE);
	free_MeasGreens(td);

	FindJ(td);
}

/*
**	I T E R A T E    S O L U T I O N
*/
void          IterateSolution(TOKAMAK * td, int *IsFirst)
{
	int           ifree;
	int           ifixed;

	printf("INFO:	Iterate solution with %d free and %d fixed iterations.\n",
		   td->MaxIterFree, td->MaxIterFixed);
	fprintf(LogFile, "INFO:	Iterate solution with %d free and %d fixed iterations.\n",
			td->MaxIterFree, td->MaxIterFixed);

	for (ifree = 1; ifree <= td->MaxIterFree; ifree++) {
		printf("\n		[Free iteration %d]\n", ifree);
		fprintf(LogFile, "\n		[Free iteration %d]\n", ifree);

		DoFreeBoundary(td, *IsFirst);
		for (ifixed = 1; ifixed <= td->MaxIterFixed; ifixed++) {
			printf("\n		[Fixed iteration %d]\n", ifixed);
			fprintf(LogFile, "\n		[Fixed iteration %d]\n", ifixed);
			DoFixedBoundary(td);
			td->IterFixed = ifixed;
		}
		td->IterFree = ifree;
		*IsFirst = FALSE;

		if ((td->PsiGrid->BoundError < td->PsiGrid->BoundThreshold) && (ifree != 1))
			break;
	}
	printf("\n\n");
	fprintf(LogFile, "\n\n");
}

/*
**	MakeEqualEq
**
*/

#define MAKE_ITERATIONS		5
#define MAKE_UNDERRELAX		0.6

void          MakeEqualEq(TOKAMAK * td)
{
	int           ie, i, it;
	double        del;
	double       *unkn0, *unkn1;
	char          Oname0[32], anam[4];

	strcpy(Oname0, td->Oname);

	GLOBAL_RAN3_SEED = -5;

	td->PsiGrid->UnderRelax1 = MAKE_UNDERRELAX;
	td->RestartUnkns = 1;

	del = DelChiSqr(td->Confidence, td->NumUnkns);

	printf("INFO:	Computing %d equal equilibria within dChiSqr = %g.\n",
		   td->NumEqualEq, del);
	fprintf(LogFile, "INFO:	Computing %d equal equilibria within dChiSqr = %g.\n",
			td->NumEqualEq, del);

	unkn0 = dvector(1, td->NumUnkns);
	unkn1 = dvector(1, td->NumUnkns);

	CopyUnknowns(td, unkn0);

	for (ie = 0; ie < td->NumEqualEq; ie++) {
		printf("\n		[Equal Eq %d]\n", ie + 1);
		fprintf(LogFile, "\n		[Equal Eq %d]\n", ie + 1);

		sprintf(anam, "_%d", ie + 1);
		strcpy(td->Oname, Oname0);
		strcat(td->Oname, anam);

		ReadRestart(td->RSname, td);

		dUnkn(td->UnknVectors, td->SValues, del, unkn0, unkn1, td->NumUnkns);
		RewriteUnknowns(td, unkn1);
		FindJ(td);

		/* Compute new equilibria */
		for (it = 0; it < MAKE_ITERATIONS; it++) {
			AddCoilJ(td);
			AddShellJ(td);
			GoPDE(td->PsiGrid);
			/*
			**	Readjust plasma boundary to new solution ??
			**	No. These variations are not correct when we NON-LINEARLY
			**	adjust the plasma boundary.  Then the scan in parameters are
			**	no longer valid.  That is PsiAxis and PsiLim change!
			**
			**	Yes. In order to compute the global parameters, we need
			**	to use our routines for flux mapping. They fail if we
			**	don't set the boundary correctly.
			*/
			PlasmaBoundary(td);
			FindJ(td);
		}

		/* F I N D   M E A S   F I T */
		LoadMeasGreens(td);
		FindMeasFit(td);
		free_MeasGreens(td);

		/* Compute new plasma parameters */
		GetPlasmaParameters(td);

		/* Write out remaining output parameters */
		printf("INFO:	Writing output files.\n", td->Stop);
		fprintf(LogFile, "INFO:	Writing output files.\n", td->Stop);
		PsiGridOutput(td);
		ConductorsOutput(td);
		PlasmaOutput(td);
		FluxProfileOutput(td);
		MeasOutput(td);
		BndMomentsOutput(td);
		EQGRUMOutput(td);
	}

	free_dvector(unkn1, 1, td->NumUnkns);
	free_dvector(unkn0, 1, td->NumUnkns);
}

#define PSIXMAX_EQ	0.95

void          MakeEqualEq2(TOKAMAK * td)
{
	int           ie, i, it;
	double        del;
	double       *unkn0, *unkn1;
	char          Oname0[32], anam[4];

	strcpy(Oname0, td->Oname);

	GLOBAL_RAN3_SEED = -5;

	/* td->Plasma->PsiXmax = PSIXMAX_EQ; */
	td->RestartUnkns = 1;

	del = DelChiSqr(td->Confidence, td->NumUnkns);

	printf("INFO:	Computing %d equal equilibria within dChiSqr = %g.\n",
		   td->NumEqualEq, del);
	fprintf(LogFile, "INFO:	Computing %d equal equilibria within dChiSqr = %g.\n",
			td->NumEqualEq, del);

	unkn0 = dvector(1, td->NumUnkns);
	unkn1 = dvector(1, td->NumUnkns);

	CopyUnknowns(td, unkn0);

	for (ie = 0; ie < td->NumEqualEq; ie++) {
		printf("\n		[Equal Eq %d]\n", ie + 1);
		fprintf(LogFile, "\n		[Equal Eq %d]\n", ie + 1);

		sprintf(anam, "_%d", ie + 1);
		strcpy(td->Oname, Oname0);
		strcat(td->Oname, anam);

		ReadRestart(td->RSname, td);

		dUnkn(td->UnknVectors, td->SValues, del, unkn0, unkn1, td->NumUnkns);
		RewriteUnknowns(td, unkn1);
		FindJ(td);

		/* Compute new equilibria */
		AddCoilJ(td);
		AddShellJ(td);
		GoPDE(td->PsiGrid);

		/* F I N D   M E A S   F I T */
		LoadMeasGreens(td);
		FindMeasFit(td);
		free_MeasGreens(td);

		/* Compute new plasma parameters */
		GetPlasmaParameters(td);

		/* Write out remaining output parameters */
		printf("INFO:	Writing output files.\n", td->Stop);
		fprintf(LogFile, "INFO:	Writing output files.\n", td->Stop);
		PsiGridOutput(td);
		ConductorsOutput(td);
		PlasmaOutput(td);
		FluxProfileOutput(td);
		MeasOutput(td);
		BndMomentsOutput(td);
		EQGRUMOutput(td);
	}

	free_dvector(unkn1, 1, td->NumUnkns);
	free_dvector(unkn0, 1, td->NumUnkns);
}

#undef PSIXMAX_EQ
#undef MAKE_ITERATIONS
#undef MAKE_UNDERRELAX

/*
**	MakeMCarloEq
**
*/

void          MakeMCarloEq(TOKAMAK * td)
{
	int           ie, i, it, IsFirst = FALSE;
	double        del;
	double       *unkn0, *unkn1;
	double		 *real_meas;	/* the actual measurements */
	MEAS		 *m;
	char          Oname0[32], anam[4];

	strcpy(Oname0, td->Oname);

	GLOBAL_RAN3_SEED = -5;

	td->RestartUnkns = 1;
	td->MaxIterFree = td->MaxIterMCarlo;

	/* C O P Y    A C T U A L    M E A S U R E M E N T S */
	real_meas = dvector(0,td->NumMeasures - 1);
	for (i = 0; i < td->NumMeasures; i++) {
		m = td->Measures[i];
		real_meas[i] = m->Value;
	}

	del = DelChiSqr(td->Confidence, td->NumUnkns);

	printf("INFO:	Computing %d Monte Carlo equilibria within dChiSqr = %g.\n",
		   td->NumMCarloEq, del);
	fprintf(LogFile, "INFO:	Computing %d Monte Carlo equilibria within dChiSqr = %g.\n",
			td->NumMCarloEq, del);

	unkn0 = dvector(1, td->NumUnkns);
	unkn1 = dvector(1, td->NumUnkns);

	CopyUnknowns(td, unkn0);

	for (ie = 0; ie < td->NumMCarloEq; ie++) {
		printf("\n		[MCarlo Eq %d]\n", ie + 1);
		fprintf(LogFile, "\n		[MCarlo Eq %d]\n", ie + 1);

		sprintf(anam, "_%d", ie + 1);
		strcpy(td->Oname, Oname0);
		strcat(td->Oname, anam);

		ReadRestart(td->RSname, td);

		dUnkn(td->UnknVectors, td->SValues, del, unkn0, unkn1, td->NumUnkns);
		RewriteUnknowns(td, unkn1);
		FindJ(td);

		/* Compute simulated Monte Carlo data set */
		LoadMeasGreens(td);
		FindMeasNow(td);
		free_MeasGreens(td);
		for (i=0; i<td->NumMeasures; i++) {
			m = td->Measures[i];
			m->Value = m->Now;
		}

		/* C O M P U T E   N E W   E Q U I L B R I A */
		IterateSolution(td, &IsFirst);

		/* Restore values */
		for (i = 0; i < td->NumMeasures; i++) {
			m = td->Measures[i];
			m->Value = real_meas[i];
		}

		/* Compute new plasma parameters */
		GetPlasmaParameters(td);

		/* Write out remaining output parameters */
		printf("INFO:	Writing output files.\n", td->Stop);
		fprintf(LogFile, "INFO:	Writing output files.\n", td->Stop);
		PsiGridOutput(td);
		ConductorsOutput(td);
		PlasmaOutput(td);
		FluxProfileOutput(td);
		MeasOutput(td);
		BndMomentsOutput(td);
		EQGRUMOutput(td);
	}

	free_dvector(unkn1, 1, td->NumUnkns);
	free_dvector(unkn0, 1, td->NumUnkns);
	free_dvector(real_meas, 0, td->NumMeasures);
}

/*
**	MakeMCarloData
**
*/

void          MakeMCarloData(TOKAMAK * td)
{
	int           ie, i, it, IsFirst = FALSE;
	double        del;
	double		 *real_meas;	/* the actual measurements */
	MEAS		 *m;
	char          Oname0[32], anam[4];

	strcpy(Oname0, td->Oname);

	GLOBAL_RAN3_SEED = -5;

	td->RestartUnkns = 1;
	td->MaxIterFree = td->MaxIterMCarlo;

	/* C O P Y    A C T U A L    M E A S U R E M E N T S */
	real_meas = dvector(0,td->NumMeasures - 1);
	for (i = 0; i < td->NumMeasures; i++) {
		m = td->Measures[i];
		real_meas[i] = m->Value;
	}

	printf("INFO:	Computing %d Monte Carlo Data.\n",td->NumMCarloData);
	fprintf(LogFile, "INFO:	Computing %d Monte Carlo Data.\n",td->NumMCarloData);

	for (ie = 0; ie < td->NumMCarloData; ie++) {
		printf("\n		[MCarlo Data %d]\n", ie + 1);
		fprintf(LogFile, "\n		[MCarlo Data %d]\n", ie + 1);

		sprintf(anam, "_%d", ie + 1);
		strcpy(td->Oname, Oname0);
		strcat(td->Oname, anam);

		ReadRestart(td->RSname, td);

		/* Compute simulated Monte Carlo data set */
		for (i=0; i<td->NumMeasures; i++) {
			m = td->Measures[i];
			m->Value = m->Value + gasdev(&GLOBAL_RAN3_SEED)*m->StdDev;
		}

		/* C O M P U T E   N E W   E Q U I L B R I A */
		IterateSolution(td, &IsFirst);

		/* Restore values */
		for (i = 0; i < td->NumMeasures; i++) {
			m = td->Measures[i];
			m->Value = real_meas[i];
		}

		/* Compute new plasma parameters */
		GetPlasmaParameters(td);

		/* Write out remaining output parameters */
		printf("INFO:	Writing output files.\n", td->Stop);
		fprintf(LogFile, "INFO:	Writing output files.\n", td->Stop);
		PsiGridOutput(td);
		ConductorsOutput(td);
		PlasmaOutput(td);
		FluxProfileOutput(td);
		MeasOutput(td);
		BndMomentsOutput(td);
		EQGRUMOutput(td);
	}

	free_dvector(real_meas, 0, td->NumMeasures);
}
