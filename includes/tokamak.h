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
** Date:		April 2, 1993
**
** Revised:
**
**		August 3, 1993		Added perfectly conducting shells
**		Oct. 3, 1993		Added covariance matrix
**		Oct. 3, 1993		Added RestartUnkns
**		Oct. 18, 1993		Added UnknVectors
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _TOKAMAK_

#define _TOKAMAK_ 1

#include "coil.h"
#include "shell.h"
#include "green.h"
#include "limiter.h"
#include "separatrix.h"
#include "measurement.h"
#include "plasma.h"
#include "psigrid.h"
#include "tokgreen.h"

/*
**
**	Constants
**
*/

#define LHGreenOK							1
#define LHGreenNotOK						0

#define MGreenOK							1
#define MGreenNotOK							0

#define ShellGreenOK						1
#define ShellGreenNotOK						0

#define ShellInductOK						1
#define ShellInductNotOK					0

#define RestartOK							1
#define RestartNotOK						0

#define TIMEBUF								60


/*
**
**	TokaMacDoc Structure
**
**	All of the data needed to define a tokamak equilibrium.
**
*/

  typedef struct tokamak {
      int MaxIterFixed;		/* The max iterat for fixed-boundary loop */
      int MaxIterFree;		/* The max iterat for boundary loop */
      int IterFixed;		/* The fixed-boudary iteratation completed */
      int IterFree;			/* The outer, boundary iteration completed */
      int LHGreenStatus;	/* if non-zero, then LH Green Functions are found */
      int MGreenStatus;		/* if non-zero, then Meas Green's Functions are found */
      int SGreenStatus;		/* if non-zero, then Shell Green's Functions are found */
      int SInductStatus;	/* if non-zero, then Shell inductance matrix is found */
      int RestartStatus;	/* if non-zero, then initialize J and Psi from restart file */
      int RestartUnkns;		/* if non-zero, then also read prior unknowns */
      int VacuumOnly;       /* if non-zero, then ignore plasma contribution to J */
      int NumEqualEq;		/* if non-zero, then generates "equal" equilibria */
      int NumMCarloEq;		/* if non-zero, then generates "Monte Carlo" equilibria */
      int NumMCarloData;	/* if non-zero, then generates "Monte Carlo" data sets */
      int MaxIterMCarlo;	/* number of iterations for each Monte Carlo Iteration */

      char Name[32];		/* string identifier of reconstruction */
      char Info[32];		/* some description of the reconstruction */
      char Iname[32];		/* The input file name */
      char Oname[32];		/* A file name for output routines */
      char LHname[32];		/* A file name for saving the LH Green Functions */
      char MGname[32];		/* A file name for saving the measure Green's Functions */
      char SGname[32];		/* A file name for saving the Shell Green's Functions */
      char SMname[32];		/* A file name for saving the Shell inductance matrix */
      char RSname[32];		/* A file name for reading and writting the restart file */
	    char Start[TIMEBUF];	/* A string containing the start time */
	    char Stop[TIMEBUF];	/* A string containing the stop time */

      int NumCoils;
      int NumShells;
      int NumLimiters;
      int NumSeps;
      int NumMeasures;
      int NumUnkns;

      double **Covar;		/* the covariance matrix for the unknowns */
      double **UnknVectors;	/* the vectors of the unknowns from SVD */
      double *SValues;		/* the singular values */
      double Confidence;	/* the confidence level for EqualEq */

      PSIGRID *PsiGrid;		/* the Psi solution */
      PLASMA *Plasma;		/* the plasma global parameters */
      COIL **Coils;			/* the coil sets */
      SHELL **Shells;		/* the perfectly conducting shells */
      LIMITER **Limiters;	/* the limiters */
      SEPARATRIX **Seps;	/* the separatrices */
      MEAS **Measures;		/* pointer to an array of measurements */
      LHARY *LHPlasmaGreen;	/* the LH Green funcs for plasma current */

  } TOKAMAK;

/*
**
**	Public subroutines
**
*/

#ifdef __cplusplus
extern "C" {
#endif

TOKAMAK *new_Tokamak();
void init_Tokamak(TOKAMAK *);
void free_Tokamak(TOKAMAK *, int full);

void SetStartTime(TOKAMAK *);
void SetStopTime(TOKAMAK *);

#ifdef __cplusplus
} // extern "C"
#endif
#endif
