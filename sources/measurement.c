/*
** TokaMac v2.0
**
** measurement.c
**
** This file define basic global creation and destruction
** subroutines.
**
** A point of clarification is needed:
** meas->Fit, is the calculated measurement using the present
** value of the plasma current, etc.
** meas->Now, is the calculated measurement using the present
** unknowns.  These may be different since the present current
** density may not be equal to the present unknowns due
** to under-relaxation.
**
**
** File:		measurement.c
** Date:		February 7, 1993
**
** Revisions:
**
**		August 6, 1993		Added "Now"-->Fit with present unknowns
**
** Routine List:
**
**	=== CREATE A NEW MEASURE:
**		MEAS *new_Measure(int mtype)
**					mtype  === which type of measurement to make
**					nsize  === the size of psigrid (nsize x nsize)
**					ncoils === the number of coil sets
**
**	=== DESTROY A MEASURE:
**		void free_Measure(MEAS *aMeas,int nsize, int ncoils)
**					nsize  === the size of psigrid (nsize x nsize)
**					ncoils === the number of coil sets
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdlib.h>
#include "VAX_Alloc.h"
#include <string.h>
#include "nrutil.h"
#include "measurement.h"
#include "meas_bp.h"
#include "meas_bpangle.h"
#include "meas_plasmacur.h"
#include "meas_coilcur.h"
#include "meas_flux.h"
#include "meas_saddle.h"
#include "meas_circle.h"
#include "meas_press.h"
#include "meas_ppsix.h"
#include "meas_pnorm.h"
#include "meas_J0.h"

/*
**
**	meas_dJdy
**
*/

double     ***meas_dJdy;		/* storage pointer for least-squares */

/*
**
**	New_Measure
**
*/

MEAS         *new_Measure(int mtype)
{
	MEAS         *aMeas;

	aMeas = (MEAS *) malloc((unsigned) sizeof(MEAS));
	if (!aMeas)
		nrerror("ERRROR: Allocation error in new_Measure.");

	aMeas->mType = mtype;
	aMeas->Value = 0.0;
	aMeas->StdDev = 1.0;
	aMeas->Fit = 0.0;
	aMeas->Now = 0.0;
	strcpy(aMeas->Name, " - ");
	aMeas->FindGreen = NULL;
	aMeas->FindFit = NULL;
	aMeas->FindNow = NULL;
	aMeas->FindL = NULL;
	aMeas->X = 1.0;
	aMeas->Z = 0.0;
	aMeas->CoilGreen = NULL;
	aMeas->ShellGreen = NULL;
	aMeas->PlasmaGreen = NULL;

	switch (mtype) {
	  case meas_bp:
		  aMeas->parm.bp.Angle = 0.0;
		  aMeas->FindGreen = (void (*)(void *, void *)) meas_bp_Green;
		  aMeas->FindFit = (void (*)(void *, void *)) meas_bp_Fit;
		  aMeas->FindNow = (void (*)(void *, void *)) meas_bp_Now;
		  aMeas->FindL = (void (*)(void *, void *, double *)) meas_bp_L;
		  break;
	  case meas_flux:
		  aMeas->FindGreen = (void (*)(void *, void *)) meas_flux_Green;
		  aMeas->FindFit = (void (*)(void *, void *)) meas_flux_Fit;
		  aMeas->FindNow = (void (*)(void *, void *)) meas_flux_Now;
		  aMeas->FindL = (void (*)(void *, void *, double *)) meas_flux_L;
		  break;
	  case meas_saddle:
		  aMeas->parm.saddle.X1 = 1.0;
		  aMeas->parm.saddle.Z1 = 0.0;
		  aMeas->parm.saddle.X2 = 1.0;
		  aMeas->parm.saddle.Z2 = 1.0;
		  aMeas->FindGreen = (void (*)(void *, void *)) meas_saddle_Green;
		  aMeas->FindFit = (void (*)(void *, void *)) meas_saddle_Fit;
		  aMeas->FindNow = (void (*)(void *, void *)) meas_saddle_Now;
		  aMeas->FindL = (void (*)(void *, void *, double *)) meas_saddle_L;
		  break;
	  case meas_circle:
		  aMeas->parm.circle.Radius = 1.0;
		  aMeas->parm.circle.Number = 24;
		  aMeas->FindGreen = (void (*)(void *, void *)) meas_circle_Green;
		  aMeas->FindFit = (void (*)(void *, void *)) meas_circle_Fit;
		  aMeas->FindNow = (void (*)(void *, void *)) meas_circle_Now;
		  aMeas->FindL = (void (*)(void *, void *, double *)) meas_circle_L;
		  break;
	  case meas_coilcur:
		  aMeas->parm.coilcur.CoilNum = 0;
		  aMeas->FindFit = (void (*)(void *, void *)) meas_coilcur_Fit;
		  aMeas->FindNow = (void (*)(void *, void *)) meas_coilcur_Now;
		  aMeas->FindL = (void (*)(void *, void *, double *)) meas_coilcur_L;
		  break;
	  case meas_plasmacur:
		  aMeas->FindFit = (void (*)(void *, void *)) meas_plasmacur_Fit;
		  aMeas->FindNow = (void (*)(void *, void *)) meas_plasmacur_Now;
		  aMeas->FindL = (void (*)(void *, void *, double *)) meas_plasmacur_L;
		  break;
	  case meas_press:
		  aMeas->FindFit = (void (*)(void *, void *)) meas_press_Fit;
		  aMeas->FindNow = (void (*)(void *, void *)) meas_press_Now;
		  aMeas->FindL = (void (*)(void *, void *, double *)) meas_press_L;
		  break;
	  case meas_ppsix:
		  aMeas->FindFit = (void (*)(void *, void *)) meas_ppsix_Fit;
		  aMeas->FindNow = (void (*)(void *, void *)) meas_ppsix_Now;
		  aMeas->FindL = (void (*)(void *, void *, double *)) meas_ppsix_L;
		  break;
	  case meas_pnorm:
		  aMeas->FindFit = (void (*)(void *, void *)) meas_pnorm_Fit;
		  aMeas->FindNow = (void (*)(void *, void *)) meas_pnorm_Now;
		  aMeas->FindL = (void (*)(void *, void *, double *)) meas_pnorm_L;
		  break;
	  case meas_bpangle:
		  aMeas->FindGreen = (void (*)(void *, void *)) meas_bpangle_Green;
		  aMeas->FindFit = (void (*)(void *, void *)) meas_bpangle_Fit;
		  aMeas->FindNow = (void (*)(void *, void *)) meas_bpangle_Now;
		  aMeas->FindL = (void (*)(void *, void *, double *)) meas_bpangle_L;
		  break;
	  case meas_J0:
		  aMeas->FindFit = (void (*)(void *, void *)) meas_J0_Fit;
		  aMeas->FindNow = (void (*)(void *, void *)) meas_J0_Now;
		  aMeas->FindL = (void (*)(void *, void *, double *)) meas_J0_L;
		  break;
	  default:
	  	  nrinfo("Unknown measurement type!\n");
	}

	return aMeas;
}

void          free_Measure(MEAS * aMeas, int nsize, int ncoils, int nsubshells)
{
	if (aMeas->PlasmaGreen)
		free_dmatrix(aMeas->PlasmaGreen, 0, nsize, 0, nsize);
	if (aMeas->CoilGreen)
		free_dvector(aMeas->CoilGreen, 0, ncoils - 1);
	if (aMeas->ShellGreen)
		free_dvector(aMeas->ShellGreen, 0, nsubshells - 1);
	free(aMeas);
	aMeas = NULL;
}
