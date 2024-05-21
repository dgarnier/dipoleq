/*
** TokaMac v2.0
**
** FindMeasFit.c
**
** Finds the computed values for each measurement.
** This assumes that FindJ has been previously run in order
** to fill the Plasma and PsiGrid arrays.
**
** File:		FindMeasFit.c
** Date:		March 24, 1993
**
** Revisions:
**
**		August 6, 1993		Added xxx_Now
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#define DEBUG_MEASURES 1

#include "tokamak.h"
#include "measurement.h"
#include "FindMeasFit.h"
#if DEBUG_MEASURES
#include <stdio.h>
#endif


/*
**	FindMeasFit
**
*/
void          FindMeasFit(TOKAMAK * td)
{
	MEAS         *m;
	int           im;

	for (im = 0; im < td->NumMeasures; im++) {
		m = td->Measures[im];
#if DEBUG_MEASURES
		if (m->Now > 1e10) {
			fprintf(stderr,"OOPS!: ");
		}
#endif
		(*(m->FindFit)) (td, m);
#if DEBUG_MEASURES

		fprintf(stderr,"MeasValNowFit [ %s ] = %g %g %g\n",m->Name,m->Value,m->Now,m->Fit);
#endif
	}
}

/*
**	FindMeasNow
**
*/
void          FindMeasNow(TOKAMAK * td)
{
	MEAS         *m;
	int           im;

	for (im = 0; im < td->NumMeasures; im++) {
		m = td->Measures[im];
		(*(m->FindNow)) (td, m);
	}
}
