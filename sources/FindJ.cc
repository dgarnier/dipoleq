
/*
** TokaMac v2.0
**
** FindJ.c
**
**
**
** File:		FindJ.c
** Date:		March 25, 1993
**
** Revisions:
**
**		August 5, 1993		Eliminated all current outside plasma boundary.
**		August 6, 1993		Added FindJ_Loc
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "plasma.h"
#include "psigrid.h"
#include "tokamak.h"
#include "J_Std.h"
#include "J_IsoNoFlow.h"
#include "J_DipoleStd.h"
#include "FindJ.h"
#include "CPlasmaModel.h"
#include "multitask.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

extern "C" FILE  *LogFile;


/*
**  Z E R O J
**
** this is a fun little progam to simulate zero plasma pressure
** ... in other words no plasma currents...
**
*/
void					ZeroJ(TOKAMAK * td)
{
	PSIGRID *pg;
	double **Cur;
	int			ix, iz, nmax;

	pg = td->PsiGrid;
	Cur = pg->Current;
	nmax = pg->Nsize;

	for (ix = 0; ix < nmax; ix++)
		for (iz = 0; iz < nmax; iz++)
			Cur[ix][iz] = 0.0;

}

/*
**	F I N D J
**
**
*/
void          FindJ(TOKAMAK * td)
{
	PSIGRID      *pg;
	PLASMA       *pl;
	double      **J, **Cur;
	int           ix, iz, nmax;
	int         **ip;
	double        ur, s = 0.0;

	pg = td->PsiGrid;
	pl = td->Plasma;
	nmax = pg->Nsize;
	ur = pg->UnderRelax1;
	Cur = pg->Current;
	ip = pg->IsPlasma;

	J = dmatrix(0, nmax, 0, nmax);

	printf("INFO:	FindJ\n");
	fprintf(LogFile, "INFO:	FindJ\n");

	switch (pl->ModelType) {
	  case Plasma_Std:
		  if (!td->VacuumOnly) J_Std(td, J, pl->Pp[1], pl->G2p[1]);
		  break;
	  case Plasma_IsoNoFlow:
		  if (!td->VacuumOnly) J_IsoNoFlow(td, J, pl->Pp, pl->G2p);
		  break;
	  case Plasma_IsoFlow:

		  break;
	  case Plasma_AnisoNoFlow:

		  break;
	  case Plasma_AnisoFlow:

		  break;
	  default:
	  	pl->Model->FindJ(td,J);
	  	break;

	}
    MULTI;


	/*	U N D E R    R E L A X A T I O N */
	/* 		When fUnderRelax = 0, then no underrelaxation	*/
	/*		When fUnderRelax = 1, then no advancement of solution */
	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (ip[ix][iz] && !td->VacuumOnly)
				Cur[ix][iz] = (1.0 - ur) * J[ix][iz] + ur * Cur[ix][iz];
			else
				Cur[ix][iz] = 0.0;

	free_dmatrix(J, 0, nmax, 0, nmax);

    MULTI;

	/*	C H E C K   F I N A L   R E S U L T */
	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			s += Cur[ix][iz];
	s *= pg->dx * pg->dz;

	pl->Ip = s / MU0;

	printf("		[Ip = %g (A)]\n", pl->Ip);
	fprintf(LogFile, "		[Ip = %g (A)]\n", pl->Ip);
}

/*
**	F I N D J _ L O C
**
**
*/
double        FindJ_Loc(TOKAMAK * td, int ix, int iz)
{
	PLASMA       *pl;
	double        J = 0.0;

	pl = td->Plasma;

	switch (pl->ModelType) {
	  case Plasma_Std:
		  J = J_Std_Loc(td, ix, iz, pl->Pp[1], pl->G2p[1]);
		  break;
	  case Plasma_IsoNoFlow:
		  J = J_IsoNoFlow_Loc(td, ix, iz, pl->Pp, pl->G2p);
		  break;
	  case Plasma_IsoFlow:

		  break;
	  case Plasma_AnisoNoFlow:

		  break;
	  case Plasma_AnisoFlow:

		  break;
		default:
	  	J = pl->Model->FindJ_Loc(td, ix, iz);

	}

	return J;
}
