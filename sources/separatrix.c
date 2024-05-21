
/*
** TokaMac v2.0
**
** separatrix.c
**
** This file define basic global creation and destruction
** subroutines.
**
** File:		separatrix.c
** Date:		February 12, 1993
**
** Routine list:
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
#include "separatrix.h"

#define SQUARE(x)	((x) * (x))

SEPARATRIX   *new_Separatrix(void)
{
	SEPARATRIX   *s;

	s = (SEPARATRIX *) malloc((unsigned) sizeof(SEPARATRIX));
	if (!s)
		nrerror("ERROR: 	Allocation error in new_Separatrix.");

	s->Enabled = Separatrix_Off;
	s->IsSeparatrix = NoSeparatrix;
	s->X1 = 1.1;
	s->Z1 = -0.4;
	s->X2 = 1.1;
	s->Z2 = 0.4;
	s->PsiSep = 0.0;
	s->Xs = 0.0;
	s->Zs = 0.0;
	strcpy(s->Name, "default");

	return s;
}

void          free_Separatrix(SEPARATRIX * s)
{
	if (s)
		free(s);
}

/*
**	IsValidSeparatrix
**
**	Returns 0 if (x,z) is outside valid region within which
** 	we can look for a separatrix.
**
**	Returns 1 if (x,z) is within valid region.
**
*/
int           IsValidSeparatrix(SEPARATRIX * s, double x, double z)
{
	double        rs, xc, zc;

	xc = 0.5 * (s->X1 + s->X2);
	zc = 0.5 * (s->Z1 + s->Z2);
	rs = sqrt(SQUARE(s->X2 - s->X1) + SQUARE(s->Z2 - s->Z1)) / 2.0;

	if (rs > sqrt(SQUARE(x - xc) + SQUARE(z - zc)))
		return 1;
	else
		return 0;
}

/*
**
** IsPtDivertor
**
** (x,z)   is the point in question.
** (xa,za) is the magnetic axis.
**
** Basically, we excluse the entire half-space that is bounded by a line
** perpendicular to a line segment passing from the magnetic axis (xa,za)
** to the separatix.
**
** We calculate this by defining a vector from (xa, za) to (xs, zs)...
**		a == (xs - xa, zs - za)
** and a vector from (xa,za) to (x,z)....
**		b == ( x - xa, z - za)
** Then, if a.b >= a.a, then (x,z) must be within a divertor.
**
**
** Should we only erase within the sep's domain? Answer: No.
**
** DTG 12/7/98
**
** changed this algorithm to use a special "center" for better results
** with dipoles.
**
*/
int           IsPtDivertor(SEPARATRIX * s, double x, double z, double xa, double za)
{
	double        a_b, a_a;

	xa = s->XC;
	za = s->ZC;

	if (s->Enabled) {
		a_b = (s->Xs - xa) * (x - xa) + (s->Zs - za) * (z - za);
		a_a = SQUARE(s->Xs - xa) + SQUARE(s->Zs - za);
		return (a_b >= a_a);
	} else
		return 0;
}
