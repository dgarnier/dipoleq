/*
** TokaMac v2.0
**
** limiter.c
**
** This file define basic global creation and destruction
** subroutines.
**
** File:		limiter.c
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
#include "nrutil.h"
#include "limiter.h"

/*
**
** new_Limiter
**
*/

LIMITER      *new_Limiter(void)
{
	LIMITER      *lim;

	lim = (LIMITER *) malloc((unsigned) sizeof(LIMITER));
	if (!lim)
		nrerror("ERROR: Allocation error in new_Limiter.");

	lim->Enabled = 0;
	lim->X1 = 0.0;
	lim->Z1 = 0.0;
	lim->X2 = 0.0;
	lim->Z2 = 0.0;
	lim->PsiMin = 1.0;
	lim->Xmin = 0.0;
	lim->Zmin = 0.0;
	strcpy(lim->Name, "default");

	return lim;
}

/*
**
** free_Limiter
**
*/

void          free_Limiter(LIMITER * lim)
{
	if (lim)
		free(lim);
}
