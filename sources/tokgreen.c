
/*
** TokaMac v2.0
**
** tokgreen.c
**
** This file contains routines which allocate and free
**	memory for the LHGreen functions.
**
** File:		tokgreen.c
** Date:		February 29, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <math.h>
#include <stdlib.h>
#include "VAX_Alloc.h"
#include "nrutil.h"
#include "tokgreen.h"

/***
 *
 * CREATION and DESTRUCTION of LHVEC and LHARY
 *
 * LHVEC --> 	[0,..,top,..nmax]
 *				[0,..,bot,..nmax]
 *      		[0,..,ins,..nmax]
 *				[0,..,out,,,nmax]
 * The amount of storage required for an LHVEC = 4 x (nmax+1).
 * For nmax = 64, size(LHVEC) = 260 * size(double).
 *
 * LHARY --> 	[[LHVEC_0],..,ins,..,[LHVEC_nmax/2]]
 *				[[LHVEC_0],..,out,..,[LHVEC_nmax/2]]
 *				[[LHVEC_0],..,bot,..,[LHVEC_nmax]  ]
 * The amount of storage required for an LHARY is VERY BIG
 * LHARY = 8 x (nmax/2+1) x (nmax+1) + 4 x (nmax+1)^2.
 * For nmax = 64, size(LHARY) = 34060 * size(double).
 *
 ***/

LHVEC        *new_LHvec(int nmax)
{
	LHVEC        *lhv;
	double       *gvec;

	/* alloc space for new LHVEC */
	lhv = (LHVEC *) malloc((unsigned) sizeof(LHVEC));
	if (!lhv)
		nrerror("ERROR: Allocation failure in new_LHvec.");

	gvec = dvector(0, nmax);
	lhv->Top = gvec;

	gvec = dvector(0, nmax);
	lhv->Bot = gvec;

	gvec = dvector(0, nmax);
	lhv->In = gvec;

	gvec = dvector(0, nmax);
	lhv->Out = gvec;

	return lhv;
}

void          free_LHvec(LHVEC * lhv, int nmax)
{
	free_dvector(lhv->Top, 0, nmax);
	free_dvector(lhv->Bot, 0, nmax);
	free_dvector(lhv->In, 0, nmax);
	free_dvector(lhv->Out, 0, nmax);

	free(lhv);
	lhv = NULL;
}

LHARY        *new_LHary(int nmax)
{
	LHARY        *lha;
	int           i;

	/* alloc space for new LHARY */
	lha = (LHARY *) malloc((unsigned) sizeof(LHARY));
	if (!lha)
		nrerror("ERROR: Allocation failure in new_LHary.");

	/* alloc space for inside array of LHVECs */
	lha->In = (LHVEC **) malloc((unsigned) (nmax / 2 + 1) * sizeof(LHVEC *));
	if (!lha->In)
		nrerror("ERROR: Allocation failure in new_LHary.");
	for (i = 0; i <= nmax / 2; i++)
		lha->In[i] = new_LHvec(nmax);

	/* alloc space for outside array of LHVECs */
	lha->Out = (LHVEC **) malloc((unsigned) (nmax / 2 + 1) * sizeof(LHVEC *));
	if (!lha->Out)
		nrerror("ERROR: Allocation failure in new_LHary.");
	for (i = 0; i <= nmax / 2; i++)
		lha->Out[i] = new_LHvec(nmax);

	/* alloc space for bottom array of LHVECs */
	lha->Bot = (LHVEC **) malloc((unsigned) (nmax + 1) * sizeof(LHVEC *));
	if (!lha->Bot)
		nrerror("ERROR: Allocation failure in new_LHary.");
	for (i = 0; i <= nmax; i++)
		lha->Bot[i] = new_LHvec(nmax);

	return lha;
}

void          free_LHary(LHARY * lha, int nmax)
{
	int           i;

	for (i = nmax / 2; i >= 0; i--)
		free_LHvec(lha->In[i], nmax);
	for (i = nmax / 2; i >= 0; i--)
		free_LHvec(lha->Out[i], nmax);
	for (i = nmax; i >= 0; i--)
		free_LHvec(lha->Bot[i], nmax);

	free(lha->In);
	free(lha->Out);
	free(lha->Bot);

	free(lha);
	lha = NULL;
}
