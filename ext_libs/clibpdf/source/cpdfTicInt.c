/* cpdfTicInt.c -- determines recommended values for axis/domain min and max values
		   and tick/mesh values and intervals.
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

Standalone test compile:
cc -Wall -DMAINDEF -o TicInt cpdfTicInt.c

1999-08-17 [IO] 
	## MULTI_THREADING_OK ##

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "cpdflib.h"
#include "cpdf_private.h"

#ifndef _WIN32
#define		max(a, b)	((a)>(b)?(a):(b))
#endif

typedef struct {
	double LEcut;		/* less or equal cut off value */
	double minorInt;	/* minor interval */
	double majorInt;	/* major interval */
} CPDFticLU;

#ifdef CPDF_DEBUG
/* This function is already in cpdfUtil.c */
/* Convert a float into M and E, where M * 10 ^ E, and 1 <= M < 10 */
float getMantissaExp(float v, int *iexp)
{
int ie = 0;
float vv = fabs(v);
    if(v == 0.0) {
	*iexp = 0;
	return(0.0);
    }

    while( vv >= 10.0) {
	vv /= 10.0;
	ie++;
    }
    while( vv < 1.0) {
	vv *= 10.0;
	ie--;
    }
    *iexp = ie;
    if(v>=0.0) return(vv);
    else       return(-vv);
}

#endif


static CPDFticLU ticTable[] = {
/*	cutval, minor, major */
	{10.0,  1.0,  5.0},
	{12.0,  1.0,  4.0},
	{15.0,  1.0,  5.0},
	{16.0,  1.0,  4.0},
	{20.0,  1.0,  5.0},
	{25.0,  1.0,  5.0},
	{30.0,  2.0, 10.0},
	{35.0,  5.0, 10.0},
	{40.0,  5.0, 10.0},
	{50.0,  5.0, 10.0},
	{60.0,  10.0, 20.0},
	{80.0,  10.0, 20.0},
	{90.0,  10.0, 30.0},
	{100.0, 10.0, 50.0}
};

#define	NticTable	sizeof(ticTable)/sizeof(CPDFticLU)
/* #define NticTable	14 */



void cpdf_suggestLinearDomainParams(float vmin, float vmax, float *recmin, float *recmax,
		float *tic1ValMajor, float *intervalMajor,
		float *tic1ValMinor, float *intervalMinor)
{
double absvmin = fabs(vmin);
double absvmax = fabs(vmax);
double maxmax = max(absvmin, absvmax);
double absMantMax;
double revFac;
double drecmin, drecmax, dtic1ValMajor, dintervalMajor, dtic1ValMinor, dintervalMinor;
int expMax, idx=0;
int i;

/* Take care of cases like [3340 - 3944] -- subtract some common offset, e.g. 3000 */

	absMantMax = getMantissaExp(maxmax, &expMax);
	if(maxmax != 0.0) {
	    absMantMax *= 10.0;
	    expMax--;
	}
	absMantMax = ceil(absMantMax);
#ifdef CPDF_DEBUG
	printf("mantMax= %g, exp= %d\n", absMantMax, expMax);
#endif
	for(i=0; i< NticTable; i++) {
	    if(absMantMax <= ticTable[i].LEcut) {
		idx = i; break;
	    }
	}
	revFac = pow(10.0, (double)expMax);

#ifdef CPDF_DEBUG
	printf("suggested max= %g,  minorInt= %g, majorInt= %g, revFac= %g\n",
		ticTable[idx].LEcut*revFac,
		ticTable[idx].minorInt*revFac, ticTable[idx].majorInt*revFac, revFac);
#endif
	if(vmax > 0.0 && vmin >= 0.0) {			/* both positive */
	    drecmax = ticTable[idx].LEcut*revFac;		/* recommended max */
	    dintervalMinor = ticTable[idx].minorInt*revFac;
	    dintervalMajor = ticTable[idx].majorInt*revFac;
	    dtic1ValMajor = drecmax;
	    while( vmin < dtic1ValMajor)
		dtic1ValMajor -= dintervalMajor;
	    drecmin = dtic1ValMajor;
	    dtic1ValMinor = dtic1ValMajor;
	}
	else if(vmax >= 0.0 && vmin < 0.0) {
	    if(vmax > absvmin) {	/* positive side is larger */
		/* drecmax = ticTable[idx].LEcut*revFac; */
		dintervalMinor = ticTable[idx].minorInt*revFac;
		dintervalMajor = ticTable[idx].majorInt*revFac;
		drecmin = 0.0;			/* include 0 as tick */
		while( vmin < drecmin)
		    drecmin -= dintervalMajor;
		dtic1ValMinor = dtic1ValMajor = drecmin;
		drecmax = 0.0;
		while( vmax > drecmax)
		    drecmax += dintervalMajor;
	    }
	    else {			/* negative side is larger */
		/* drecmin = -ticTable[idx].LEcut*revFac; */
		dintervalMinor = ticTable[idx].minorInt*revFac;
		dintervalMajor = ticTable[idx].majorInt*revFac;
		drecmax = 0.0;
		while( vmax > drecmax)
		    drecmax += dintervalMajor;
		drecmin = 0.0;
		while( vmin < drecmin)
		    drecmin -= dintervalMajor;
		dtic1ValMinor = dtic1ValMajor = drecmin;
	    }
	}
	else {	/* both must be negative */
	    drecmin = -ticTable[idx].LEcut*revFac;
	    dintervalMinor = ticTable[idx].minorInt*revFac;
	    dintervalMajor = ticTable[idx].majorInt*revFac;

	    dtic1ValMinor = dtic1ValMajor = drecmin;
	    drecmax = drecmin;
	    while( vmax > drecmax)
		drecmax += dintervalMajor;
	}
	*recmin = drecmin;
	*recmax = drecmax;
	*tic1ValMajor = dtic1ValMajor;
	*intervalMajor = dintervalMajor;
	*tic1ValMinor = dtic1ValMinor;
	*intervalMinor = dintervalMinor;
}


void cpdf_suggestMinMaxForLinearDomain(float vmin, float vmax, float *recmin, float *recmax)
{
float a, b, c, d;	/* unused return values for this call */
	cpdf_suggestLinearDomainParams(vmin, vmax, recmin, recmax, &a, &b, &c, &d);
}


#ifdef MAINDEF

void main(void)
{
float vmax, vmin;
float recmax=1.0, recmin=0.0;
float tic1Major=0.0, intvMajor=1.0;
float tic1Minor=0.0, intvMinor=1.0;
/*
	printf("Size of table: %d\n", NticTable);
	for(i=0; i<NticTable; i++) {
	    printf("%g %g %g\n", ticTable[i].LEcut, ticTable[i].minorInt, ticTable[i].majorInt);
	}
*/

again:
	printf("Min: "); fflush(stdout);
	scanf("%f", &vmin);
	printf("Max: "); fflush(stdout);
	scanf("%f", &vmax);
	cpdf_suggestLinearDomainParams(vmin, vmax, &recmin, &recmax,
		&tic1Major, &intvMajor, &tic1Minor, &intvMinor);

	printf("Axis (%g .. %g)\n", recmin, recmax);
	printf("Major ticks:   first= %g,  interval= %g\n", tic1Major, intvMajor);
	printf("Minor ticks:   first= %g,  interval= %g\n\n", tic1Minor, intvMinor);
	goto again;
}

#endif

