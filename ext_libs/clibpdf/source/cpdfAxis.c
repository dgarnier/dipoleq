/* cpdfAxis.c  -- axis at any orientation (horizontal, vertical, and oblique)
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf.
 *  or http://www.fastio.com/LICENSE.pdf

## MULTI_THREADING_OK ##

1998-08-24 [IO]
*/

#include "config.h"
#include "version.h"

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cpdflib.h"		/* This must be included before all other local include files */
#include "cpdf_private.h"
#include "cglobals.h"

typedef struct {
	float Tcut;
	int   MinorVar;
	int   MajorVar;
	int   MinorBump;
	int   MajorBump;
} CPDFtimeLU;



/* ------------------------------------------------------------------------------------------ */

/* A plot domain must exist before calling this function. */
/* ## MULTI_THREADING_OK ## */
CPDFaxis *cpdf_createAxis(float angle, float axislength, int typeflag, float valL, float valH)
{
CPDFaxis *anAx;
float sL, sH;
	anAx = (CPDFaxis *) malloc((size_t)sizeof(CPDFaxis));
	_cpdf_malloc_check((void *)anAx);
	/* It will be attached later by cpdf_attachAxisToDomain() to a real plot domain.
	   But, for now just attach it to the default domain.
	*/
	anAx->magic = (unsigned long)AXIS_MAGIC_NUMBER;
	/* anAx->plotDomain = pdf->defaultDomain; */
	anAx->plotDomain = NULL;		/* Intially attached to no domain */
	anAx->angle = angle;
	anAx->type = typeflag;
	anAx->xloc = 0.0;
	anAx->yloc = 0.0;
	anAx->length = axislength;
	if(typeflag == LINEAR ) {
	    cpdf_suggestLinearDomainParams(valL, valH, &sL, &sH,
		&anAx->valFirstTicLinMajor, &anAx->ticIntervalLinMajor,
		&anAx->valFirstTicLinMinor, &anAx->ticIntervalLinMinor);
	}
	anAx->valL = valL;
	anAx->valH = valH;

	anAx->ticEnableMajor = 1;
	anAx->ticEnableMinor = 1;
	anAx->ticLenMajor = 10.0;
	anAx->ticLenMinor = 5.0;
	anAx->axisLineWidth = 1.0;
	anAx->tickWidthMajor = 2.0;
	anAx->tickWidthMinor = 1.0;
	if(angle > 60.0 && angle < 120.0) {
	    anAx->ticPosition = 0;	/* almost vertical axis */
	    anAx->numPosition = 2;
	}
	else {
	    anAx->ticPosition = 2;
	    anAx->numPosition = 0;
	}

	anAx->numEnable = 1;		/* put numbers */
	anAx->ticNumGap = 7.0;
	anAx->numFontSize = 14.0;
	anAx->horizNumber = 1;
	anAx->numStyle = 0;	 /* 0=Regular, 1=exponential, 2=free style (e.g. Jan98  Feb98 ..) */
	anAx->numFormat = NULL;
	anAx->numFontName = NULL;
	anAx->numEncoding = NULL;

	anAx->numLabelGap = 7.0; /* gap (in points) between number and axis label */
	anAx->labelFontSize = 18.0;
	anAx->horizLabel = 0;	 /* if true, label will be plotted at horizontal orientation */
	anAx->labelEncoding = NULL;
	anAx->labelFontName = NULL;
	anAx->axisLabel = NULL;  /* Axis label string, if NULL, no label is shown */
	/* ------------------------------------------------------------ */
	/* num   min   1    2    3    4    5    6    7    8    9   max  */
	/* bit    0    1    2    3    4    5    6    7    8    9    10  */
	/* ------------------------------------------------------------ */
	anAx->numSelectorLog = LOGAXSEL_1;
	anAx->ticSelectorLog = LOGAXSEL_123456789;

	cpdf_setAxisNumberFormat(anAx, "%g", "Helvetica", 12.0);

	return anAx;
}

/* ## MULTI_THREADING_OK ## */
/* Create Time Axis */
CPDFaxis *cpdf_createTimeAxis(float angle, float axislength, int typeflag, struct tm *vTL, struct tm *vTH)
{
CPDFaxis *anAx;
	anAx = (CPDFaxis *) malloc((size_t)sizeof(CPDFaxis));
	_cpdf_malloc_check((void *)anAx);
	/* It will be attached later by cpdf_attachAxisToDomain() to a real plot domain.
	   But, for now just attach it to the default domain.
	*/
	anAx->magic = (unsigned long)AXIS_MAGIC_NUMBER;
	/* anAx->plotDomain = defaultDomain; */
	anAx->plotDomain = NULL;
	anAx->angle = angle;
	anAx->type = typeflag;
	anAx->xloc = 0.0;
	anAx->yloc = 0.0;
	anAx->length = axislength;
	cpdf_mktime(vTL); cpdf_mktime(vTH);	/* This sets all fields of struct tm correctly */
	memcpy(&anAx->vTL, vTL, sizeof(struct tm));
	memcpy(&anAx->vTH, vTH, sizeof(struct tm));
	anAx->valL = 0.0;			/* keep also float representation for date */
	anAx->valH = tm_to_NumDays(vTL, vTH);	 /* in float number of days */
	anAx->ticEnableMajor = 1;
	anAx->ticEnableMinor = 1;
	anAx->ticLenMajor = 10.0;
	anAx->ticLenMinor = 5.0;
	anAx->axisLineWidth = 1.0;
	anAx->tickWidthMajor = 2.0;
	anAx->tickWidthMinor = 1.0;
	if(angle > 60.0 && angle < 120.0) {
	    anAx->ticPosition = 0;	/* almost vertical axis */
	    anAx->numPosition = 2;
	}
	else {
	    anAx->ticPosition = 2;
	    anAx->numPosition = 0;
	}

	anAx->numEnable = 1;		/* put numbers */
	anAx->ticNumGap = 7.0;
	anAx->numFontSize = 14.0;
	anAx->useMonthName = 1;		/* use month names */
	anAx->use2DigitYear = 0;	/* use full year number as in 1998 */
	anAx->horizNumber = 1;
	anAx->numStyle = 0;	 /* 0=Regular, 1=exponential, 2=free style (e.g. Jan98  Feb98 ..) */
	anAx->numFormat = NULL;
	anAx->numFontName = NULL;
	anAx->numEncoding = NULL;

	anAx->numLabelGap = 7.0; /* gap (in points) between number and axis label */
	anAx->labelFontSize = 18.0;
	anAx->horizLabel = 0;	 /* if true, label will be plotted at horizontal orientation */
	anAx->labelEncoding = NULL;
	anAx->labelFontName = NULL;
	anAx->axisLabel = NULL;  /* Axis label string, if NULL, no label is shown */
	/* ------------------------------------------------------------ */
	/* num   min   1    2    3    4    5    6    7    8    9   max  */
	/* bit    0    1    2    3    4    5    6    7    8    9    10  */
	/* ------------------------------------------------------------ */
	anAx->numSelectorLog = LOGAXSEL_1;
	anAx->ticSelectorLog = LOGAXSEL_123456789;

	anAx->lastMin=99;
	anAx->lastHour=99;
	anAx->lastDay=99;
	anAx->lastMonth=99;
	anAx->lastYear=0;

	cpdf_setTimeAxisNumberFormat(anAx, MONTH_NAME, YEAR_2DIGIT, "Helvetica", 12.0);

	return anAx;
}


/* ## MULTI_THREADING_OK ## */
void cpdf_freeAxis(CPDFaxis *anAx)
{
	if(anAx->axisLabel) free(anAx->axisLabel);
	if(anAx->labelFontName) free(anAx->labelFontName);
	if(anAx->labelEncoding) free(anAx->labelEncoding);
	if(anAx->numFormat) free(anAx->numFormat);
	if(anAx->numFontName) free(anAx->numFontName);
        if(anAx->numEncoding) free(anAx->numEncoding); /* alex_V !!!!! */
	if(anAx) free(anAx);
}


/* This function draws the axis -- Because it draws, it needs PDF doc context. */

/* ## MULTI_THREADING_OK ## */
void cpdf_drawAxis(CPDFaxis *anAx)
{
float a, b, c, d, e, f, angle, vcos, vsin;
CPDFdoc *pdf;

    if(anAx->plotDomain == NULL) {
	cpdf_GlobalError(1, "ClibPDF", "Axis is not attached to a plot domain");
	return;
    }
    pdf = anAx->plotDomain->pdfdoc;	/* find the PDF document via axis'es attached domain */
    pdf->anAx2 = anAx;
    angle = PI*(anAx->angle)/180.0;
    vcos = cos(angle);
    vsin = sin(angle);
    a =  vcos;
    b =  vsin;
    c = -vsin;
    d =  vcos;
    /* axis origin is relative to the plotDomain to which it is attached */
    e = (anAx->xloc) + (anAx->plotDomain)->xloc;
    f = (anAx->yloc) + (anAx->plotDomain)->yloc;
    cpdf_comments(pdf, "\n% Axis starts here \n");
    cpdf_gsave(pdf);
    cpdf_rawConcat(pdf, a, b, c, d, e, f);	/* This moves the origin to axis origin, and rotates. */

    /* precompute these scaling factors */
    pdf->xa2points = (anAx->length)/((anAx->valH) - (anAx->valL));
    if(anAx->type == 1) {
	pdf->xaLlog = log10((double)(anAx->valL));
	pdf->xaHlog = log10((double)(anAx->valH));
    }

    /* OK, now draw the main axis line. */
    cpdf_setlinewidth(pdf, anAx->axisLineWidth);
    /* Project out line ends by 1/2 of Major tick linewidth
       This eliminates a notch at the axis ends when ticks are one-sided. */
    cpdf_rawMoveto(pdf, -0.5*anAx->tickWidthMajor, 0.0);
    cpdf_rawLineto(pdf, anAx->length + 0.5*anAx->tickWidthMajor, 0.0);
    cpdf_stroke(pdf);

    if(anAx->type == 0) {
	_do_linearTics(pdf, anAx);
	_do_linearNumbers(pdf, anAx);
    }
    else if(anAx->type == 1) {
	_do_logTics(pdf, anAx);
	_do_logNumbers(pdf, anAx);
    }
    else if(anAx->type == 2) {
	_do_timeTics(pdf, anAx);
	_do_timeNumbers(pdf, anAx);
    }

    _do_axisLabel(pdf, anAx);

    cpdf_grestore(pdf);
}

/* ## MULTI_THREADING_OK ## */
void cpdf_attachAxisToDomain(CPDFaxis *anAx, CPDFplotDomain *domain, float x, float y)
{
	anAx->plotDomain = domain;	/* domain knows the PDF document, so Axis doesn't have to know. */
	anAx->xloc = x;
	anAx->yloc = y;
}

void cpdf_setTicNumEnable(CPDFaxis *anAx, int ticEnableMaj, int ticEnableMin, int numEnable)
{
	anAx->ticEnableMajor = ticEnableMaj;
	anAx->ticEnableMinor = ticEnableMin;
	anAx->numEnable = numEnable;		/* put numbers */
}

void cpdf_setAxisTicNumLabelPosition(CPDFaxis *anAx, int ticPos, int numPos, int horizNum, int horizLab)
{
	anAx->horizNumber = horizNum;
	anAx->horizLabel = horizLab;	/* if true, label will be plotted at horizontal orientation */
	anAx->ticPosition = ticPos;
	anAx->numPosition = numPos;
}


/* ----------------------------------------------------------------------- */
/*  Axis Numbering fuctions */
/* ----------------------------------------------------------------------- */

/* This shouldn't be necessary, but fixes problems with %g formatting with
   some compilers/libraries. Specifically to fix stupid Metrowerks
   ( which prints trailing zeros with %g.
   E.g., fix "0.300000" to correct "0.3" */
char *fix_trailingZeros(char *sstr)
{
char *p, *pe;
int ch;
	pe = sstr + strlen(sstr);	/* on NULL */
	p = strchr(sstr, '.');
	if(p != NULL) {
	    /* we have a decimal point */
	    pe--;		/* on the last char */
	    while(((ch = *pe) == '0') && pe > p)
		pe--;
	    pe++;
	    *pe = '\0';
	}
	return(sstr);
}

void _do_oneNumber(CPDFdoc *pdf, CPDFaxis *anAx, float v, float ticlen)
{
float startX, startY, ticstartN, ticstartP, vt, ang;
float slen, hslen;
char str[32];
	sprintf(str, anAx->numFormat, v);
	if(anAx->type == 1)
	    fix_trailingZeros(str);	/* If LOG axis. (Only MacOS Metrowerks CW needs it, IDE 1.7.4) */
	slen = cpdf_stringWidth(pdf, (unsigned char *)str);	/* number string width for centring */
	hslen = slen*0.5;
	ang = PI*anAx->angle/180.0;
	vt = vAxis2Points(pdf, v);		/* convert domain value for number to points */
	if(anAx->ticPosition == 0) { ticstartN = -ticlen; ticstartP = 0.0;}
	else if(anAx->ticPosition == 1) { ticstartN = -0.5*ticlen; ticstartP = 0.5*ticlen;}
	else { ticstartN = 0.0; ticstartP = ticlen;}
	if(anAx->horizNumber) {
	    /* number stays always horizontal */
	    if(anAx->numPosition == 0) {
		startY = ticstartN - (anAx->ticNumGap) - (anAx->numFontSize)*0.6;
		pdf->numEdgeY = startY - 2.5;	/* store edge of numbers for positioning label */
	    }
	    else {
		startY = ticstartP + (anAx->ticNumGap) + slen;
		pdf->numEdgeY = startY + 4.0;
		/* FIXME: note a potential problem -- this function is called multiple times for
			  for a given axis, pdf->numEdgeY will store the last value.  If this is not
			  the longest string, there is a problem.
		*/
	    }
	    startX = vt - ((anAx->numFontSize)*0.3*sin(ang) + hslen*cos(ang));
	    cpdf_rawText(pdf, startX, startY, -(anAx->angle), str);
	}
	else {
	    /* number is parallel to axis and rotate with axis */
	    if(anAx->numPosition == 0) {
		startY = ticstartN - (anAx->ticNumGap) - (anAx->numFontSize)*0.6;
		pdf->numEdgeY = startY;	/* store edge of numbers for positioning label */
	    }
	    else {
		startY = ticstartP + (anAx->ticNumGap);
		pdf->numEdgeY = startY + (anAx->numFontSize);
	    }
	    startX = vt - hslen;
	    cpdf_rawText(pdf, startX, startY, 0.0, str);
	}
}

char *_yearFormat(int year, int flag, char *stbuf)
{
	if(flag)
	    sprintf(stbuf, "%02d", year%100);
	else
	    sprintf(stbuf, "%d", year);
	return(stbuf);
}

int cpdf_setMonthNames(CPDFdoc *pdf, char *mnArray[])
{
int i;
    _cpdf_freeMonthNames(pdf);
    for(i=0; i<12; i++) {
	pdf->monthName[i] = (char *)malloc((size_t)(strlen(mnArray[i]) + 1));
	_cpdf_malloc_check((void *)pdf->monthName[i]);
	strcpy(pdf->monthName[i], mnArray[i]);
    }
    return(0);
}

int _cpdf_freeMonthNames(CPDFdoc *pdf)
{
int i;
    for(i=0; i<12; i++) {
     	if(pdf->monthName[i])
	    free(pdf->monthName[i]);
    }
    return(0);
}

/* now in CPDFaxis struct and initialized by createTimeAxis ... */
/* static int lastMin=99, lastHour=99, lastDay=99, lastMonth=99, lastYear=0; */

void _do_oneTimeNumber(CPDFdoc *pdf, CPDFaxis *anAx, float v, struct tm *vtm, int majorBumpVar, float ticlen)
{
float startX, startY, ticstartN, ticstartP, vt, ang;
float slen, hslen;
char str[32], stryear[8];
	/* _printfTime(vtm); */		/* debugging */
	switch(majorBumpVar) {
	    case SECOND:
			if(anAx->lastMin != vtm->tm_min)
			    sprintf(str, "%d:%d", vtm->tm_min, vtm->tm_sec);
			else
			    sprintf(str, "%d", vtm->tm_sec);
			break;
	    case MINUTE:
			if(anAx->lastHour != vtm->tm_hour)
			    sprintf(str, "%d:%d", vtm->tm_hour, vtm->tm_min);
			else
			    sprintf(str, "%d", vtm->tm_min);
			break;
	    case HOUR:
			if(anAx->lastDay != vtm->tm_mday)
			    sprintf(str, "%d [%d]", vtm->tm_hour, vtm->tm_mday);
			else
			    sprintf(str, "%d", vtm->tm_hour);
			break;
	    case DAY:
			if(anAx->lastMonth != vtm->tm_mon) {
			    if(anAx->useMonthName)
				sprintf(str, "%s %d", pdf->monthName[vtm->tm_mon], vtm->tm_mday);
			    else
				sprintf(str, "%d/%d", vtm->tm_mon, vtm->tm_mday);
			}
			else if(vtm->tm_mday < 30)
			    sprintf(str, "%d", vtm->tm_mday);
			break;
	    case MONTH:
			if(anAx->lastYear != vtm->tm_year) {
			    _yearFormat(vtm->tm_year + 1900, anAx->use2DigitYear, stryear);
			    if(anAx->useMonthName)
				sprintf(str, "%s %s", pdf->monthName[vtm->tm_mon], stryear);
			    else
				sprintf(str, "%d/%s", vtm->tm_mon+1, stryear );
			}
			else {
			    if(anAx->useMonthName)
				sprintf(str, "%s", pdf->monthName[vtm->tm_mon]);
			    else
				sprintf(str, "%d", vtm->tm_mon+1);
			}
			break;
	    case YEAR:
			_yearFormat(vtm->tm_year + 1900, anAx->use2DigitYear, stryear);
			sprintf(str, "%s", stryear);
			break;
	}
	/* sprintf(str, anAx->numFormat, v); */
	if(anAx->type == 1)
	    fix_trailingZeros(str);	/* If LOG axis. (Only MacOS Metrowerks CW needs it, IDE 1.7.4) */
	slen = cpdf_stringWidth(pdf, (unsigned char *)str);	/* number string width for centring */
	hslen = slen*0.5;
	ang = PI*anAx->angle/180.0;
	vt = vAxis2Points(pdf, v);		/* convert domain value for number to points */
	if(anAx->ticPosition == 0) { ticstartN = -ticlen; ticstartP = 0.0;}
	else if(anAx->ticPosition == 1) { ticstartN = -0.5*ticlen; ticstartP = 0.5*ticlen;}
	else { ticstartN = 0.0; ticstartP = ticlen;}
	if(anAx->horizNumber) {
	    /* number stays always horizontal */
	    if(anAx->numPosition == 0) {
		startY = ticstartN - (anAx->ticNumGap) - (anAx->numFontSize)*0.6;
		pdf->numEdgeY = startY - 2.5;	/* store edge of numbers for positioning label */
	    }
	    else {
		startY = ticstartP + (anAx->ticNumGap) + slen;
		pdf->numEdgeY = startY + 4.0;
		/* FIXME: note a potential problem -- this function is called multiple times for
			  for a given axis, pdf->numEdgeY will store the last value.  If this is not
			  the longest string, there is a problem.
		*/
	    }
	    /* startX = vt - ((anAx->numFontSize)*0.3*sin(ang) + hslen*cos(ang)); */
	    /* do not center date/time values on tick */
	    startX = vt - ((anAx->numFontSize)*0.3*sin(ang) + (anAx->numFontSize)*0.25);
	    cpdf_rawText(pdf, startX, startY, -(anAx->angle), str);
	}
	else {
	    /* number is parallel to axis and rotate with axis */
	    if(anAx->numPosition == 0) {
		startY = ticstartN - (anAx->ticNumGap) - (anAx->numFontSize)*0.6;
		pdf->numEdgeY = startY;	/* store edge of numbers for positioning label */
	    }
	    else {
		startY = ticstartP + (anAx->ticNumGap);
		pdf->numEdgeY = startY + (anAx->numFontSize);
	    }
	    /* startX = vt - hslen; */
	    startX = vt - (anAx->numFontSize)*0.25;  /* do not center date/time values on tick */
	    cpdf_rawText(pdf, startX, startY, 0.0, str);
	}
	/* Save time fields for showing next field up when it rolls up */
	anAx->lastMin = vtm->tm_min;
	anAx->lastHour = vtm->tm_hour;
	anAx->lastDay = vtm->tm_mday;
	anAx->lastMonth = vtm->tm_mon;
	anAx->lastYear = vtm->tm_year;
}



void _do_linearNumbers(CPDFdoc *pdf, CPDFaxis *anAx)
{
float vL, vH, v;
	cpdf_beginText(pdf, 0);
	cpdf_setFont(pdf, anAx->numFontName, anAx->numEncoding, anAx->numFontSize);
	vL = anAx->valL;
	vH = anAx->valH;
	/* do numbers at major ticks */
	for(v = anAx->valFirstTicLinMajor; v <= vH*1.0001; v += anAx->ticIntervalLinMajor) {
	    _do_oneNumber(pdf, anAx, v, anAx->ticLenMajor);
	}
	cpdf_endText(pdf);
}

void _do_logNumbers(CPDFdoc *pdf, CPDFaxis *anAx)
{
float vL, vH, v;
int m, ep;
int mL, expL, mH, expH;	/* mantissa MSD and exponents */
	cpdf_beginText(pdf, 0);
	/* Encoding shouldn't matter for numbers on the axis. */
	cpdf_setFont(pdf, anAx->numFontName, anAx->numEncoding, anAx->numFontSize);
	vL = anAx->valL;
	vH = anAx->valH;
	mL = (int)getMantissaExp(vL*1.0001, &expL);
	mH = (int)getMantissaExp(vH*1.0001, &expH);
	/* fprintf(stderr, "%dE%d -- %dE%d\n", mL, expL, mH, expH); */
	/* do numbers at major ticks */
	for(m=mL, ep=expL; (v = (float)m * pow(10.0, (double)ep)) <= vH*1.0001; ) {
	    if( _bittest(anAx->numSelectorLog, m) )	/* Test if we want this number */
	        _do_oneNumber(pdf, anAx, v, anAx->ticLenMajor);
	    m++;
	    if(m >= 10) {	/* cycle m through 1 .. 9 */
		m = 1;
		ep++;
	    }
	}
	cpdf_endText(pdf);
}


void _do_timeNumbers(CPDFdoc *pdf, CPDFaxis *anAx)
{
float fndays, v;
int minorBump = 1, majorBump = 1;
int minorBumpVar = MINUTE, majorBumpVar = HOUR;
struct tm vtm;
	anAx->lastMin=99; anAx->lastHour=99; anAx->lastDay=99; anAx->lastMonth=99; anAx->lastYear=0;
	cpdf_beginText(pdf, 0);
	/* Encoding shouldn't matter for numbers on the axis. */
	cpdf_setFont(pdf, anAx->numFontName, anAx->numEncoding, anAx->numFontSize);

	/* float number of days spanned by the axis */
	fndays = tm_to_NumDays(&anAx->vTL, &anAx->vTH);
	_setDefaultTimeBumpVar(fndays, &minorBumpVar, &majorBumpVar, &minorBump, &majorBump);
	/* do numbers at major ticks */
	memcpy(&vtm, &anAx->vTL, sizeof(struct tm));	/* initialize vtm to start of axis */
	for(v = tm_to_NumDays(&anAx->vTL, &vtm) ; v <= fndays*1.0001 ;
		v = _bump_tm_Time(&anAx->vTL, &vtm, majorBumpVar, majorBump)) {
	    _do_oneTimeNumber(pdf, anAx, v, &vtm, majorBumpVar, anAx->ticLenMajor);
	}
	cpdf_endText(pdf);
}

/*                              1   2   3   4   5   6   7   8   9  10  11  12 */
static int daysInMon[][12] = {{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
			      {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}};

float _bump_tm_Time(struct tm *rT, struct tm *vT, int bumpVar, int bump)
{
float fndays;
int leapY = 0;
	cpdf_mktime(vT);		/* correct all fields of struct tm */
        leapY = isLeapYear(vT->tm_year + 1900);
	switch(bumpVar) {
	    case SECOND:
			vT->tm_sec += bump;
			break;
	    case MINUTE:
			vT->tm_min += bump;
			break;
	    case HOUR:
			vT->tm_hour += bump;
			break;
	    case DAY:
			vT->tm_mday += bump;
			if(vT->tm_mday > daysInMon[leapY][vT->tm_mon]) {
			    vT->tm_mday = 1;
			    vT->tm_mon++;
			}
			else {
			    /* do: [1, 5, 10, 15, 20, 25, 30] or [1, 10, 20, 30] */
			    if(( (bump == 5 && vT->tm_mday % 5 == 1) ||
				 (bump == 10 && vT->tm_mday % 10 == 1) ) && vT->tm_mday > 5 )
				vT->tm_mday--;
			}
			break;
	    case MONTH:
			vT->tm_mon += bump;
			break;
	    case YEAR:
			vT->tm_year += bump;
			break;
	    default:
			break;
	}
	fndays = tm_to_NumDays(rT, vT);		/* This calls cpdf_mktime() for both internally */
	return(fndays);
}


/*  You may edit/add entries in this table to customize Time domain mesh and
	axis tick spacings.
	Tcut 	   : max # of days (float value) that the entry applies.
	MinorVar   : variable that defines minor tick intervals
	MajorVar   : variable that defines major tick intervals
	MinorBump  : Minor tick interval in the unit of MinorVar (interger)
	MajorBump  : Major tick interval in the unit of MajorVar (interger)
*/

static CPDFtimeLU timeLU[] = {
    /*	  Tcut	 MinorVar, 	MajorVar	MinorBump	MajorBump */
	{ 0.042,   MINUTE,	MINUTE,		 1,		10 },	/* 1 hour or less */
	{ 0.084,   MINUTE,	HOUR,		10,		1 },	/* 2 hours or less */
	{  0.51,   MINUTE,	HOUR,		30,		1 },	/* 12 hours or less */
	{   1.1,   HOUR,	HOUR,		1,		6 },	/* 1 day or less */
	{   8.0,   HOUR,	DAY,		6,		1 },	/* 8 days or less */
	{  12.5,   HOUR,	DAY,		12,		1 },	/* 12 days or less */
	{  21.5,   HOUR,	DAY,		12,		5 },	/* 21 days or less */
	{  42.0,   DAY,		DAY,		1,		10 },	/* 42 days or less */
	{ 190.0,   DAY,		MONTH,		10,		1 },	/* 190 days or less */
	{ 367.0,   MONTH,	MONTH,		1,		3 },	/* 1 year of less */
	{ 735.0,   MONTH,	MONTH,		1,		6 },	/* 2 years or less */
	{ 2193.0,  MONTH,	YEAR,		6,		1 },	/* 6 years or less */
	{ 7306.0,  YEAR,	YEAR,		1,		5 },	/* 20 years or less */
	{ 21916.0, YEAR,	YEAR,		1,		10 },	/* 60 years or less */
	{ 54788.0, YEAR,	YEAR,		5,		10 },	/* 150 years or less */
	{ 99999.0, YEAR,	YEAR,		10,		100 }	/* long time */
};

#define	NtimeLU		(sizeof(timeLU)/sizeof(CPDFtimeLU))

void _setDefaultTimeBumpVar(float fndays, int *minorBumpVar, int *majorBumpVar, int *minorBump, int *majorBump)
{
int i, idx =-1;

	*minorBumpVar = MINUTE;
	*majorBumpVar = HOUR;
	*minorBump = 10;
	*majorBump = 1;

	for(i=0; i < NtimeLU; i++) {
	    if(  fndays <= timeLU[i].Tcut ) {
		*minorBumpVar = timeLU[i].MinorVar;
		*majorBumpVar = timeLU[i].MajorVar;
		*minorBump = timeLU[i].MinorBump;
		*majorBump = timeLU[i].MajorBump;
		idx = i;
		break;
	    }
	}

	if(idx < 0) {
	    *minorBumpVar = YEAR;
	    *majorBumpVar = YEAR;
	    *minorBump = 10;
	    *majorBump = 100;
	}
}


void cpdf_suggestTimeDomainParams(struct tm *xTL, struct tm *xTH, struct tm *recTL, struct tm *recTH)
{
float fndays;
int minorBump = 1, majorBump = 1;
int minorBumpVar = MINUTE, majorBumpVar = HOUR;
int minorBumpSave = 1, majorBumpSave = 1;
int minorBumpVarSave = MINUTE, majorBumpVarSave = HOUR;
struct tm tempTL, tempTH;

    memcpy(recTL, xTL, sizeof(struct tm));	/* copy given limits to recTL */
    memcpy(recTH, xTH, sizeof(struct tm));	/* copy given limits to recTH */
    fndays = tm_to_NumDays(xTL, xTH);
    _setDefaultTimeBumpVar(fndays, &minorBumpVar, &majorBumpVar, &minorBump, &majorBump);
    minorBumpSave = minorBump;
    majorBumpSave = majorBump;
    minorBumpVarSave = minorBumpVar;
    majorBumpVarSave = majorBumpVar;

    switch(majorBumpVar) {
	case SECOND:
		    break;
	case MINUTE:
		    recTL->tm_min = (xTL->tm_min/majorBump) * majorBump;
		    if(recTL->tm_sec) {
			recTL->tm_sec = 0;
		    }
		    if(recTH->tm_sec) {
			recTH->tm_sec = 0;
			recTH->tm_min = (xTH->tm_min/majorBump +1) * majorBump;
		    }
		    break;
	case HOUR:
		    recTL->tm_hour = (xTL->tm_hour/majorBump) * majorBump;
		    if(recTL->tm_min) {
			recTL->tm_sec = 0;
			recTL->tm_min = 0;
		    }
		    if(recTH->tm_min) {
			recTH->tm_sec = 0;
			recTH->tm_min = 0;
			recTH->tm_hour = (xTH->tm_hour/majorBump +1) * majorBump;
		    }
		    break;
	case DAY:
		    if((xTL->tm_mday - 1)/majorBump)
		        recTL->tm_mday = ((xTL->tm_mday - 1)/majorBump) * majorBump;
		    else
		        recTL->tm_mday = ((xTL->tm_mday - 1)/majorBump) * majorBump +1;
		    if(recTL->tm_hour) {
			recTL->tm_sec = 0;
			recTL->tm_min = 0;
			recTL->tm_hour = 0;
		    }
		    if(recTH->tm_hour) {
			recTH->tm_sec = 0;
			recTH->tm_min = 0;
			recTH->tm_hour = 0;
			recTH->tm_mday = ((xTH->tm_mday - 1)/majorBump +1) * majorBump +1;
		    }
		    break;
	case MONTH:
		    recTL->tm_mon = (xTL->tm_mon/majorBump) * majorBump;
		    if(recTL->tm_mday != 1) {
			recTL->tm_sec = 0;
			recTL->tm_min = 0;
			recTL->tm_hour = 0;
			recTL->tm_mday = 1;	/* day of month starts at 1 */
		    }
		    if(recTH->tm_mday != 1) {
			recTH->tm_sec = 0;
			recTH->tm_min = 0;
			recTH->tm_hour = 0;
			recTH->tm_mday = 1;	/* day of month starts at 1 */
			recTH->tm_mon = (xTH->tm_mon/majorBump +1) * majorBump;
		    }
		    break;
	case YEAR:
		    recTL->tm_year = (xTL->tm_year/majorBump) * majorBump;
		    if(recTL->tm_mon) {
			recTL->tm_sec = 0;
			recTL->tm_min = 0;
			recTL->tm_hour = 0;
			recTL->tm_mday = 1;	/* day of month starts at 1 */
			recTL->tm_mon  = 0;
		    }
		    if(recTH->tm_mon) {
			recTH->tm_sec = 0;
			recTH->tm_min = 0;
			recTH->tm_hour = 0;
			recTH->tm_mday = 1;	/* day of month starts at 1 */
			recTH->tm_mon = 0;
			recTH->tm_year = (xTH->tm_year/majorBump +1) * majorBump;
		    }
		    break;
	default:
		    break;
    }  /* end of switch() */

    cpdf_mktime(recTL); cpdf_mktime(recTH);	/* This sets all fields of struct tm correctly */
    /* check if new values don't change the major and minor variables and bumps */
    fndays = tm_to_NumDays(recTL, recTH);
    _setDefaultTimeBumpVar(fndays, &minorBumpVar, &majorBumpVar, &minorBump, &majorBump);
    if(minorBumpVar == minorBumpVarSave && majorBumpVar == majorBumpVarSave &&
	minorBump == minorBumpSave && majorBump == majorBumpSave)
	return;		/* all are the same, so OK to return */
    /* Ooops, something has changed by that little extension, so redo it */
    memcpy(&tempTL, recTL, sizeof(struct tm));	/* copy given limits to recTL */
    memcpy(&tempTH, recTH, sizeof(struct tm));	/* copy given limits to recTH */
    cpdf_suggestTimeDomainParams(&tempTL, &tempTH, recTL, recTH);	/* recursive call */
}


/* return non-zero if aNumber (in binary representation) has bitpos ON */
/* bitpos must be 0 .. MSB of int */
int _bittest(int aNumber, int bitpos)
{
int tmask = 1;
	tmask <<= bitpos;
	if(aNumber & tmask) return(1);
	else		    return(0);
}

/* ## MULTI_THREADING_OK ## */
void cpdf_setAxisNumberFormat(CPDFaxis *anAx, const char *format, const char *fontName, float fontSize)
{
    cpdf_setAxisNumberFormat2(anAx, format, fontName, fontSize, NULL);
}

void cpdf_setAxisNumberFormat2(CPDFaxis *anAx, const char *format,
		const char *fontName, float fontSize, const char *encoding)
{
char *encodeName;
	if(encoding) encodeName = (char *)encoding;
	else encodeName = "MacRomanEncoding";		/* for backward compatibility */
	if(anAx->numFormat) free(anAx->numFormat);
	if(anAx->numFontName) free(anAx->numFontName);
	if(anAx->numEncoding) free(anAx->numEncoding);
	anAx->numFormat = (char *)malloc((size_t)(strlen(format) + 1));
	_cpdf_malloc_check((void *)anAx->numFormat);
	anAx->numFontName = (char *)malloc((size_t)(strlen(fontName) + 1));
	_cpdf_malloc_check((void *)anAx->numFontName);
	anAx->numEncoding = (char *)malloc((size_t)(strlen(encodeName) + 1));
	_cpdf_malloc_check((void *)anAx->numEncoding);

	strcpy(anAx->numFormat, format);
	strcpy(anAx->numFontName, fontName);
	strcpy(anAx->numEncoding, encodeName);
	anAx->numFontSize = fontSize;
}

/* ## MULTI_THREADING_OK ## */
void cpdf_setTimeAxisNumberFormat(CPDFaxis *anAx, int useMonName, int use2DigYear,
		const char *fontName, float fontSize)
{
    cpdf_setTimeAxisNumberFormat2(anAx, useMonName, use2DigYear, fontName, fontSize, NULL);
}

void cpdf_setTimeAxisNumberFormat2(CPDFaxis *anAx, int useMonName, int use2DigYear,
		const char *fontName, float fontSize, const char *encoding)
{
char *encodeName;
	if(encoding) encodeName = (char *)encoding;
	else encodeName = "MacRomanEncoding";		/* for backward compatibility */
	if(anAx->numFontName) free(anAx->numFontName);
	if(anAx->numEncoding) free(anAx->numEncoding);
	anAx->numFontName = (char *)malloc((size_t)(strlen(fontName) + 1));
	_cpdf_malloc_check((void *)anAx->numFontName);
	anAx->numEncoding = (char *)malloc((size_t)(strlen(encodeName) + 1));
	_cpdf_malloc_check((void *)anAx->numEncoding);

	strcpy(anAx->numFontName, fontName);
	strcpy(anAx->numEncoding, encodeName);
	anAx->numFontSize = fontSize;
	anAx->use2DigitYear = use2DigYear;
	anAx->useMonthName = useMonName;
}


/* ----------------------------------------------------------------------- */
/*  Axis Tic fuctions */
/* ----------------------------------------------------------------------- */

/* ## MULTI_THREADING_OK ## */
void cpdf_setAxisLineParams(CPDFaxis *anAx, float axLineWidth, float ticLenMaj, float ticLenMin,
				float tickWidMaj, float tickWidMin)
{
	anAx->ticLenMajor = ticLenMaj;
	anAx->ticLenMinor = ticLenMin;
	anAx->axisLineWidth = axLineWidth;
	anAx->tickWidthMajor = tickWidMaj;
	anAx->tickWidthMinor = tickWidMin;
}


void _do_oneTick(CPDFdoc *pdf, CPDFaxis *anAx, float vt, float ticlen)
{
float tstart, tend;
	if(anAx->ticPosition == 0) {
	     tstart = -ticlen;
	     tend = 0.0;
	     /* tend = 0.4 * anAx->axisLineWidth; */	/* Adjust 1/2 linewidth to avoid a notch at axis ends */
	}
	else if(anAx->ticPosition == 1) {
	    tstart = -0.5*ticlen;
	    tend = 0.5*ticlen;
	}
	else {
	    /* tstart = -0.4 * anAx->axisLineWidth; */ /* Adjust 1/2 linewidth to avoid a notch at axis ends */
	    tstart = 0.0;
	    tend = ticlen;
	}
	cpdf_rawMoveto(pdf, vt, tstart);
	cpdf_rawLineto(pdf, vt, tend);
}


void _do_linearTics(CPDFdoc *pdf, CPDFaxis *anAx)
{
float vL, vH, v, vt;
    vL = anAx->valL;
    vH = anAx->valH;

    /* do minor ticks */
    if(anAx->ticEnableMinor) {
        cpdf_setlinewidth(pdf, anAx->tickWidthMinor);
	for(v = anAx->valFirstTicLinMinor; v <= vH*1.0001; v += anAx->ticIntervalLinMinor) {
	    vt = vAxis2Points(pdf, v);
	    _do_oneTick(pdf, anAx, vt, anAx->ticLenMinor);
	}
	cpdf_stroke(pdf);
    }

    /* do major ticks */
    if(anAx->ticEnableMajor) {
        cpdf_setlinewidth(pdf, anAx->tickWidthMajor);
	for(v = anAx->valFirstTicLinMajor; v <= vH*1.0001; v += anAx->ticIntervalLinMajor) {
	    vt = vAxis2Points(pdf, v);
	    _do_oneTick(pdf, anAx, vt, anAx->ticLenMajor);
	}
	cpdf_stroke(pdf);
    }
}


void _do_logTics(CPDFdoc *pdf, CPDFaxis *anAx)
{
float vL, vH, v, vt;
int m, ep;
int mL, expL, mH, expH;	/* mantissa MSD and exponents */
	vL = anAx->valL;
	vH = anAx->valH;
	mL = (int)getMantissaExp(vL*1.0001, &expL);
	mH = (int)getMantissaExp(vH*1.0001, &expH);
	/* fprintf(stderr, "%dE%d -- %dE%d\n", mL, expL, mH, expH); */

	/* do ticks at 1 sig digits */
	for(m=mL, ep=expL; (v = (float)m * pow(10.0, (double)ep)) <= vH*1.0001; ) {
	    /* fprintf(stderr, "tick %g\n", v); */
	    vt = vAxis2Points(pdf, v);
	    if( _bittest(anAx->ticSelectorLog, m) ) {	/* Test if we want this tick */
		if(m==1) {
        	    cpdf_setlinewidth(pdf, anAx->tickWidthMajor);
		    _do_oneTick(pdf, anAx, vt, anAx->ticLenMajor);
		    cpdf_stroke(pdf);
		}
		else {
        	    cpdf_setlinewidth(pdf, anAx->tickWidthMinor);
		    _do_oneTick(pdf, anAx, vt, anAx->ticLenMinor);
		    cpdf_stroke(pdf);
		}
	    }
	    m++;
	    if(m >= 10) {
		m = 1;
		ep++;
	    }
	}
	/* cpdf_stroke(pdf); */
}


void _do_timeTics(CPDFdoc *pdf, CPDFaxis *anAx)
{
float v, vt;
float fndays;
int minorBump = 1, majorBump = 1;
int minorBumpVar = MINUTE, majorBumpVar = HOUR;
struct tm vtm;

    fndays = tm_to_NumDays(&anAx->vTL, &anAx->vTH);
    _setDefaultTimeBumpVar(fndays, &minorBumpVar, &majorBumpVar, &minorBump, &majorBump);

    /* do minor ticks */
    if(anAx->ticEnableMinor) {
	cpdf_setlinewidth(pdf, anAx->tickWidthMinor);
	memcpy(&vtm, &anAx->vTL, sizeof(struct tm));
	for(v = tm_to_NumDays(&anAx->vTL, &vtm) ; v <= fndays*1.0001 ;
			v = _bump_tm_Time(&anAx->vTL, &vtm, minorBumpVar, minorBump)) {
	    /* _printfTime(&vtm); */
	    vt = vAxis2Points(pdf, v);
	    _do_oneTick(pdf, anAx, vt, anAx->ticLenMinor);
	}
	cpdf_stroke(pdf);
    }

    /* do major ticks */
    if(anAx->ticEnableMajor) {
	cpdf_setlinewidth(pdf, anAx->tickWidthMajor);
	memcpy(&vtm, &anAx->vTL, sizeof(struct tm));
	for(v = tm_to_NumDays(&anAx->vTL, &vtm) ; v <= fndays*1.0001 ;
			v = _bump_tm_Time(&anAx->vTL, &vtm, majorBumpVar, majorBump)) {
	    /* _printfTime(&vtm); */
	    vt = vAxis2Points(pdf, v);
	    _do_oneTick(pdf, anAx, vt, anAx->ticLenMajor);
	}
	cpdf_stroke(pdf);
    }
}

void cpdf_setLinearAxisParams(CPDFaxis *anAx, float tic1ValMajor, float intervalMajor,
			float tic1ValMinor, float intervalMinor)
{
	anAx->valFirstTicLinMajor = tic1ValMajor;
	anAx->valFirstTicLinMinor = tic1ValMinor;
	anAx->ticIntervalLinMajor = intervalMajor;
	anAx->ticIntervalLinMinor = intervalMinor;
}


void cpdf_setLogAxisTickSelector(CPDFaxis *anAx, int ticselect)
{
	anAx->ticSelectorLog = ticselect;
}

void cpdf_setLogAxisNumberSelector(CPDFaxis *anAx, int numselect)
{
	anAx->numSelectorLog = numselect;
}


/* ----------------------------------------------------------------------- */
/*  Axis label fuctions */
/* ----------------------------------------------------------------------- */

void cpdf_setAxisLabel(CPDFaxis *anAx, const char *labelstring, const char *fontName,
	const char *encoding, float fontSize)
{
	if(anAx->axisLabel) free(anAx->axisLabel);
	anAx->axisLabel = (char *)malloc((size_t)(strlen(labelstring) + 1));
	_cpdf_malloc_check((void *)anAx->axisLabel);
	strcpy(anAx->axisLabel, labelstring);

	if(anAx->labelFontName) free(anAx->labelFontName);
	anAx->labelFontName = (char *)malloc((size_t)(strlen(fontName) + 1));
	_cpdf_malloc_check((void *)anAx->labelFontName);
	strcpy(anAx->labelFontName, fontName);

	if(anAx->labelEncoding) free(anAx->labelEncoding);
	anAx->labelEncoding = (char *)malloc((size_t)(strlen(encoding) + 1));
	_cpdf_malloc_check((void *)anAx->labelEncoding);
	strcpy(anAx->labelEncoding, encoding);

	anAx->labelFontSize = fontSize;
}


void _do_axisLabel(CPDFdoc *pdf, CPDFaxis *anAx)
{
float startX, startY, v, vt, ang;
float slen, hslen;
	if(anAx->axisLabel == NULL) return;		/* no axis label */

	cpdf_beginText(pdf, 0);
	cpdf_setFont(pdf, anAx->labelFontName, anAx->labelEncoding, anAx->labelFontSize);
	slen = cpdf_stringWidth(pdf, (unsigned char *)(anAx->axisLabel));	/* number string width for centring */
	hslen = slen*0.5;
	ang = PI*anAx->angle/180.0;
	/* find the value for half point of axis to center label */
	if(anAx->type == 1)
	    v = sqrt(anAx->valH / anAx->valL) * anAx->valL;	/* log axis mid value */
	else
	    v = 0.5* (anAx->valH - anAx->valL) + anAx->valL;	/* linear axis mid value */
	vt = vAxis2Points(pdf, v);		/* X pos of the center of axis in points */
	if(anAx->horizLabel) {
	    /* number stays always horizontal */
	    if(anAx->numPosition == 0)		/* numPosition also controls label position */
		startY = pdf->numEdgeY - (anAx->numLabelGap) - (anAx->labelFontSize)*0.6;
	    else
		startY = pdf->numEdgeY + (anAx->numLabelGap) + slen;
	    startX = vt - ((anAx->labelFontSize)*0.3*sin(ang) + hslen*cos(ang));
	    cpdf_rawText(pdf, startX, startY, -(anAx->angle), anAx->axisLabel);
	}
	else {
	    /* number is parallel to axis and rotate with axis */
	    if(anAx->numPosition == 0)		/* numPosition also controls label position */
		startY = pdf->numEdgeY - (anAx->numLabelGap) - (anAx->labelFontSize)*0.6;
	    else
		startY = pdf->numEdgeY + (anAx->numLabelGap);
	    startX = vt - hslen;
	    cpdf_rawText(pdf, startX, startY, 0.0, anAx->axisLabel);
	}
	cpdf_endText(pdf);
}


/* ----------------------------------------------------------------------- */
/*  Utility fuctions */
/* ----------------------------------------------------------------------- */


float vAxis2Points(CPDFdoc *pdf, float x)
{
float xrval = 0.0;
double xvlog = 0.0, fraction=0.0;

    /* xrval = pdf->anAx2->xloc; */	/* min point value for X */
    if(pdf->anAx2->type == 0) {
	/* linear */
	xrval += (x - (pdf->anAx2->valL)) * pdf->xa2points;
    }
    else if(pdf->anAx2->type == 1) {
	/* logarithmic */
	if(x > 0.0) {
	    xvlog = log10((double)x);
	    fraction = (xvlog - pdf->xaLlog)/(pdf->xaHlog - pdf->xaLlog); /* xHlog, xLlog precomputed */
	    xrval += (pdf->anAx2->length) * fraction;
	}
    }
    else if(pdf->anAx2->type == 2) {
	/* time */
	xrval += (x - (pdf->anAx2->valL)) * pdf->xa2points;
    }
    return xrval;
}



