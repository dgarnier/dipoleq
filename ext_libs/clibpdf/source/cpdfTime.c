/* cpdfTime.c -- Time/Date string management.
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

cc -Wall -DCPDF_DEBUG -o cpdfTime cpdfTime.c


1999-08-17 [IO] 
	## MULTI_THREADING_OK ##

*/

#include <stdio.h>
#include <time.h>
#include "cpdflib.h"		/* This must be included before all other local include files */
#include "cpdf_private.h"

/* To make it thread-safe, we let the caller pass char *timebuf now */

char *timestring(int tformat, char *TimeBuf)
{
time_t t;
struct tm tmbuf;	/* buffer for thread-safe localtime_r() */
struct tm *tm;

    time(&t);
    tm = cpdf_localtime(&t, &tmbuf);

    switch(tformat){
	case 0:		/* 19981231235959 */
	    sprintf(TimeBuf, "%04d%02d%02d%02d%02d%02d",
		(tm->tm_year) +1900,
		(tm->tm_mon) +1,
		tm->tm_mday,
		tm->tm_hour,
		tm->tm_min,
		tm->tm_sec);
		break;

	case 1:		/* 1998-12-31 23:59:59 */
	default:
	    sprintf(TimeBuf, "%04d-%02d-%02d %02d:%02d:%02d",
		(tm->tm_year) +1900,
		(tm->tm_mon) +1,
		tm->tm_mday,
		tm->tm_hour,
		tm->tm_min,
		tm->tm_sec);
		break;
    } /* end switch() */
  return(TimeBuf);
}

float tm_to_NumDays(struct tm *fromDate, struct tm *toDate)
{
float fndays = 0.0, fndaysFrom=0.0, fndaysTo =0.0;
int i, fromYear, toYear;
long ndays = 0;
	cpdf_mktime(fromDate);	/* This sets all fields of struct tm correctly */
	cpdf_mktime(toDate);
	toYear = 1900 + toDate->tm_year;
	fromYear = 1900 + fromDate->tm_year;
	fndaysFrom = (float)(fromDate->tm_yday);	/* # of days since Jan 1 of fromDate */
	fndaysFrom += (float)(fromDate->tm_hour)/24.0;		/* add hour contribution */
	fndaysFrom += (float)(fromDate->tm_min)/1440.0;		/* add minute contribution */
	fndaysFrom += (float)(fromDate->tm_sec)/86400.0;	/* add second contribution */

	fndaysTo = (float)(toDate->tm_yday);		/* # of days since Jan 1 of ToDate */
	fndaysTo += (float)(toDate->tm_hour)/24.0;		/* add hour contribution */
	fndaysTo += (float)(toDate->tm_min)/1440.0;		/* add minute contribution */
	fndaysTo += (float)(toDate->tm_sec)/86400.0;		/* add second contribution */

	for(i=fromYear; i < toYear; i++) {
	     if(isLeapYear(i)) ndays += 366;
	     else ndays += 365;
	}
	fndays = (float)ndays + fndaysTo;
	fndays -= fndaysFrom;
	return(fndays);
}

int isLeapYear(int year)
{
    if( ((year % 4 == 0) && !(year % 100 == 0)) || (year % 400 == 0) )
	 return(1);
    else return(0);
}

void _printfTime(struct tm *vtm)
{
    printf("%04d-%02d-%02d %02d:%02d:%02d\n",
	(vtm->tm_year) +1900,
	(vtm->tm_mon) +1,
	vtm->tm_mday,
	vtm->tm_hour,
	vtm->tm_min,
	vtm->tm_sec);
}

#ifdef CPDF_DEBUG
#include <string.h>

void main(void)
{
time_t t;
struct tm tmbuf;	/* buffer for thread-safe localtime_r() */
struct tm tm;
struct tm tm2;
char timebuf[64];

    time(&t);
    memcpy(&tm, cpdf_localtime(&t, &tmbuf), sizeof(struct tm));
    memcpy(&tm2, &tm, sizeof(struct tm));
    tm2.tm_year += 2;		/* +2 years */
    tm2.tm_hour += 6;		/* +6 hours */
    tm2.tm_mon += 12;		/* +1 year */
    printf("Before - tm2.tm_mon= %d\n", tm2.tm_mon);
    printf("Number of days (3 years and 6 hours from now): %g\n", tm_to_NumDays(&tm, &tm2));
    printf("After - tm2.tm_mon= %d\n", tm2.tm_mon);
    printf("From: %04d-%02d-%02d %02d:%02d:%02d\n",
	(tm.tm_year) +1900,
	(tm.tm_mon) +1,
	tm.tm_mday,
	tm.tm_hour,
	tm.tm_min,
	tm.tm_sec);
    printf("To:   %04d-%02d-%02d %02d:%02d:%02d\n\n",
	(tm2.tm_year) +1900,
	(tm2.tm_mon) +1,
	tm2.tm_mday,
	tm2.tm_hour,
	tm2.tm_min,
	tm2.tm_sec);
    printf("Test of timestring() function:\n");
    printf("%s\n", timestring(0, timebuf));
    printf("%s\n", timestring(1, timebuf));

}

#endif
