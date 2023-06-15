/* timeaxis.c -- Test program for time axis
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

cc -Wall -o timeaxis timeaxis.c -lcpdfpm

1999-08-22  [io] for v2.00
1998-09-12  [IO]
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cpdflib.h>


int main(int argc, char *argv[])
{
CPDFdoc *pdf;
int i;
CPDFplotDomain *oldDomain, *timeDomain;
CPDFaxis  *tAxis;
time_t t;
struct tm tm; 
struct tm tm2;
float respmax = 10.0;
float ybump, yoffset;
int yearformat = YEAR_2DIGIT;

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);
    cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */
    cpdf_translate(pdf, 0.0, 1.0);
    /* ---------------------------------------------------------------------------------- */

    /* Time domain version */
    t = time(NULL);
    /* memcpy(&tm, localtime(&t), sizeof(struct tm)); */
    tm.tm_sec = 0;
    tm.tm_min = 0;
    tm.tm_hour = 0;
    tm.tm_mday = 1;		/* day starts at 1 */
    tm.tm_mon = 0;		/* month starts at 0 */
    tm.tm_year = 98;		/* year since 1900 */
    tm.tm_isdst = -1;		/* This is needed. */
    mktime(&tm);		/* this rebuilds all the fields */
    ybump = 0.5*inch;

    /* -- start the first axis with 4 hour span */
    memcpy(&tm2, &tm, sizeof(struct tm));
    tm2.tm_hour += 4;		/* 4 hour */
    mktime(&tm2);		/* this rebuilds all the fields */
    yoffset = -0.2*inch;

    /* Creat time plot domain and use previously created time axis */
    timeDomain = cpdf_createTimePlotDomain( 1.0*inch, 8.0*inch, 4.0*inch, 0.5*inch,
			    &tm, &tm2, 0.0, respmax, TIME, LINEAR, 0);
    oldDomain = cpdf_setPlotDomain(pdf, timeDomain);	/* save oldDomain for later restore */
    cpdf_fillDomainWithGray(timeDomain, 0.85);				/* white */
    /* cpdf_drawMeshForDomain(timeDomain); */

    cpdf_setgray(pdf, 0.0);
    tAxis = cpdf_createTimeAxis(0.0, 4.0*inch, TIME, &tm, &tm2);
    cpdf_attachAxisToDomain(tAxis, timeDomain, 0.0, yoffset );
    cpdf_setTimeAxisNumberFormat(tAxis, MONTH_NAME, YEAR_2DIGIT, "Helvetica", 11.0);	
	    /* YEAR_FULL for 1998 or YEAR_2DIGIT for 98 */
    /* cpdf_setAxisLabel(tAxis, "Date", "Helvetica-BoldOblique", "MacRomanEncoding", 14.0); */
    cpdf_drawAxis(tAxis);
    cpdf_freeAxis(tAxis);
    /* -- First axis done */

    /* -- Other axes -- */
    for(i=0; i<4; i++) {
	memcpy(&tm2, &tm, sizeof(struct tm));	/* freshen tm2 */
	switch(i) {
		default:
		case 0:
		    tm2.tm_hour += 24;		/* 24 hours */
		    break;
		case 1:
		    tm2.tm_mday += 4;		/* 4 days */
		    break;
		case 2:
		    tm2.tm_mon += 6;		/* 6 months */
		    break;
		case 3:
		    tm2.tm_year += 4;		/* 4 years */
		    yearformat = YEAR_FULL;	/* use 4 digits */
		    break;
	}
	tm2.tm_isdst = -1;		/* This is needed. */
	mktime(&tm2);		/* this rebuilds all the fields */
	yoffset -= ybump;

	/* We reuse the domain, as we don't actually plot anything other than axis */
	cpdf_setgray(pdf, 0.0);
	tAxis = cpdf_createTimeAxis(0.0, 4.0*inch, TIME, &tm, &tm2);
	cpdf_attachAxisToDomain(tAxis, timeDomain, 0.0, yoffset );
	cpdf_setTimeAxisNumberFormat(tAxis, MONTH_NAME, yearformat, "Helvetica", 11.0);	
		/* YEAR_FULL for 1998 or YEAR_2DIGIT for 98 */
	/* cpdf_setAxisLabel(tAxis, "Date", "Helvetica-BoldOblique", "MacRomanEncoding", 14.0); */
	cpdf_drawAxis(tAxis);
	cpdf_freeAxis(tAxis);
	/* -- one axis done */
    }



    cpdf_setPlotDomain(pdf, oldDomain);	/* restore old domain */
    cpdf_freePlotDomain(timeDomain);		/* deallocate the plot domain */

    /* ---------------------------------------------------------------------------------- */
    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "timeaxis.pdf");

    /* == Clean up == */
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return(0);
}

