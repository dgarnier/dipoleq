/* pdfclock.c -- This CGI generates PDF with analog+digital clock based on time of day.
 * Copyright (C) 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

NEXTSTEP/OPENSTEP
cc -Wall -s -object -DNEXTSTEP -o pdfclock.cgi pdfclock.c -I/usr/local/include -lcpdfm
cc -Wall -DNEXTSTEP -o pdfclock.cgi pdfclock.c -I/usr/local/include -lcpdfm

FreeBSD
BSDI/gcc
gcc -Wall -o pdfclock.cgi pdfclock.c -L. -lcpdfm -lm

SunOS5.4/CC-4.0
cc -o pdfclock.cgi pdfclock.c -Xa -L. -lcpdfm -lm

HP-UX B.11.00 (assuming header and lib are installed, or in examples directory below):
cc -O -Aa -Ae +Z +DA1.0 -DHPUX  -I /usr/local/include -I.. -L.. -o pdfclock pdfclock.c -lcpdfm -lm

1999-08-24  [io] for v2.00, also made it runnable in shell (non-CGI)
1999-02-07  [IO]
---
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#if defined(_WIN32) || defined(WIN32)
#include <fcntl.h>
extern int setmode(int filenum, int mode);
#endif
#include "cpdflib.h"

#define  DEGPERHOUR	30.0
#define	 DEGPERMIN	6.0
#define  DEGPERSEC	6.0
#define  XORIG		4.25
#define  YORIG		7.0
#define  FACERADIUS	2.5
#define  HOURHLENGTH	1.8
#define  MINHLENGTH	2.3
#define  SECHLENGTH	2.4
#define  HOURHWIDTH	16.0
#define  MINHWIDTH	12.0
#define  SECHWIDTH	4.0
#define  DEG2RADIAN	0.0174532925

void drawClockFace(CPDFdoc *pdf);
void drawHand(CPDFdoc *pdf, float angle, float length, float width, float r, float g, float b);

int main(int argc, char *argv[])
{
CPDFdoc *pdf;
int length;
char *bufPDF;
time_t t;
struct tm tm; 
float hourangle, minangle, secangle;
char strbuf[128];

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);
    cpdf_setTitle(pdf, "PDF clock CGI");
    cpdf_setSubject(pdf, "PDF clock CGI - generates PDF with analog and digital clock");
    cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);	/* create page */

    /* == PDF drawing here === */
 
    drawClockFace(pdf);

    /* Get current time and compute angles of clock hands */
    t = time(NULL);
    strcpy(strbuf, ctime(&t));		/* for text version */
    strbuf[24] = 0;
    memcpy(&tm, localtime(&t), sizeof(struct tm));
    hourangle = 90.0 - DEGPERHOUR * ((float)(tm.tm_hour % 12) + ((float)tm.tm_min)/60.0);
    minangle = 90.0 - DEGPERMIN * ((float)tm.tm_min + ((float)tm.tm_sec)/60.0);
    secangle = 90.0 - DEGPERSEC * tm.tm_sec;
    if(hourangle < -180.0) hourangle += 360.0;
    if(minangle < -180.0) minangle += 360.0;
    if(secangle < -180.0) secangle += 360.0;

    /* draw clock hands */
    cpdf_setlinecap(pdf, 1);			/* round cap for hands */
    drawHand(pdf, hourangle, HOURHLENGTH, HOURHWIDTH, 0.5, 0.2, 0.0);	/* hour */
    drawHand(pdf, minangle, MINHLENGTH, MINHWIDTH, 0.5, 0.2, 0.0);		/* minute */
    drawHand(pdf, secangle, SECHLENGTH, SECHWIDTH, 0.3, 0.3, 1.0);		/* second */
    cpdf_setlinecap(pdf, 0);

    /* Show Time/Date in text */
    cpdf_setrgbcolor(pdf, 1.0, 0.0, 0.0);
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Times-Roman", "MacRomanEncoding", 34.0);
    cpdf_textAligned(pdf, XORIG, 9.85, 0.0, TEXTPOS_LM, strbuf);
    cpdf_endText(pdf);

    /* some notes */ 
    cpdf_setgray(pdf, 0.0);
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Times-Italic", "MacRomanEncoding", 24.0);
    cpdf_textAligned(pdf, XORIG, 4.0, 0.0, TEXTPOS_LM, "[California Time]");
    cpdf_textAligned(pdf, XORIG, 3.55, 0.0, TEXTPOS_LM, "Dynamic PDF CGI Demo");

    if(getenv("REMOTE_HOST")) {
        strncpy(strbuf, getenv("REMOTE_HOST"), 60);
        cpdf_textAligned(pdf, XORIG, 3.1, 0.0, TEXTPOS_LM, strbuf);
    }
    else
        cpdf_textAligned(pdf, XORIG, 3.1, 0.0, TEXTPOS_LM, "Remote hostname when run as CGI");

    if(getenv("REMOTE_ADDR")) {
        strncpy(strbuf, getenv("REMOTE_ADDR"), 60);
        cpdf_textAligned(pdf, XORIG, 2.65, 0.0, TEXTPOS_LM, strbuf);
    }
    else
        cpdf_textAligned(pdf, XORIG, 2.65, 0.0, TEXTPOS_LM, "Remote host IP when run as CGI");

    cpdf_endText(pdf);

    /* == generate output and clean up == */
    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */

    /* Write to stdout with proper content specification */
    bufPDF = cpdf_getBufferForPDF(pdf, &length);
#if defined(_WIN32) || defined(WIN32)
    setmode(fileno(stdout), O_BINARY);                  /* for binary stdout in Windows */
#endif
    if(getenv("REMOTE_ADDR")) {		/* do this only when running as a CGI */
	printf("Content-Type: application/pdf%c", 10);
	printf("Content-Length: %d%c%c", length, 10, 10);
    }
    fwrite((void *)bufPDF, 1, (size_t)length, stdout);	/* Send PDF now */
    fflush(stdout);

    /* == Clean up == */
    cpdf_close(pdf);			/* shut down */
    return(0);
}

void drawHand(CPDFdoc *pdf, float angle, float length, float width, float r, float g, float b)
{
float xh, yh;
    cpdf_setrgbcolor(pdf, r, g, b);
    cpdf_setlinewidth(pdf, width);
    xh = length * cos(DEG2RADIAN*angle);
    yh = length * sin(DEG2RADIAN*angle);
    cpdf_moveto(pdf, XORIG, YORIG);
    cpdf_lineto(pdf, xh+XORIG, yh+YORIG);
    cpdf_stroke(pdf);
}

void drawClockFace(CPDFdoc *pdf)
{
int i;
float xb, yb, xe, ye;
    cpdf_setlinewidth(pdf, 5.0);
    cpdf_setgray(pdf, 0.0);
    cpdf_setrgbcolorFill(pdf, 0.85, 0.93, 0.7);
    cpdf_circle(pdf, XORIG, YORIG, FACERADIUS);
    cpdf_fillAndStroke(pdf);

    /* hour ticks */
    for(i=0; i<12; i++) {
	xe = FACERADIUS * cos(DEGPERHOUR*DEG2RADIAN*(float)i);
	ye = FACERADIUS * sin(DEGPERHOUR*DEG2RADIAN*(float)i);
	xb = xe * MINHLENGTH / FACERADIUS;
	yb = ye * MINHLENGTH / FACERADIUS;
	cpdf_moveto(pdf, xb+XORIG, yb+YORIG);
	cpdf_lineto(pdf, xe+XORIG, ye+YORIG);
    }
    cpdf_stroke(pdf);

    /* min/second ticks */
    cpdf_setlinewidth(pdf, 2.0);
    for(i=0; i<60; i++) {
	xe = FACERADIUS * cos(DEGPERSEC*DEG2RADIAN*(float)i);
	ye = FACERADIUS * sin(DEGPERSEC*DEG2RADIAN*(float)i);
	xb = xe * SECHLENGTH / FACERADIUS;
	yb = ye * SECHLENGTH / FACERADIUS;
	cpdf_moveto(pdf, xb+XORIG, yb+YORIG);
	cpdf_lineto(pdf, xe+XORIG, ye+YORIG);
    }
    cpdf_stroke(pdf);
}

