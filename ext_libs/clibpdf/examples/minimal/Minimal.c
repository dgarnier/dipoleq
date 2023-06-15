/* Minimal.c -- A simple drawing and one-line text example.
		It can't be really simpler than this without being silly.

 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

cc -Wall -o Minimal -I/usr/local/include Minimal.c -lcpdfpm

HP-UX B.11.00 (assuming header and lib are installed, or in examples directory below):
cc -O -Aa -Ae +Z +DA1.0 -DHPUX  -I /usr/local/include -I.. -L.. -o Minimal Minimal.c -lcpdfp -lm

1999-08-24  [io] for v2.00
1998-07-18  [IO]
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cpdflib.h"

int main(int argc, char *argv[])
{
CPDFdoc *pdf;
int i;
float x, y, angle;
float radius2= 1.2;
float xorig = 2.5, yorig = 2.5;

    if(argc != 2) {
	fprintf(stderr, "Usage: %s file.pdf\n    Use '-' as argument for stdout output.\n", argv[0]);
	exit(0);
    }

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);
    cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */
    cpdf_translate(pdf, 0.0, 6.0);

    /* == Simple graphics drawing example == */
    for(i=0; i<=200; i++) {
	angle = PI*(float)i/100.0;		/* 0 .. 2pi */
	x = radius2 * cos(3.0*angle) + xorig;
	y = radius2 * sin(5.0*angle) + yorig;
	if(i) cpdf_lineto(pdf, x, y);
	else  cpdf_moveto(pdf, x, y);		/* first point */
    }
    cpdf_closepath(pdf);
    cpdf_setrgbcolorFill(pdf, 1.0, 0.7, 0.7);	/* pink as fill color */
    cpdf_setrgbcolorStroke(pdf, 0.0, 0.0, 1.0);	/* blue as stroke color */
    cpdf_eofillAndStroke(pdf);			/* even-odd fill and stroke*/

    /* == Text example == */
    cpdf_setgrayFill(pdf, 0.0);				/* Black */
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Times-Italic", "MacRomanEncoding", 16.0);
    cpdf_text(pdf, 1.6, 1.0, 0.0, "x=cos(3t), y=sin(5t)");	/* cpdf_text() may be repeatedly used */
    cpdf_endText(pdf);

    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, argv[1]);

    /* == Clean up == */
    if(strcmp(argv[1], "-") != 0)
        cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return(0);
}

