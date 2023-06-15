/* Arcs.c -- Test program for arc and circle functions.
 * Copyright (C) 1998,1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

NEXTSTEP/OPENSTEP
cc -Wall -o Arcs Arcs.c -I/usr/local/include -lcpdfm
cc -Wall -o Arcs Arcs.c -I/usr/local/include -lcpdfpm

HP-UX B.11.00 (assuming header and lib are installed, or in examples directory below):
cc -O -Aa -Ae +Z +DA1.0 -DHPUX  -I /usr/local/include -I.. -L.. -o arctest Arcs.c -lcpdfm -lm
cc -O -Aa -Ae +Z +DA1.0 -DHPUX  -I /usr/local/include -I.. -L.. -o arctest Arcs.c -lcpdfpm -lm

1999-08-24  [io] for ClibPDF 2.00.
1998-09-05  [io]
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cpdflib.h"

int main(int argc, char *argv[])
{
int i;
float radius= 2.2;
float rbump = 0.1;
float xorig = 2.5, yorig = 2.5;
float sangle, eangle;
float greenblue = 1.0;
CPDFdoc *pdf;

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);
    cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */

    /* == Simple graphics drawing example == */
    cpdf_setgrayStroke(pdf, 0.0);	/* black as stroke color */
    cpdf_comments(pdf, "\n%% one circle.\n");
    cpdf_circle(pdf, xorig, yorig, radius);
    cpdf_stroke(pdf);			/* even-odd fill and stroke*/

    /* Just a bunch of blue arcs */
    cpdf_setrgbcolorStroke(pdf, 0.0, 0.0, 1.0);	/* blue as stroke color */
    cpdf_comments(pdf, "\n%% Arcs from 30 to various angles in 30 deg steps.\n");
    for(i=11; i>=0; i--) {
	radius -= rbump;
	eangle = (float)(i+1)*30.0;		/* end angle */
	sangle = 0.0;
	cpdf_arc(pdf, xorig, yorig, radius, sangle, eangle, 1);	/* moveto to starting point of arc */
	cpdf_stroke(pdf);
    }

    yorig = 7.5;
    radius= 2.2;
    cpdf_setgrayStroke(pdf, 0.0);	/* black as stroke color */
    cpdf_comments(pdf, "\n%% one circle again.\n");
    cpdf_circle(pdf, xorig, yorig, radius);
    cpdf_stroke(pdf);			/* even-odd fill and stroke*/

    /* Now, do progressively redder pie shapes from large ones to small */
    cpdf_setrgbcolorStroke(pdf, 0.7, 0.7, 0.0);	/* yellow as stroke color */
    cpdf_comments(pdf, "\n%% Pie-shapes from 30 to various angles in 30 deg steps.\n");
    for(i=11; i>=0; i--) {
	radius -= rbump;
	greenblue = (float)i/12.0;
	cpdf_setrgbcolorFill(pdf, 1.0, greenblue, greenblue);
	eangle = (float)(i+1)*30.0;		/* end angle */
	sangle = 0.0;
	cpdf_moveto(pdf, xorig, yorig);
	cpdf_arc(pdf, xorig, yorig, radius, sangle, eangle, 0);	/* lineto to start of arc */
	cpdf_closepath(pdf);
        cpdf_fillAndStroke(pdf);			/* fill and stroke */
    }

    /* Demonstration for non-zero winding rule for fill operator */
    xorig = 6.0;
    yorig = 2.5;
    cpdf_setgrayStroke(pdf, 0.0);				/* Black */
    cpdf_setrgbcolorFill(pdf, 0.6, 0.6, 1.0);
    cpdf_arc(pdf, xorig, yorig, 0.5, 0.0, 360.0, 1); cpdf_closepath(pdf);  /* CCW */
    cpdf_arc(pdf, xorig, yorig, 1.0, 0.0, 360.0, 1); cpdf_closepath(pdf);  /* CCW */
    cpdf_fillAndStroke(pdf);
    cpdf_setgrayFill(pdf, 0.0);
    cpdf_pointer(pdf, xorig+0.5, yorig, PTR_UP, 8.0);
    cpdf_pointer(pdf, xorig-0.5, yorig, PTR_DOWN, 8.0);
    cpdf_pointer(pdf, xorig+1.0, yorig, PTR_UP, 8.0);
    cpdf_pointer(pdf, xorig-1.0, yorig, PTR_DOWN, 8.0);

    xorig = 6.0;
    yorig = 7.5;
    cpdf_setrgbcolorFill(pdf, 0.6, 0.6, 1.0);
    cpdf_arc(pdf, xorig, yorig, 0.5, 360.0, 0.0, 1); cpdf_closepath(pdf);  /* CW */
    cpdf_arc(pdf, xorig, yorig, 1.0, 0.0, 360.0, 1); cpdf_closepath(pdf);  /* CCW */
    cpdf_fillAndStroke(pdf);
    cpdf_setgrayFill(pdf, 0.0);
    cpdf_pointer(pdf, xorig+0.5, yorig, PTR_DOWN, 8.0);
    cpdf_pointer(pdf, xorig-0.5, yorig, PTR_UP, 8.0);
    cpdf_pointer(pdf, xorig+1.0, yorig, PTR_UP, 8.0);
    cpdf_pointer(pdf, xorig-1.0, yorig, PTR_DOWN, 8.0);

    /* == Text examples == */
    /* cpdf_setgrayFill(pdf, 0.0); */				/* Black */
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Times-Italic", "MacRomanEncoding", 16.0);
    cpdf_text(pdf, 1.6, 2.0, 0.0, "Test of arcs and circles");	/* cpdf_text() may be repeatedly used */
    cpdf_text(pdf, 1.6, 7.0, 0.0, "Color filled pie shapes");	/* cpdf_text() may be repeatedly used */
    cpdf_text(pdf, 4.7, 5.0, 0.0, "Non-zero winding rule for fill");	/* cpdf_text() may be repeatedly used */
    cpdf_endText(pdf);

    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "arctest.pdf");

    /* == Clean up == */
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return 0;
}

