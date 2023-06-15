/* beziertest.c -- Test program for dashes
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

cc -Wall -o beziertest beziertest.c -lcpdfm

1998-09-12  [IO]
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cpdflib.h>


int main(int argc, char *argv[])
{
CPDFdoc *pdf;
float x0 = 1.0, y0 = 1.0;
float x1 = 1.2, y1 = 1.5;
float x2 = 2.8, y2 = 2.2;
float x3 = 3.0, y3 = 1.2;
float fontsize = 12.0;
float markersize = 4.0;

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);
    cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */
    cpdf_translate(pdf, 0.0, 6.5);

    /* First Bezier curve example */
    cpdf_setgray(pdf, 0.0);
    cpdf_setlinewidth(pdf, 1.5);
    cpdf_moveto(pdf, x0, y0);
    cpdf_curveto(pdf, x1, y1, x2, y2, x3, y3);
    cpdf_stroke(pdf);

    cpdf_setdash(pdf, "[2 3] 0");
    cpdf_setlinewidth(pdf, 0.3);
    cpdf_moveto(pdf, x0, y0);
    cpdf_lineto(pdf, x1, y1);
    cpdf_moveto(pdf, x3, y3);
    cpdf_lineto(pdf, x2, y2);
    cpdf_stroke(pdf);
    cpdf_nodash(pdf);

    cpdf_marker(pdf, x1, y1, 0, markersize);		/* control points */
    cpdf_marker(pdf, x2, y2, 0, markersize);

    cpdf_setgrayFill(pdf, 1.0);
    cpdf_marker(pdf, x0, y0, 3, markersize);		/* end points */
    cpdf_marker(pdf, x3, y3, 3, markersize);

    cpdf_setgrayFill(pdf, 0.0);
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Helvetica-Oblique", "MacRomanEncoding", fontsize);
    cpdf_text(pdf, x0-0.15, y0-0.25, 0.0, "x0, y0");	/* cpdf_text() may be repeatedly used */
    cpdf_text(pdf, x3-0.15, y3-0.25, 0.0, "x3, y3");	/* cpdf_text() may be repeatedly used */
    cpdf_text(pdf, x1-0.15, y1+0.1, 0.0, "x1, y1");	/* cpdf_text() may be repeatedly used */
    cpdf_text(pdf, x2-0.15, y2+0.1, 0.0, "x2, y2");	/* cpdf_text() may be repeatedly used */
    cpdf_endText(pdf);

    /* The second Bezier example */
    x0 = 4.0; y0 = 1.0;
    x1 = 4.3; y1 = 2.2;
    x2 = 5.8; y2 = 0.7;
    x3 = 6.0; y3 = 1.7;

    cpdf_setgray(pdf, 0.0);
    cpdf_nodash(pdf);
    cpdf_setlinewidth(pdf, 1.5);
    cpdf_moveto(pdf, x0, y0);
    cpdf_curveto(pdf, x1, y1, x2, y2, x3, y3);
    cpdf_stroke(pdf);

    cpdf_setdash(pdf, "[2 3] 0");
    cpdf_setlinewidth(pdf, 0.3);
    cpdf_moveto(pdf, x0, y0);
    cpdf_lineto(pdf, x1, y1);
    cpdf_moveto(pdf, x3, y3);
    cpdf_lineto(pdf, x2, y2);
    cpdf_stroke(pdf);
    cpdf_nodash(pdf);

    cpdf_marker(pdf, x1, y1, 0, markersize);		/* control points */
    cpdf_marker(pdf, x2, y2, 0, markersize);

    cpdf_setgrayFill(pdf, 1.0);
    cpdf_marker(pdf, x0, y0, 3, markersize);		/* end points */
    cpdf_marker(pdf, x3, y3, 3, markersize);

    cpdf_setgrayFill(pdf, 0.0);
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Helvetica-Oblique", "MacRomanEncoding", fontsize);
    cpdf_text(pdf, x0-0.15, y0-0.25, 0.0, "x0, y0");	/* cpdf_text() may be repeatedly used */
    cpdf_text(pdf, x1-0.15, y1+0.1, 0.0, "x1, y1");	/* cpdf_text() may be repeatedly used */
    cpdf_text(pdf, x2-0.15, y2-0.25, 0.0, "x2, y2");	/* cpdf_text() may be repeatedly used */
    cpdf_text(pdf, x3-0.15, y3+0.1, 0.0, "x3, y3");	/* cpdf_text() may be repeatedly used */
    cpdf_endText(pdf);

    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "bezier.pdf");

    /* == Clean up == */
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return(0);
}

