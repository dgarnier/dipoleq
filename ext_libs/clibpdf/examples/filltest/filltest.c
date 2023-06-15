/* filltest.c -- non-zero winding rule for fill
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

cc -Wall -o filltest filltest.c -lcpdfm

1999-08-22 [io] for v2.00
1998-09-05 [io]
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cpdflib.h>

int main(int argc, char *argv[])
{
CPDFdoc *pdf;
float r1=0.3, r2 = 0.6;
float xorig = 1.5, yorig = 1.5;
float psize = 6.0;

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);
    cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */
    cpdf_translate(pdf, 0.0, 6.5);


    /* Demonstration for non-zero winding rule for fill operator */
    xorig = 1.5;
    yorig = 1.5;
    cpdf_setgrayStroke(pdf, 0.0);				/* Black */
    cpdf_setrgbcolorFill(pdf, 0.7, 0.7, 1.0);
    cpdf_arc(pdf, xorig, yorig, r1, 0.0, 360.0, 1); cpdf_closepath(pdf);  /* CCW */
    cpdf_arc(pdf, xorig, yorig, r2, 0.0, 360.0, 1); cpdf_closepath(pdf);  /* CCW */
    cpdf_fillAndStroke(pdf);
    cpdf_setgrayFill(pdf, 0.0);
    cpdf_pointer(pdf, xorig+r1, yorig, PTR_UP, psize);
    cpdf_pointer(pdf, xorig-r1, yorig, PTR_DOWN, psize);
    cpdf_pointer(pdf, xorig+r2, yorig, PTR_UP, psize);
    cpdf_pointer(pdf, xorig-r2, yorig, PTR_DOWN, psize);

    xorig = 3.0;
    yorig = 1.5;
    cpdf_setrgbcolorFill(pdf, 0.7, 0.7, 1.0);
    cpdf_arc(pdf, xorig, yorig, r1, 360.0, 0.0, 1); cpdf_closepath(pdf);  /* CW */
    cpdf_arc(pdf, xorig, yorig, r2, 0.0, 360.0, 1); cpdf_closepath(pdf);  /* CCW */
    cpdf_fillAndStroke(pdf);
    cpdf_setgrayFill(pdf, 0.0);
    cpdf_pointer(pdf, xorig+r1, yorig, PTR_DOWN, psize);
    cpdf_pointer(pdf, xorig-r1, yorig, PTR_UP, psize);
    cpdf_pointer(pdf, xorig+r2, yorig, PTR_UP, psize);
    cpdf_pointer(pdf, xorig-r2, yorig, PTR_DOWN, psize);


    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "filltest.pdf");

    /* == Clean up == */
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return(0);
}

