/* MarkerTest.c -- Draws all marker symbols, filled and open.
		   Also draws pointer arrows.
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

cc -Wall -o MarkerTest -I/usr/local/include MarkerTest.c -lcpdfpm

1999-08-24  [io] for v2.00
1998-08-16  [IO]
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cpdflib.h"

#define NUMMARKERS	11
#define XBUMP		0.35

int main(int argc, char *argv[])
{
CPDFdoc *pdf;
float x, y1, y2, y0, x0 = 1.0;
char strbuf[128];
int i;
float markersize = 12.0;
	pdf = cpdf_open(0, NULL);
	cpdf_enableCompression(pdf, YES);
	cpdf_setCreator(pdf, "ClibPDF Demo Program: MarkerTest.c");
	cpdf_init(pdf);
	cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */
        cpdf_translate(pdf, 0.0, 6.5);

	/* Your drawing code in PDF below */
	/* ----------------------------------------------------------------------------------------- */

	/* draw dashed lines to define reference points at which markers are positioned. */
	cpdf_setdash(pdf, "[2 2] 0");
	cpdf_setlinewidth(pdf, 0.12);
	y1 = 1.3;
	cpdf_moveto(pdf, 1.1, y1);
	cpdf_lineto(pdf, 6.9, y1);
	y2 = 1.6;
	cpdf_moveto(pdf, 1.1, y2);
	cpdf_lineto(pdf, 6.9, y2);

	y0 = 1.45;
	for(i=0; i< (NUMMARKERS + 4); i++) {
	    if(i < NUMMARKERS)
	        x = XBUMP * (float)i + XBUMP + x0;
	    else
	        x = XBUMP * (float)i + XBUMP*2.0 + x0;
	    cpdf_moveto(pdf, x, y0+0.45);
	    cpdf_lineto(pdf, x, y0-0.45);
	}
	cpdf_stroke(pdf);
	cpdf_nodash(pdf);

	/* Draw markers and pointers at intersections defined by dashed lines */
	/* First do filled (black) symbols */
	cpdf_setgrayStroke(pdf, 0.0);
	cpdf_setgrayFill(pdf, 0.0);
	for(i=0; i< (NUMMARKERS + 4); i++) {
	    if(i < NUMMARKERS) {
	        x = XBUMP * (float)i + XBUMP + x0;
	        cpdf_marker(pdf, x, y1, i, markersize);	/* filled */
	    }
	    else {
	        x = XBUMP * (float)i + XBUMP*2 + x0;
	        cpdf_pointer(pdf, x, y1, i-NUMMARKERS, markersize);	/* filled */
	    }
	}

	/* Do open symbols */
	cpdf_setgrayStroke(pdf, 0.0);
	cpdf_setgrayFill(pdf, 1.0);
	for(i=0; i< (NUMMARKERS + 4); i++) {
	    if(i < NUMMARKERS) {
	        x = XBUMP * (float)i + XBUMP + x0;
	        cpdf_marker(pdf, x, y2, i, markersize);	/* open */
	    }
	    else {
	        x = XBUMP * (float)i + XBUMP*2 + x0;
	        cpdf_pointer(pdf, x, y2, i-NUMMARKERS, markersize);	/* open */
	    }
	}

	/* Show symbol number above corresponding symbols */
	cpdf_setgrayFill(pdf, 0.0);		/* black */
	cpdf_beginText(pdf, 0);
	cpdf_setFont(pdf, "Helvetica", "MacRomanEncoding", 12.0);
	for(i=0; i< (NUMMARKERS + 4); i++) {
	    if(i < NUMMARKERS) {
	        x = XBUMP * (float)i + XBUMP + x0;
		if(i == 10) x -= 0.06;
	        sprintf(strbuf, "%d", i);
	    }
	    else {
	        x = XBUMP * (float)i + XBUMP*2 + x0;
	        sprintf(strbuf, "%d", i-NUMMARKERS);
	    }
	    cpdf_text(pdf, x-0.05, y0+0.52, 0.0, strbuf);
	}
	cpdf_endText(pdf);

	cpdf_beginText(pdf, 0);
	cpdf_setFont(pdf, "Times-Roman", "MacRomanEncoding", 14.0);
	cpdf_text(pdf, 1.2, 2.6, 0.0, "MarkerTest.c  --   markersize = 12 points");
	cpdf_text(pdf, 2.55, 2.24, 0.0, "cpdf_marker()");
	cpdf_text(pdf, 5.48, 2.24, 0.0, "cpdf_pointer()");
	cpdf_endText(pdf);

	/* ----------------------------------------------------------------------------------------- */
	/* Your drawing code in PDF above */
        cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
	cpdf_savePDFmemoryStreamToFile(pdf, "markers.pdf");

	cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
	cpdf_close(pdf);			/* shut down */
	return(0);
}

