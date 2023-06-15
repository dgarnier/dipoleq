/* textalign.c -- alignment mode for text for manual
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

cc -Wall -o textalign -I/usr/local/include textalign.c -lcpdfpm

1999-08-22 [io] for v2.00

*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cpdflib.h>

static char *tposstr[] = {
  "TEXTPOS_LL", "TEXTPOS_LM", "TEXTPOS_LR",
  "TEXTPOS_ML", "TEXTPOS_MM", "TEXTPOS_MR",
  "TEXTPOS_UL", "TEXTPOS_UM", "TEXTPOS_UR"
};

int main(int argc, char *argv[])
{
CPDFdoc *pdf;
int i;
float xlb = 0.5, xle = 6.5;
float ylb = 0.5, yle = 2.3;
float x = xlb+0.5;
float y = ylb+0.5;
float xbump = 2.5;
float ybump = 0.4;
float xp, yp;
float fontsize = 14.0;

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);
    cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */
    cpdf_translate(pdf, 0.0, 7.0);

    cpdf_setlinewidth(pdf, 0.15);
    cpdf_setdash(pdf, "[1 2] 0");
    for(i=0; i<3; i++) {
	cpdf_moveto(pdf, xlb, y+ybump*(float)i);
	cpdf_lineto(pdf, xle, y+ybump*(float)i);
    }
    cpdf_stroke(pdf);
    for(i=0; i<3; i++) {
	cpdf_moveto(pdf, x+xbump*(float)i, ylb);
	cpdf_lineto(pdf, x+xbump*(float)i, yle);
    }
    cpdf_stroke(pdf);
    cpdf_nodash(pdf);

    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Helvetica", "MacRomanEncoding", fontsize);
    for(i=0; i< 9; i++) {
	xp = x + (float)(i%3)*xbump;
	yp = y + (float)(i/3)*ybump;
        cpdf_textAligned(pdf, xp, yp, 0.0, i, tposstr[i]);	/* cpdf_text() may be repeatedly used */
    }
    cpdf_endText(pdf);

    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "textalign.pdf");

    /* == Clean up == */
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return(0);
}
