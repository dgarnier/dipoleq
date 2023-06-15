/* dashtest.c -- Test program for dashes
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

cc -Wall -o dashtest -I/usr/local/include dashtest.c -lcpdfm

1999-08-24 [io] for v2.00
1998-09-11 [io]
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cpdflib.h>

static char *dashpattern[] = {
	"[2  2] 0",
	"[4  4] 0",
	"[8  8] 0",
	"[8  8] 4",
	"[8  8] 8",
	"[12  4] 0",
	"[16  3  4  3] 0",
	"[13  3  2  3  2  3] 0",
	"[ ] 0"
};

#define Ndash	9

int main(int argc, char *argv[])
{
CPDFdoc *pdf;
int i;
float ybump = 0.22;
float x = 1.0, y = 3.0;
float fontsize = 12.0;

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);
    cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */
    cpdf_translate(pdf, 0.0, 6.5);

    cpdf_setgray(pdf, 0.0);
    for(i=0; i< Ndash; i++) {
	cpdf_setdash(pdf, dashpattern[i]);
	cpdf_moveto(pdf, x, y);
	cpdf_lineto(pdf, x+2.0, y);
	cpdf_stroke(pdf);
	y -= ybump;
    }
    
    x = 3.1; y = 3.0 - fontsize*0.3/inch;
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Helvetica", "MacRomanEncoding", fontsize);
    for(i=0; i< Ndash; i++) {
        cpdf_text(pdf, x, y, 0.0, dashpattern[i]);	/* cpdf_text() may be repeatedly used */
	y -= ybump;
    }
    cpdf_endText(pdf);

    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "dashtest.pdf");

    /* == Clean up == */
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return(0);
}

