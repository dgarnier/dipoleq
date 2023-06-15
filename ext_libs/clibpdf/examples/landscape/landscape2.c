/* landscape2.c -- How to do LANDSCAPE mode page with 90 deg rotation.
 * Copyright (C) 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

This example rotates the page the other way from the one in landscape.c.
This ALTLANDSCAPE mode makes the bottom-right corner of the original PORTRAIT
page the new origin.

cc -Wall -o landscape2 -I/usr/local/include landscape2.c -lcpdfpm

*/

#include <stdio.h>
#include <stdlib.h>
#include <cpdflib.h>

int main(int argc, char *argv[])
{
CPDFdoc *pdf;

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);

    /*  The following 3 lines will cause the page to be in ALTLANDSCAPE orientation.
	"ALTLANDSCAPE" as the third argument to cpdf_pageInit() essentially instructs
	the PDF viewer (Acrobat Reader) to display the page in LANDSCAPE mode.
	However, it does not change what goes on to the page -- the short dimension
	of the sheet is still the X axis.  To change that, you must actually change
	the coordinate systems by explicitly translating the origin and rotating the
	page.  For non-LETTER size page, you must change the 2-nd argument in
	cpdf_translate() to the length of the short dimension of the page.

    */
    cpdf_pageInit(pdf, 1, ALTLANDSCAPE, LETTER, LETTER);		/* page orientation */
    cpdf_translate(pdf, 8.5, 0.0);	/* 8.5 must be changed for non-LETTER size */
    cpdf_rotate(pdf, 90.0);

    cpdf_rect(pdf, 1.0, 1.5, 5.0, 2.5);
    cpdf_stroke(pdf);

    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Times-Italic", "WinAnsiEncoding", 48.0);
    cpdf_text(pdf, 1.0, 3.0, 0.0, "Hello Euro \200");
    cpdf_endText(pdf);

    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "landscape2.pdf");

    /* == Clean up == */
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return 0;
}
