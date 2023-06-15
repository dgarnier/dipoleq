/* charwidth.c -- Test program for all width data for char codes > 128
		  and for font encoding look-up tables.
 * Copyright (C) 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

   Right-aligned text will reveal most problems of wrong encoding or width data.

cc -Wall -o charwidth -I/usr/local/include charwidth.c -lcpdfpm

*/

#include <stdio.h>
#include <stdlib.h>
#include <cpdflib.h>

int main(int argc, char *argv[])
{
CPDFdoc *pdf;
float x = 7.0;
float y = 9.0;
float ybump = 0.5;
int i, j, ibase=128;
unsigned char strbuf[20];
char *fontname = "AvantGarde-BookOblique";

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);
    cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */

    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, fontname, "WinAnsiEncoding", 32.0);
    cpdf_textAligned(pdf, x, y, 0.0, TEXTPOS_LR, fontname);
     
    for(j=12; j<26; j++) {
	ibase = j*10;
	for(i=0; i<10; i++)
	    strbuf[i] = i + ibase;
	strbuf[10] = 0;
        y -= ybump;
        cpdf_textAligned(pdf, x, y, 0.0, TEXTPOS_LR, strbuf);	/* 10 chars at a time */
    }

    cpdf_endText(pdf);

    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "charwidth.pdf");

    /* == Clean up == */
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return 0;
}
