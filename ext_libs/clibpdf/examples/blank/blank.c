/* blank.c -- creates a blank PDF file.
 * Copyright (C) 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

cc -Wall -o blank -I/usr/local/include blank.c -lcpdfpm

*/

#include <stdio.h>
#include <stdlib.h>
#include <cpdflib.h>


int main(int argc, char *argv[])
{
CPDFdoc *pdf;

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    //cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);
    cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */

    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "blank.pdf");

    /* == Clean up == */
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return 0;
}
