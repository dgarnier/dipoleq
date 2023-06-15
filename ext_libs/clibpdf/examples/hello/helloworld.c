/* helloworld.c -- For a bare minimum example.
 * Copyright (C) 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

cc -Wall -o helloworld -I/usr/local/include helloworld.c -lcpdfpm

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
    cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */

    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Times-Italic", "WinAnsiEncoding", 48.0);
    cpdf_text(pdf, 1.0, 8.0, 0.0, "Hello Euro \200");
    cpdf_endText(pdf);

    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "helloeuro.pdf");

    /* == Clean up == */
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return 0;
}
