/* barcode.c -- demonstrates uses of a bar code font.
 * Copyright (C) 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

cc -Wall -o barcode -I/usr/local/include barcode.c -lcpdfpm

*/

#include <stdio.h>
#include <stdlib.h>
#include <cpdflib.h>

#define WEBSITE_STR	"http://www.azalea.com/"

int main(int argc, char *argv[])
{
CPDFdoc *pdf;
float barcodefontsize = 72.0;
float fontsize = 14.0;
float y = 7.0, strwidth1;
float ybump = 0.25;
char *barstring = "*452486646698*";

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);
    cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */

    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "C39ONE__", NULL, barcodefontsize);	/* bar code font */
    cpdf_textAligned(pdf, 4.25, y, 0.0, TEXTPOS_LM, barstring);
    cpdf_endText(pdf);

    /* human readable version */
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Ocrbcr__", NULL, fontsize);		/* OCR-B font */
    cpdf_textAligned(pdf, 4.25, y-0.5, 0.0, TEXTPOS_LM, barstring);
    cpdf_endText(pdf);

    y -= 1.0;
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Helvetica", "MacRomanEncoding", fontsize);
    cpdf_textAligned(pdf, 4.25, y, 0.0, TEXTPOS_LM,
	"This example uses a Type-1 bar code font (C39) from Azalea Software.");
    y -= ybump;
    cpdf_textAligned(pdf, 4.25, y, 0.0, TEXTPOS_LM, WEBSITE_STR);
    strwidth1 = cpdf_stringWidth(pdf, WEBSITE_STR)/72.0;
    cpdf_endText(pdf);

    /* Add Web and email links over the address text above. */
    cpdf_setActionURL(pdf, 4.25 - strwidth1/2.0, y - 0.05,
			 4.25 + strwidth1/2.0, y + 0.2,
			 WEBSITE_STR, NULL);

    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "barcode0.pdf");

    /* == Clean up == */
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return 0;
}
