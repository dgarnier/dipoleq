/* linkpdfpage.c -- to produce HTML embedded PDF for linking to a given page of
		target PDF without named destination.

 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

cc -Wall -o linkpdfpage -I/usr/local/include linkpdfpage.c -lcpdfpm

1999-08-24  [io]  for v2.00
1998-07-18  [IO]
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cpdflib.h"

int main(int argc, char *argv[])
{
CPDFdoc *pdf;
char *act_dict;
    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);
    cpdf_pageInit(pdf, 1, PORTRAIT, "0 0 144 20", "0 0 144 20");		/* page orientation */

    cpdf_rawRect(pdf, 0.0, 0.0, 144.0, 20.0);
    cpdf_setrgbcolorFill(pdf, 1.0, 1.0, 0.0);
    cpdf_fill(pdf);

    /* == Text example == */
    cpdf_setgrayFill(pdf, 0.0);				/* Black */
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Times-Roman", "MacRomanEncoding", 14.0);
    cpdf_rawTextAligned(pdf, 72.0, 10.0, 0.0, TEXTPOS_MM, "GoToR p.48"); /* cpdf_text() may be repeatedly used */
    cpdf_endText(pdf);


    act_dict = "/Type /Action /S /GoToR /D [47 /FitH 780] /F (../cpdfman110.pdf)";
    cpdf_rawSetLinkAction(pdf, 0.0, 0.0, 144.0, 20.0, act_dict, NULL);

    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "linkpdfpage.pdf");

    /* == Clean up == */
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return(0);
}

