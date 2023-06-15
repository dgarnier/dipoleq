/* text2pdf.c -- Yet another TEXT2PDF program.
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

This example is nearly identical to examples/textbox/textbox.c.
It just shows the flexibility of cpdf_textBox() function.

OPENSTEP
cc -Wall -DUNIX -o text2pdf text2pdf.c -I/usr/local/include -lcpdfpm
NEXTSTEP
cc -Wall -DUNIX -s -object -o text2pdf text2pdf.c -lcpdfpm

BSDI/FreeBSD/gcc
gcc -Wall -DUNIX -o text2pdf text2pdf.c -L/usr/local/lib -lcpdfm -lm

HP-UX B.11.00 (assuming header and lib are installed, or in examples directory below):
cc -O -Aa -Ae +Z +DA1.0 -DHPUX -DUNIX -I /usr/local/include -I.. -L.. -o text2pdf text2pdf.c -lcpdfm -lm

1999-09-14 [io] for v2.00
1999-06-05 [IO] First version.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cpdflib.h"

#define NUMCOLS		1	/* === Number of text columns per page (1 or 2) */


#define LLM		72.0
#define RLM		315.0
#define BM		60.0
#define COLUMN1WIDTH	468.0
#define COLUMN2WIDTH	225.0
#define COLUMNHEIGHT	690.0

int main(int argc, char *argv[])
{
CPDFdoc *pdf;
FILE *fp;
char *textdata, *currtext;
long filesize, nbread;
int page =1;
int nColumns = NUMCOLS;
float wColumn;
float fontsize = 12.0;
float linespace = -1.1;	/* scale line space with font size at this factor */
CPDFtboxAttr tboxatr;
char numbuf[64];

    if(argc != 2) {
	fprintf(stderr, "Usage: text2pdf textfile\n");
	exit(1);
    }

    if(nColumns == 2)  wColumn = COLUMN2WIDTH;
    else               wColumn = COLUMN1WIDTH;

    filesize = getFileSize(argv[1]);
    if(filesize <= 0) {
	fprintf(stderr, "Zero-size file or no such file: %s\n", argv[1]);
	exit(1);
    }
    /* Open the text file in its respective native format, sigh. */
    if((fp = fopen(argv[1], "r")) == NULL) {
	fprintf(stderr, "Can't open file: %s\n", argv[1]);
	exit(1);
    }

    /* Get the whole text from file into a memory buffer.  Memory map if you have decent VM. */
    textdata = (char *)malloc(filesize+1024);
    if(textdata == NULL) {
	fprintf(stderr, "Can't malloc() memory\n");
	exit(1);
    }
    nbread = fread((void *)textdata, 1, (size_t)filesize, fp);
    fclose(fp);
    textdata[nbread] = '\0';	/* terminate as string */
    currtext = textdata;	/* initialize pointer */

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_setOutputFilename(pdf, "text2pdf.pdf");
    cpdf_init(pdf);

    /* PDF Info object */
    cpdf_setCreator(pdf, "text2pdf.c");
    cpdf_setTitle(pdf, "Example for cpdf_textBox() function");
    cpdf_setSubject(pdf, "ANSI C Library for Direct PDF Generation");
    cpdf_setKeywords(pdf, "ClibPDF, FastIO Systems, ANSI C Library, Acrobat, PDF, text2pdf, Dynamic PDF, CGI");

    tboxatr.alignmode = TBOX_LEFT;
    tboxatr.NLmode = 1;			/* honor new lines */
    tboxatr.paragraphSpacing = 0.0;
    tboxatr.noMark = 0;			/* reformat */

    /* ====== MAIN TEXT2PDF LOOP ===================================================== */
    while(currtext) {
        cpdf_pageInit(pdf, page, PORTRAIT, LETTER, LETTER);			/* create page */
	cpdf_beginText(pdf, 0);
	cpdf_setFont(pdf, "CPDF-Monospace", "MacRomanEncoding", fontsize);
	/* left column */
	currtext = cpdf_rawTextBox(pdf, LLM, BM, wColumn, COLUMNHEIGHT, 0.0, linespace, &tboxatr, currtext);
	/* right column */
	if(currtext && nColumns > 1)
	    currtext = cpdf_rawTextBox(pdf, RLM, BM, wColumn, COLUMNHEIGHT, 0.0, linespace, &tboxatr, currtext);
	sprintf(numbuf, "page %d", page);
	cpdf_rawTextAligned(pdf,  540.0, 40.0, 0.0, TEXTPOS_LR, numbuf);	/* page number */
        cpdf_endText(pdf);
	cpdf_finalizePage(pdf, page);					/* close page */
	page++;
    }
    /* ====== That's it for multi-column, multi-page text formatting. Cool? === */

    /* === All drawing done.  Finalize the content === */
    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    /* cpdf_savePDFmemoryStreamToFile(pdf, "text2pdf.pdf"); */

    /* == Clean up == */
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return(0);
}


