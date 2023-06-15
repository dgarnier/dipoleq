/* shiplabel.c [v1.02] -- This CGI generates a shipping label PDF with bar code
 *                dynamically from HTML form data.
 * Copyright (C) 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

#### NOTE ####
You will need to obtain CGIHTML library to compile and link this CGI.
CGIHTML is avaialble freely from: http://www.eekim.com/software/cgihtml/
On Windows 95/98/NT, you must name the CGI executable as  *.exe, not *.cgi.

FreeBSD/Linux/BSDI (gcc)
gcc -Wall -o shiplabel.cgi shiplabel.c -L. -lcpdfpm -lcgihtml -lm
1999-10-03  [io] for v2.01-r2 or later of ClibPDF
*/

#include <stdio.h>
#include <string.h>
#if defined(_WIN32) || defined(WIN32)
#include <fcntl.h>
extern int setmode(int filenum, int mode);
#endif
#include "cpdflib.h"    /* ClibPDF */
#include "cgi-lib.h"    /* CGIHTML */
#include "string-lib.h" /* CGIHTML */

static  llist entries;  /* for CGIHTML */
void generatePDF(void);

int main(int argc, char *argv[])
{
    read_cgi_input(&entries);  /* parse the form data and add it to the list */
    generatePDF();             /* This is where things happen. */
    list_clear(&entries);      /* free up the pointers in the linked list */
    return(0);
}

void generatePDF(void)
{
CPDFdoc *pdf;
int length;
char *bufPDF;
char linebuf[32], TimeBuf[64], textbuf[1024];
float linespace = 13.0;
float width = 1.2*inch, height = 1.02*inch, xscal = 0.0, yscal = 0.0;
CPDFtboxAttr tba = {TBOX_LEFT, 1, 0.0, 0 }; /* for text box */

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);    /* use Flate/ZIP compression */
    cpdf_init(pdf);
    cpdf_setTitle(pdf, "Shipping Label CGI with Barcode AirBill Number");
    cpdf_pageInit(pdf, 1, PORTRAIT, "0 0 360 288", "0 0 360 288");  /* create page */

    /* pink box at the top */
    cpdf_rect(pdf, 0.0, 2.6, 5.0, 1.4);    /* rectangluar path */
    cpdf_setrgbcolor(pdf, 1.0, 0.9, 0.9);  /* pink */
    cpdf_fill(pdf);          /* fill the rectangle */
    cpdf_setgray(pdf, 0.0);  /* back to black */

    /* image -- you might put a logo here. */
    cpdf_importImage(pdf,  "image1.jpg", JPEG_IMG, 0.4, 2.8, 0.0, &width, &height, &xscal, &yscal, 1);
    cpdf_setActionURL(pdf, 0.4, 2.8, 1.6, 3.82, "http://www.fastio.com/", NULL);

    cpdf_beginText(pdf, 0);	/* Text only from here to cpdf_endText() */

    /* From Address */
    cpdf_setFont(pdf, "Helvetica", "WinAnsiEncoding", 12.0);
    strncpy(textbuf, cgi_val(entries, "From"), 1000);
    strcat(textbuf, "\n");  /* Text Box need to have text terminated by \n. */
    cpdf_textBox(pdf, 2.2, 2.8, 2.4, 1.0, 0.0, linespace, &tba, textbuf);

    /* ShipTo Address */
    cpdf_setFont(pdf, "Helvetica-Bold", "WinAnsiEncoding", 14.0);
    cpdf_textAligned(pdf, 1.2, 2.0, 0.0, TEXTPOS_LR, "SHIP TO:");
    linespace = 15.0;
    strncpy(textbuf, cgi_val(entries, "ShipTo"), 1000);
    strcat(textbuf, "\n");
    cpdf_textBox(pdf, 1.5, 0.9, 3.1, 1.3, 0.0, linespace, &tba, textbuf);

    /* Do tracking number in barcode */
    cpdf_setFont(pdf, "C39HALF_", NULL, 36.0);	/* Code 39 bar code font from www.azalea.com */
    strncpy(linebuf, cgi_val(entries, "TrackingNum"), 12);	/* limit to 12 digits */
    linebuf[12] = '\0';
    sprintf(textbuf, "*%s*", linebuf);	/* Code 39 needs '*' for start/stop */
    cpdf_textAligned(pdf, 3.2, 0.5, 0.0, TEXTPOS_LM, textbuf);	/* bar code */
    cpdf_setFont(pdf, "Courier-Bold", "WinAnsiEncoding", 10.0);
    cpdf_textAligned(pdf, 3.2, 0.28, 0.0, TEXTPOS_UM, textbuf);	/* for human */

    cpdf_setFont(pdf, "Courier-Bold", "WinAnsiEncoding", 8.0); /* for time stamp */
    cpdf_textAligned(pdf, 0.1, 0.12, 0.0, TEXTPOS_LL, timestring(1, TimeBuf));
    cpdf_endText(pdf);	/* text only to this point from cpdf_beginText() */

    /* == generate output  == */
    cpdf_finalizeAll(pdf);
    bufPDF = cpdf_getBufferForPDF(pdf, &length);        /* PDF is in memory buffer */
#if defined(_WIN32) || defined(WIN32)
    setmode(fileno(stdout), O_BINARY);                  /* for binary stdout in Windows */
#endif
    printf("Content-Type: application/pdf%c", 10);      /* correct MIME type for HTTP */
    printf("Content-Length: %d%c%c", length, 10, 10);   /* size of PDF */
    fwrite((void *)bufPDF, 1, (size_t)length, stdout);  /* send PDF now */
    fflush(stdout);

    /* == Clean up == */
    cpdf_close(pdf);			/* shut down ClibPDF */
}

