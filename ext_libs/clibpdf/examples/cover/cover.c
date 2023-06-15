/* cover.c -- Program to generate a cover for the ClibPDF API Reference Manual
 * Copyright (C) 1998,1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

OPENSTEP
cc -Wall -o cover cover.c -I/usr/local/include -lcpdfm
NEXTSTEP
cc -Wall -s -object -o cover cover.c -lcpdfpm

BSDI/FreeBSD/gcc
gcc -Wall -o cover cover.c -L/usr/local/lib -lcpdfm -lm

HP-UX B.11.00 (assuming header and lib are installed, or in examples directory below):
cc -O -Aa -Ae +Z +DA1.0 -DHPUX  -I /usr/local/include -I.. -L.. -o cover cover.c -lcpdfm -lm

1998-10-18 [IO]
*/

/* Enable this line if you want RGB Color In-line image at the bottom of page 1. 
   Commenting it out produces gray scale sine wave for it */
#define COLOR_INLINE_IMAGE 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <cpdflib.h>

#define WIDTH 		6.0
#define HEIGHT		4.0
#define XBASE		1.5
#define YBASE		4.0
#define NCYCLES		10.0
#define WEBSITE_STR	"http://www.fastio.com/"
#define EMAIL_STR	"clibpdf@fastio.com"
#define PAGE2TEXT	"page2.txt"

float random01(void);
void do_page_2(CPDFdoc *pdf);


float random01(void)
{
float v;
    v = (float)rand();
    v /= (float)RAND_MAX;
    return(v);
}

int main(int argc, char *argv[])
{
CPDFdoc *pdf;
int i;
float radius;
float xorig, yorig;
float xr, yr;
float sangle, eangle;
float r, g, b;
float strwidth1, strwidth2;
unsigned char sinewave[1024*3];
#ifdef COLOR_INLINE_IMAGE
float rfreq = 10.0;
float gfreq = 6.0;
float bfreq = 16.0;
#else
float angle;
#endif

    /* == Initialization == */
    srand(23456);
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);

    /* PDF Info object */
    cpdf_setCreator(pdf, "FastIO Systems: cover.c");
    cpdf_setTitle(pdf, "ClibPDF Reference Manual");
    cpdf_setSubject(pdf, "ANSI C Library for Direct PDF Generation");
    cpdf_setKeywords(pdf, "ClibPDF, ANSI C Library, Acrobat, PDF, Dynamic Web, Graph, Plot");

    cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */

    /* == Simple graphics drawing example == */
    cpdf_setgrayStroke(pdf, 0.0);	/* black as stroke color */
    cpdf_setlinewidth(pdf, 0.2);

    /* Install clipping path for the rectangle */
    cpdf_gsave(pdf);			/* needed to remove clip path later */
    cpdf_newpath(pdf);
    cpdf_rect(pdf, XBASE, YBASE, WIDTH, HEIGHT);
    cpdf_clip(pdf);
    cpdf_newpath(pdf);			/* clip does not consume path */

    /* Just a bunch of random pie/pacmans */
    cpdf_comments(pdf, "\n%% random pie/pacmans.\n");
    for(i=0; i<160; i++) {
	radius = 0.2 + random01();
	eangle = 360.0*random01();	/* end angle */
	sangle = 360.0*random01();	/* start angle */
	xr = WIDTH * random01();	/* position within the box */
	yr = HEIGHT * random01();
	r = 0.7*random01() + 0.3;	/* random RGB but desaturate a bit */
	g = 0.7*random01() + 0.3;
	b = 0.7*random01() + 0.3;
	cpdf_setrgbcolorFill(pdf, r, g, b);
	xorig = XBASE + xr;		/* center of arc */
	yorig = YBASE + yr;
	cpdf_moveto(pdf, xorig, yorig);
	cpdf_arc(pdf, xorig, yorig, radius, sangle, eangle, 0);	/* lineto to starting point of arc */
	cpdf_closepath(pdf);
        cpdf_fillAndStroke(pdf);		/* fill and stroke */
    }
    cpdf_grestore(pdf);			/* remove clipping path */

    cpdf_rect(pdf, XBASE, YBASE, WIDTH, HEIGHT);	/* stroke the domain border of pie/pacman field */
    cpdf_setlinewidth(pdf, 0.8);
    cpdf_stroke(pdf);

    /* == Text examples == */
    cpdf_setrgbcolorFill(pdf, 0.0, 0.0, 1.0);				/* Blue */
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Times-BoldItalic", "MacRomanEncoding", 38.0);
    cpdf_textAligned(pdf, 4.5, 9.0, 0.0, TEXTPOS_LM, "ClibPDF Reference Manual");
    cpdf_endText(pdf);

    cpdf_setrgbcolorFill(pdf, 0.4, 0.3, 0.2);
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Times-Roman", "MacRomanEncoding", 20.0);
    cpdf_textAligned(pdf, 4.5, 8.4, 0.0, TEXTPOS_LM, "ANSI C Library for Direct PDF Generation");
    cpdf_endText(pdf);

    cpdf_setgray(pdf, 0.0);
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Helvetica-Bold", "MacRomanEncoding", 18.0);
    cpdf_setHorizontalScaling(pdf, 140.0);
    cpdf_setCharacterSpacing(pdf, 1.0);
    cpdf_textAligned(pdf, 4.5, 2.8, 0.0, TEXTPOS_LM, "FastIO Systems");
    cpdf_endText(pdf);

    cpdf_setgray(pdf, 0.0);
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Helvetica", "MacRomanEncoding", 14.0);
    cpdf_setHorizontalScaling(pdf, 120.0);
    cpdf_setCharacterSpacing(pdf, 0.5);
    cpdf_textAligned(pdf, 4.5, 2.4, 0.0, TEXTPOS_LM, WEBSITE_STR);
    cpdf_textAligned(pdf, 4.5, 2.1, 0.0, TEXTPOS_LM, EMAIL_STR);
    strwidth1 = cpdf_stringWidth(pdf, WEBSITE_STR);
    strwidth2 = cpdf_stringWidth(pdf, EMAIL_STR);
    cpdf_endText(pdf);

    /* Add Web and email links over the address text above. */
    cpdf_rawSetActionURL(pdf, 4.5*inch - strwidth1/2.0, 2.4*inch - 3.0,
			 4.5*inch + strwidth1/2.0, 2.4*inch + 12.0,
			 WEBSITE_STR, NULL);
    cpdf_rawSetActionURL(pdf, 4.5*inch - strwidth2/2.0, 2.1*inch - 3.0,
			 4.5*inch + strwidth2/2.0, 2.1*inch + 12.0,
			 "mailto:clibpdf@fastio.com", NULL);

#ifdef COLOR_INLINE_IMAGE
    /* colored sine wave density plot with 24-bit RGB (8 bits/component) */
    for(i=0; i<1024*3; i+=3) {
        float rangle, gangle, bangle;
	rangle = (float)i/512.0 * PI * 2.0 * rfreq;
	gangle = (float)i/512.0 * PI * 2.0 * gfreq;
	bangle = (float)i/512.0 * PI * 2.0 * bfreq;
	if(i < 512*3) {
	    sinewave[i] = (unsigned char)(127.0 + 60.0*sin(rangle));
	    sinewave[i+1] = (unsigned char)(127.0 + 60.0*sin(gangle));
	    sinewave[i+2] = (unsigned char)(127.0 + 60.0*sin(bangle));
	}
	else {
	    sinewave[i] = (unsigned char)(127.0 - 60.0*sin(rangle));
	    sinewave[i+1] = (unsigned char)(127.0 - 60.0*sin(gangle));
	    sinewave[i+2] = (unsigned char)(127.0 - 60.0*sin(bangle));
	}
    }

    cpdf_placeInLineImage(pdf, (void *)sinewave, 1024*3,
		XBASE, 1.0, 0.0, WIDTH*72.0, 0.2*72.0,
		512, 2, 8, CS_RGB, 1);

#else
    /* Gray scale sine wave density plot */
    for(i=0; i<1024; i++) {
	angle = (float)i/512.0 * PI * 2.0 * NCYCLES;
	if(i < 512)
	    sinewave[i] = (unsigned char)(127.0 + 60.0*sin(angle));
	else
	    sinewave[i] = (unsigned char)(127.0 - 60.0*sin(angle));
    }

    cpdf_placeInLineImage(pdf, (void *)sinewave, 1024,
		XBASE, 1.0, 0.0, WIDTH*72.0, 0.2*72.0,
		512, 2, 8, CS_GRAY, 1);
#endif

    cpdf_rect(pdf, XBASE, 1.0, WIDTH, 0.2);	/* put a box around it */
    cpdf_setgray(pdf, 0.0);
    cpdf_setlinewidth(pdf, 0.2);
    cpdf_stroke(pdf);

    /* === Do page 2, mostly typical fine print stuff === */
    do_page_2(pdf);

    /* === All drawing done.  Finalize the content === */
    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "cover.pdf");

    /* == Clean up == */
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return(0);
}


void do_page_2(CPDFdoc *pdf)
{
FILE *fp;
char linebuf[1024];
    cpdf_pageInit(pdf, 2, PORTRAIT, LETTER, LETTER);		/* page orientation */
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Helvetica", "MacRomanEncoding", 7.0);
    cpdf_setHorizontalScaling(pdf, 130.0);
    cpdf_setTextLeading(pdf, 11.0);
    cpdf_setTextPosition(pdf, 1.0, 10.0);
    cpdf_textShow(pdf, "ClibPDF Library Reference Manual");

    if((fp = fopen(PAGE2TEXT, "r")) != NULL) {
	while( ( fgets(linebuf, 1022, fp) ) != NULL )
	    cpdf_textCRLFshow(pdf, linebuf);
	fclose(fp);
    }
    cpdf_endText(pdf);
}

