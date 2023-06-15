/* textboxfit.c -- TextBox that will fit given text into a box by shrinking.
 * Copyright (C) 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

cc -Wall -o textboxfit -I/usr/local/include textboxfit.c -lcpdfpm

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cpdflib.h>

/* Make sure to end the text with at least one '\n'.  If you do justification, use 2 consecutive
   new lines "\n\n" to end text.  This is necessary for correct last line handling.
*/
char *textdata = "We the People of the United States, in Order to form a more perfect Union,\
 establish Justice, insure domestic Tranquility, provide for the common\
 defence, promote the general Welfare, and secure the Blessings of Liberty to\
 ourselves and our Posterity, do ordain and establish this Constitution for the\
 United States of America.\n\n";

int main(int argc, char *argv[])
{
CPDFdoc *pdf;
float fontsize = 14.0, fontsizeused;
float angle = 0.0;
float x = 1.0, y = 7.0;		/* bottom left of box */
float w = 3.0, h = 2.0;		/* size of box */
float xr, yr;
char textbuf[4096], strbuf[128];
CPDFtboxAttr tbA;

    tbA.alignmode = TBOX_JUSTIFY;
    tbA.NLmode = 0;		/* zero reformat, nonzero = honor NL */
    tbA.paragraphSpacing = 0.0;
    tbA.noMark = 0;		/* if nonzero, do not actually paint text (invisible) */

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);
    cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */

    /* Upright text box */
    strncpy(textbuf, textdata, 4095);	/* cpdf_textBox() modifies text, so use a copy */
    /* show the box */
    cpdf_setrgbcolor(pdf, 0.5, 0.5, 1.0);
    cpdf_rect(pdf, x, y, w, h);
    cpdf_stroke(pdf);
    cpdf_setgray(pdf, 0.0);
    /* do text box at (x, y) */
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Helvetica", "MacRomanEncoding", fontsize);
    fontsizeused = cpdf_textBoxFit(pdf, x, y, w, h, angle, fontsize, 0.5, -1.2, &tbA, textbuf);
    sprintf(strbuf, "Font size: asked=%g,  used=%g", fontsize, fontsizeused);
    cpdf_setrgbcolor(pdf, 1.0, 0.0, 0.0);
    cpdf_textAligned(pdf, x, y-0.03, angle, TEXTPOS_UL, strbuf);
    cpdf_endText(pdf);
    
    /* Rotated text box with the same locaton for the bottom-left corner of the box as above */
    strncpy(textbuf, textdata, 4095);
    /* show the box */
    angle = 25.0;
    x = 5.0; y = 7.0;		/* bottom left of box */
    w = 1.5; h = 1.5;		/* size of box */
    cpdf_setrgbcolor(pdf, 0.5, 0.5, 1.0);
    cpdf_rectRotated(pdf, x, y, w, h, angle);
    cpdf_moveto(pdf, x, y);
    cpdf_lineto(pdf, x + w, y);
    cpdf_arc(pdf, x, y, 0.35, 0.0, angle, 1);
    cpdf_stroke(pdf);
    cpdf_setgray(pdf, 0.0);

    /* do text box at the new rotated origin */
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Helvetica", "MacRomanEncoding", fontsize);
    fontsizeused = cpdf_textBoxFit(pdf, x, y, w, h, angle, fontsize, 0.5, -1.2, &tbA, textbuf);

    /* legends */
    sprintf(strbuf, "Font size: asked=%g,  used=%g", fontsize, fontsizeused);
    cpdf_setrgbcolor(pdf, 1.0, 0.0, 0.0);
    cpdf_textAligned(pdf, x, y-0.05, 0.0, TEXTPOS_UL, strbuf);
    cpdf_setrgbcolor(pdf, 0.0, 0.0, 1.0);
    cpdf_setFont(pdf, "Helvetica-BoldOblique", "MacRomanEncoding", 8.0);
    cpdf_text(pdf, x + 0.37, y+ 0.05, 0.0, "angle");
    cpdf_textAligned(pdf, x-0.05, y, 0.0, TEXTPOS_MR, "(xl, yl)");
    xr = 0.75; yr = -0.025;
    rotate_xyCoordinate(xr, yr, angle, &xr, &yr);
    cpdf_textAligned(pdf, x+xr, y+yr, angle, TEXTPOS_UM, "width");
    xr = 0.75; yr = 0.05;
    rotate_xyCoordinate(xr, yr, angle+90.0, &xr, &yr);
    cpdf_textAligned(pdf, x+xr, y+yr, angle+90.0, TEXTPOS_LM, "height");
    cpdf_endText(pdf);

    cpdf_finalizeAll(pdf);		/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "textboxfit.pdf");

    /* == Clean up == */
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return 0;
}
