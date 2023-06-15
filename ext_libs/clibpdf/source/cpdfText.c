/* cpdfText.c
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

1999-08-17 [IO] 
	## MULTI_THREADING_OK ##

*/

/* #define CPDF_DEBUG 1 */

#include "config.h"
#include "version.h"

#include <string.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "cpdflib.h"		/* This must be included before all other local include files */
#include "cpdf_private.h"
#include "cglobals.h"

/* Current Text Matrix.  give Text Matrix concatenation behavior as Tm operator
   resets it fresh every time */

/* 2 x 3 matrix for CTM : matrix mult is done by adding 0 0 1 as a column (see cpdflib.h)
typedef struct {
    float a; float b;
    float c; float d;
    float x; float y;
} CPDFctm;
*/


/* ======================================================================= */
void _cpdf_resetTextCTM(CPDFdoc *pdf)
{
    pdf->textCTM.a = 1.0;
    pdf->textCTM.b = 0.0;
    pdf->textCTM.c = 0.0;
    pdf->textCTM.d = 1.0;
    pdf->textCTM.x = 0.0;
    pdf->textCTM.y = 0.0;
}


void cpdf_beginText(CPDFdoc *pdf, int clipmode)
{
int count, bufsize;
char *mbuff;
    cpdf_clearMemoryStream(pdf->scratchMem);
    cpdf_memPutc('\n', pdf->scratchMem);
    if(pdf->textClipMode) {
	cpdf_memPuts("q\n", pdf->scratchMem);
    }
    cpdf_memPuts("BT\n", pdf->scratchMem);
    cpdf_getMemoryBuffer(pdf->scratchMem, &mbuff, &count, &bufsize);

    if(pdf->useContentMemStream)
	cpdf_writeMemoryStream(pdf->currentMemStream, mbuff, count);
    else
	fputs(mbuff, pdf->fpcontent);

    pdf->inTextObj = 1;
    _cpdf_resetTextCTM(pdf);	/* FIXME: This probably should not be here. */
}

void cpdf_endText(CPDFdoc *pdf)
{
int count, bufsize;
char *mbuff;
    cpdf_clearMemoryStream(pdf->scratchMem);
    cpdf_memPuts("ET\n", pdf->scratchMem);
    if(pdf->textClipMode)
	cpdf_memPuts("Q\n", pdf->scratchMem);
    cpdf_getMemoryBuffer(pdf->scratchMem, &mbuff, &count, &bufsize);
    if(pdf->useContentMemStream)
	cpdf_writeMemoryStream(pdf->currentMemStream, mbuff, count);
    else
	fputs(mbuff, pdf->fpcontent);

    pdf->inTextObj = 0;
    pdf->textClipMode = 0;
}

/* Position and draw text in current domain using its coordinate system.
   This function will set a new Text CTM.
   Lower left corner of text string is place at (x, y).
*/
void cpdf_text(CPDFdoc *pdf, float x, float y, float orientation, const char *textstr)
{
float xc, yc;
    xc = x_Domain2Points(pdf, x);
    yc = y_Domain2Points(pdf, y);
    cpdf_rawText(pdf, xc, yc, orientation, textstr);
}

/* Position  and draw text in raw domain in points with rotation.
   This function will set a new Text CTM.
   Lower left corner of text string is place at (x, y).
 */
void cpdf_rawText(CPDFdoc *pdf, float x, float y, float orientation, const char *textstr)
{
float a, b, c, d, angle, vcos, vsin;

    angle = PI*orientation/180.0;
    vcos = cos(angle);
    vsin = sin(angle);
    a =  vcos;
    b =  vsin;
    c = -vsin;
    d =  vcos;
    cpdf_setTextMatrix(pdf, a, b, c, d, x, y);
    cpdf_textShow(pdf, textstr);
}


/* Position  and draw text in raw domain in points with rotation.
   This function will set a new Text CTM.
   Text string is placed at (x, y) with various justification mode.

#define	TEXTPOS_LL	0
#define	TEXTPOS_LM	1
#define	TEXTPOS_LR	2
#define	TEXTPOS_ML	3
#define	TEXTPOS_MM	4
#define	TEXTPOS_MR	5
#define	TEXTPOS_UL	6
#define	TEXTPOS_UM	7
#define	TEXTPOS_UR	8
 */

void cpdf_textAligned(CPDFdoc *pdf, float x, float y, float orientation, int centmode, const char *textstr)
{
    cpdf_rawTextAligned(pdf, x_Domain2Points(pdf, x), y_Domain2Points(pdf, y), orientation, centmode, textstr);
}

void cpdf_rawTextAligned(CPDFdoc *pdf, float x, float y, float orientation, int centmode, const char *textstr)
{
float a, b, c, d, angle, vcos, vsin;
float wd=0.0, hd=0.0, swidth, sheight;
float xc, yc;
int xcent, ycent;
    angle = PI*orientation/180.0;
    vcos = cos(angle);
    vsin = sin(angle);

    xcent = centmode % 3;	/* 0=left, 1=middle, 2=right */
    ycent = centmode / 3;	/* 0=lower,  1=middle, 2=upper */
    swidth = cpdf_stringWidth(pdf, (unsigned char *)textstr);	/* string width according to current font */
    /* sheight = FONTSIZE2HEIGHT * pdf->font_size; */		/* current font size */
    sheight = cpdf_capHeight(pdf);
    /* necessary offset along the string's width dimension */
    if(xcent == 2)
	wd = swidth;
    else if(xcent == 1)
	wd = 0.5 * swidth;
    else
	wd = 0.0;
    /* necessary offset along the string's height dimension */
    if(ycent == 2)
	hd = sheight;
    else if(ycent == 1)
	hd = 0.5 * sheight;
    else
	hd = 0.0;

    xc = x - (wd*vcos - hd*vsin);	/* correction to the X position in unroted coord */
    yc = y - (wd*vsin + hd*vcos);
    a =  vcos;
    b =  vsin;
    c = -vsin;
    d =  vcos;
    cpdf_setTextMatrix(pdf, a, b, c, d, xc, yc);
    cpdf_textShow(pdf, textstr);
}

void cpdf_setWordSpacing(CPDFdoc *pdf, float spacing)
{
    pdf->word_spacing = spacing;
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.3f Tw\n", spacing);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.3f Tw\n", spacing);
}


void cpdf_setCharacterSpacing(CPDFdoc *pdf, float spacing)
{
    pdf->char_spacing = spacing;
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.3f Tc\n", spacing);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.3f Tc\n", spacing);
}

void cpdf_setHorizontalScaling(CPDFdoc *pdf, float scale)
{
    pdf->horiz_scaling = scale;
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.2f Tz\n", scale);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.2f Tz\n", scale);
}

void cpdf_setTextLeading(CPDFdoc *pdf, float leading)
{
    pdf->text_leading = leading;
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.2f TL\n", leading);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.2f TL\n", leading);
}

/* Only in Text mode */
void cpdf_setTextMatrix(CPDFdoc *pdf, float a, float b, float c, float d, float x, float y)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.4f %.4f %.4f %.4f %.4f %.4f Tm\n", a, b, c, d, x, y);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.4f %.4f %.4f %.4f %.4f %.4f Tm\n", a, b, c, d, x, y);

     pdf->textCTM.a = a;
     pdf->textCTM.b = b;
     pdf->textCTM.c = c;
     pdf->textCTM.d = d;
     pdf->textCTM.x = x;
     pdf->textCTM.y = y;
}

/* We keep a concatenated Text CTM in C, update it and set it using "Tm" operator.
*/
void cpdf_concatTextMatrix(CPDFdoc *pdf, float a, float b, float c, float d, float x, float y)
{
CPDFctm  sCTM;
    sCTM.a = a;
    sCTM.b = b;
    sCTM.c = c;
    sCTM.d = d;
    sCTM.x = x;
    sCTM.y = y;
    multiplyCTM(&(pdf->textCTM), &sCTM);
    cpdf_setTextMatrix(pdf, pdf->textCTM.a, pdf->textCTM.b, pdf->textCTM.c,
		pdf->textCTM.d, pdf->textCTM.x, pdf->textCTM.y);
}

/* Only in Text mode */
void cpdf_rotateText(CPDFdoc *pdf, float degrees)
{
float a, b, c, d, e, f, angle, vcos, vsin;

    angle = PI*degrees/180.0;
    vcos = cos(angle);
    vsin = sin(angle);
    a =  vcos;
    b =  vsin;
    c = -vsin;
    d =  vcos;
    e = 0.0; f = 0.0;
    cpdf_concatTextMatrix(pdf, a, b, c, d, e, f);
}

/* Only in Text mode */
void cpdf_skewText(CPDFdoc *pdf, float alpha, float beta)
{
float a, b, c, d, e, f;
    a =  1.0;
    b =  tan(PI*alpha/180.0);
    c =  tan(PI*beta/180.0);
    d =  1.0;
    e = 0.0; f = 0.0;
    cpdf_concatTextMatrix(pdf, a, b, c, d, e, f);
}

/* Only in Text mode */
void cpdf_setTextPosition(CPDFdoc *pdf, float x, float y)
{
    cpdf_rawSetTextPosition(pdf, x_Domain2Points(pdf, x), y_Domain2Points(pdf, y));
}

/* Changed at version 1.10-beta5 to keep current skew, scaling, and rotation
   as per suggestion by Alex Valyayev. [1999-04-07]
*/
void cpdf_rawSetTextPosition(CPDFdoc *pdf, float x, float y)
{
/* float a, b, c, d, e, f; */
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.4f %.4f %.4f %.4f %.4f %.4f Tm\n",
			pdf->textCTM.a, pdf->textCTM.b, pdf->textCTM.c, pdf->textCTM.d, x, y);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.4f %.4f %.4f %.4f %.4f %.4f Tm\n",
			pdf->textCTM.a, pdf->textCTM.b, pdf->textCTM.c, pdf->textCTM.d, x, y);
    pdf->textCTM.x = x;
    pdf->textCTM.y = y;

    /* Hmm, should I use textCTM.a here???? */
/*
    a =  1.0;		
    b =  0.0;
    c =  0.0;
    d =  1.0;
    e = x; f = y;
    cpdf_setTextMatrix(pdf, a, b, c, d, e, f);
*/
}


/* In order to change Text Rendering Mode, you must first end text, and
   begin new text.
#define	TEXT_FILL		0
#define	TEXT_STROKE		1
#define	TEXT_FILL_STROKE	2
#define	TEXT_INVISIBLE		3
#define	TEXT_FILL_CLIP		4
#define	TEXT_STROKE_CLIP	5
#define	TEXT_FILL_STROKE_CLIP	6
#define	TEXT_CLIP		7

*/
void cpdf_setTextRenderingMode(CPDFdoc *pdf, int mode)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%d Tr\n", mode);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%d Tr\n", mode);
}

void cpdf_setTextRise(CPDFdoc *pdf, float rise)
{
    pdf->text_rise = rise;
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.3f Ts\n", rise);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.3f Ts\n", rise);
}

/* This sends Td operator: only in Text mode */
void cpdf_rawSetNextTextLineOffset(CPDFdoc *pdf, float x, float y)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.3f %.3f Td\n", x, y);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.3f %.3f Td\n", x, y);
}

/* Only in Text mode */
void cpdf_setNextTextLineOffset(CPDFdoc *pdf, float x, float y)
{
    cpdf_rawSetNextTextLineOffset(pdf, x_Domain2Points(pdf, x), y_Domain2Points(pdf, y));
}


/* This sends a T* operator: only in Text mode */
void cpdf_textCRLF(CPDFdoc *pdf)
{
    if(pdf->useContentMemStream)
	cpdf_memPuts("T*\n", pdf->currentMemStream);
    else
	fprintf(pdf->fpcontent, "T*\n");
}

/* Tj operator: only in Text mode */
void cpdf_textShow(CPDFdoc *pdf, const char *txtstr)
{
    if(pdf->hex_string) {	/* Hex encoded Unicode string */
	long length, lenEsc;
	char *escUnicode;	/* will point to escaped Unicode binary */
	char *tempbuf = (char *)malloc((size_t)(strlen(txtstr)/2 + 3));
	_cpdf_malloc_check((void *)tempbuf);
	cpdf_convertHexToBinary(txtstr, (unsigned char *)tempbuf, &length);
	escUnicode = cpdf_escapeSpecialCharsBinary(tempbuf, length, &lenEsc);	/* free(escUnicode) ! */
	free(tempbuf);
	
	if(pdf->useContentMemStream) {
	    cpdf_memPutc('(', pdf->currentMemStream);
	    cpdf_writeMemoryStream(pdf->currentMemStream, escUnicode, lenEsc);
	    cpdf_memPuts(") Tj\n", pdf->currentMemStream);
	}
	else {
	    fputc('(', pdf->fpcontent);
	    fwrite(escUnicode, 1, lenEsc, pdf->fpcontent);
	    fprintf(pdf->fpcontent, ") Tj\n");
	}
	free(escUnicode);
    }
    else {	/* Normal ASCII string, and not a hexed Unicode string */
        char *fixedstr = cpdf_escapeSpecialChars(txtstr);
	if(pdf->useContentMemStream) {
	    sprintf(pdf->spbuf, "(%s) Tj\n", fixedstr);
	    cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
	}
	else {    
	    fprintf(pdf->fpcontent, "(%s) Tj\n", fixedstr);
	}
	if(fixedstr)
	    free(fixedstr);
    }
}


/* ' operator: only in Text mode */
void cpdf_textCRLFshow(CPDFdoc *pdf, const char *txtstr)
{
    if(pdf->hex_string) {	/* Hex encoded Unicode string */
	long length, lenEsc;
	char *escUnicode;	/* will point to escaped Unicode binary */
	char *tempbuf = (char *)malloc((size_t)(strlen(txtstr)/2 + 3));
	_cpdf_malloc_check((void *)tempbuf);
	cpdf_convertHexToBinary(txtstr, (unsigned char *)tempbuf, &length);
	escUnicode = cpdf_escapeSpecialCharsBinary(tempbuf, length, &lenEsc);	/* free(escUnicode) ! */
	free(tempbuf);
	
	if(pdf->useContentMemStream) {
	    cpdf_memPutc('(', pdf->currentMemStream);
	    cpdf_writeMemoryStream(pdf->currentMemStream, escUnicode, lenEsc);
	    cpdf_memPuts(") '\n", pdf->currentMemStream);
	}
	else {
	    fputc('(', pdf->fpcontent);
	    fwrite(escUnicode, 1, lenEsc, pdf->fpcontent);
	    fprintf(pdf->fpcontent, ") '\n");
	}
	free(escUnicode);
    }
    else {	/* Normal ASCII string, and not a hexed Unicode string */
	char *fixedstr = cpdf_escapeSpecialChars(txtstr);
	if(pdf->useContentMemStream) {
	    sprintf(pdf->spbuf, "(%s) '\n", fixedstr);
	    cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
	}
	else
	    fprintf(pdf->fpcontent, "(%s) '\n", fixedstr);
	free(fixedstr);		/* MUST */
    }
}

/* This function will return a new string (in malloced memory)
   that contains appropriately escaped version of the input string.
   It will escape '\', '(', and ')'.
   ### IMPORTANT: You must FREE the returned string after use. ###
   Otherwise, memory leak will result.
*/

char *cpdf_escapeSpecialChars(const char *instr)
{
char *ptr, *ptr2, *buf, ch;
int escapecount = 0;
	ptr = (char *)instr;
	while( (ch = *ptr++) != '\0') {
	    if( (ch == '(') || (ch == ')') || (ch == '\\') || (ch == '\r'))
		escapecount++;
	}

	ptr = (char *)instr;
	/* Enough string length for new escaped string. */
	buf = (char *)malloc((size_t)(strlen(instr) + escapecount + 1));
	_cpdf_malloc_check((void *)buf);
	ptr2 = buf;
	while( (ch = *ptr++) != '\0') {
	    if( (ch == '\\') || (ch == '(') || (ch == ')') ) {
		*ptr2++ = '\\';
		*ptr2++ = ch;
	    }
	    else if(ch == '\r') {
		*ptr2++ = '\\';
		*ptr2++ = 'r';
	    }
	    else
		*ptr2++ = ch;
	}
	*ptr2 = '\0';	/* properly terminate the new string */
	return(buf);
}

/* This is a binary version of PDF special char escape function.
   For use with Unicode and Encrypted (may not be needed for encryption?) strings.
   This function will return a new binary data (in newly malloced memory)
   that contains appropriately escaped version of the input binary data.
   It will escape '\', '(', and ')'.
   ### IMPORTANT: You must FREE the returned buffer after use. ###
   Otherwise, memory leak will result.
*/

char *cpdf_escapeSpecialCharsBinary(const char *instr, long lengthIn, long *lengthOut)
{
char *ptr, *ptr2, *buf, ch;
long i;
int escapecount = 0;
	ptr = (char *)instr;
	for(i=0; i<lengthIn; i++) {
	    ch = *ptr++;
	    if( (ch == '(') || (ch == ')') || (ch == '\\') || (ch == '\r') )
		escapecount++;
	}
	*lengthOut = lengthIn + escapecount;

	ptr = (char *)instr;
	/* Enough string length for new escaped string. */
	buf = (char *)malloc((size_t)(lengthIn + escapecount + 1));
	_cpdf_malloc_check((void *)buf);
	ptr2 = buf;
	for(i=0; i<lengthIn; i++) {
	    ch = *ptr++;
	    if( (ch == '\\') || (ch == '(') || (ch == ')')  ) {
		*ptr2++ = '\\';
		*ptr2++ = ch;
	    }
	    else if(ch == '\r') {
		*ptr2++ = '\\';
		*ptr2++ = 'r';
	    }
	    else
		*ptr2++ = ch;
	}
	*ptr2 = '\0';	/* not really needed */
	return(buf);
}


