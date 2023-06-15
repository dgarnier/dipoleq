/* cpdfTextBox.c  - Columnar text formatting function for ClibPDF with column chaining capability.
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

1999-11-18 [io]
	TextBox now honors '\f'  Form Feed character for text2pdf support for legacy
	ASCII files.
1999-08-17 [IO] 
	## MULTI_THREADING_OK ##

1999-06-06 [io] FastIO Systems
--
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

/* Now in cpdflib.h: For passing as "align" to cpdf_textBox() and cpdf_rawTextBox() */
/*
#define	TBOX_LEFT		0
#define	TBOX_CENTER		1
#define	TBOX_RIGHT		2
#define	TBOX_JUSTIFY		3
*/

/*
  This TextBox function returns the bottom Y position of the last line taking into account
  the Descent value of the current font.  The Y position is with respect to the base line
  of the TextBox, not with respect to the coordinate system in which TextBox sits.
  Before using this function, at a minimum, cpdf_beginText(), cpdf_setFont() must be called.
  After the call, cpdf_endText() must be called at some point. 
  ##FIXME (for the future): This won't work with CJK languages where there may be
  no space characters in text.
*/

char *cpdf_rawTextBoxY(CPDFdoc *pdf, float xl, float yl, float width, float height, float angle,
			float linespace, float *yend, CPDFtboxAttr *tbattr, char *text)
{
int  nspace=0;
int  alignLocal, isNL = 0, isNL2 = 0;
int  align = TBOX_LEFT;
int  NLisNL = 0;		/* if non-zero, NL is is a line break, if 0 reformatted */
int  nomark = 0;		/* if non-zero, will not actually draw text (for fitting text into box) */
int  isPageBreak = 0;
float x, y, wslac=0.0, swid=0.0, swidlast=0.0;
float xR, yR;			/* rotated x y position */
float paraSpace = 0.0;
float lineSpaceLocal = linespace;
char *currline, *p, *plast, *retptr;

	/* If CPDFtboxAttr *tbattr is NON-NULL, use the values */
	if(tbattr) {
	    align = tbattr->alignmode;
	    NLisNL = tbattr->NLmode;	/* if non-zero, NL is is a line break, if 0 reformatted */
	    paraSpace = tbattr->paragraphSpacing;
	    nomark = tbattr->noMark;
	}
	/* If linespace is negative, its absolute value is taken as the multiplying factor
	   for font size where the product is used as the line spaceing.
	   This way, line spaceing scales with font size.
	*/
	if(linespace < 0.0)
	    lineSpaceLocal = -linespace * pdf->font_size;
	/* same for paraSpace */
	if(paraSpace < 0.0)
	    paraSpace = -paraSpace * pdf->font_size;

	if(NLisNL && align == TBOX_JUSTIFY)
		align = TBOX_LEFT;		/* if NLisNL is ON, do not allow justification */

	retptr = NULL;
	y = height - cpdf_ascent(pdf);	/* first line's location */

	/* Adjust string length so that it will fit on the page */
	currline = text;	/* initialize currline pointer to the head of text */

	/* Loop while there is more text that fits. */
	while( *currline != '\0') {	/* 1-st while */
	    nspace = 0;		/* number of space chars */
	    p = plast = currline;
	    alignLocal = align;
	    if(align == TBOX_JUSTIFY)
		pdf->word_spacing = 0.0;		/* ClibPDF internal data */

	    /* look for next line break */
	    while( *p != '\0') {	/* 2-nd while */
		if( *p != ' ' && *p != '\n' && *p != '\f') {
		    p++;
		    continue;	/* Regular non-space char move to next char */
		}
		/* 'p' is now at ' ' or '\n' or '\f' */
		/* Treat LF LF as a line break or if NLisNL is true, treat single NL as a NL */
		if(((*p == '\n' || *p == '\f') && NLisNL) || (!NLisNL && *p == '\n' && *(p+1) == '\n')) {
		    char chsave = *p;
		    if(*p == '\f') {
			isPageBreak = 1;		/* form feed char */
			plast = p;
			break;
		    }
		    if(!NLisNL && *p == '\n' && *(p+1) == '\n') {
			while( *(p+1) == '\n') {
			    *p = ' ';
			    p++;
			}
		    }
		    *p = '\0';
		    isNL = 1;
		    if((swid = cpdf_stringWidth(pdf, (unsigned char *)currline)) <= width) {
			/* Still fits */
			plast = p;
			swidlast = swid;
			if(align == TBOX_JUSTIFY) {
		            alignLocal = TBOX_LEFT;	/* don't justify the last line of para. */
			    if( !nomark )
			        cpdf_setWordSpacing(pdf, 0.0);	/* this writes to stream */
			    else
				pdf->word_spacing = 0.0;
			}
			isNL = 1;
			break;		/* BREAK #1: we are within the right margin for NL */
		    }
		    else if(!NLisNL) {
			*(p-1) = '\n';
		    }
		    *p = chsave;
		}

		if(*p == '\n') isNL2 = 1;
		else           isNL2 = 0;
		*p = '\0';	/* break string at this space just for string width test */
		if((swid = cpdf_stringWidth(pdf, (unsigned char *)currline)) > width) {
		    if(isNL2) *p = '\n';
		    else      *p = ' ';
		    isNL = 0;
		    break;		/* BREAK #2 : we are over the right margin, so use previous break. */
		}
		/* still fits, so find next space */
		plast = p;	/* save the last word break at which the string still fits */
		swidlast = swid;
		nspace++;
		*p = ' ';
		p++;	/* move over to next char */
	    } /* end 2-nd while */

	    *plast = '\0';	/* use the last break */
	    switch(alignLocal) {
		default:
		case TBOX_LEFT:
			x = 0.0;
			break;
		case TBOX_CENTER:
			x =  0.5*(width - swidlast);
			break;
		case TBOX_RIGHT:
			x =  (width - swidlast);
			break;
		case TBOX_JUSTIFY:
			x = 0.0;
			if(nspace > 1) {
			    wslac = (width - swidlast)/(float)(nspace-1);
			    if( !nomark )
			        cpdf_setWordSpacing(pdf, wslac);	/* this writes to stream */
			    else
				pdf->word_spacing = wslac;
			}
			break;
	    }
	    /* ------------ This is where one text line is painted ----------- */
	    if( !nomark ) {
		rotate_xyCoordinate(x, y, angle, &xR, &yR);
	        cpdf_rawText(pdf, xR + xl, yR + yl, angle, currline);
	    }
	    *yend = y;			/* return the last Y position to caller [1999-10-15] */
	    y -= lineSpaceLocal;
	    if(isNL) {
		y -= paraSpace;		/* extra space for paragraph */
		*plast = '\n';
	    }
	    else *plast = ' ';
	    currline = ++plast;	/* start point of next line */
	    nspace = 0;		/* reset space count */

	    /* check for the bottom of the text box considering space for descenders */
	    if( ((y + cpdf_descent(pdf)) < 0.0 || isPageBreak) && *currline != '\0') {
		/* text did not fit, so give the pointer to the new currline as the return value */
		retptr = currline;
		break;
	    }
	}  /* end 1-st while */

	/* undo the effect of the line feed to return the last Y value */
	*yend += cpdf_descent(pdf);	/* Descent is negative */
			/* If text did not fit in the box, pointer to the remainder is returned. */
	return(retptr);	/* NULL return value indicates all text contained in the box */
}

/* The original cpdf_rawTextBox() -- ignores returned yend value */
char *cpdf_rawTextBox(CPDFdoc *pdf, float xl, float yl, float width, float height, float angle,
			float linespace, CPDFtboxAttr *tbattr, char *text)
{
float yend;
    return( cpdf_rawTextBoxY(pdf, xl, yl, width, height, angle, linespace, &yend, tbattr, text) );
}


float cpdf_rawTextBoxFit(CPDFdoc *pdf, float xl, float yl, float width, float height, float angle,
		float inifontsize, float fsdecrement, float linespace, CPDFtboxAttr *tbattr, char *text)
{
int noMarkSave = 0;
char *txbuf;
float ydummy, fsize = inifontsize;
CPDFtboxAttr tbA;

	if(tbattr) {
	    noMarkSave = tbattr->noMark;		/* save state for later */
	    tbA.alignmode = tbattr->alignmode;
	    tbA.NLmode = tbattr->NLmode;		/* if non-zero, NL is is a line break, if 0 reformatted */
	    tbA.paragraphSpacing = tbattr->paragraphSpacing;
	}
	else {
	    tbA.alignmode = TBOX_LEFT;
	    tbA.NLmode = 0;	/* if non-zero, NL is is a line break, if 0 reformatted */
	    tbA.paragraphSpacing = 0.0;
	}
	tbA.noMark = 1;	    

	txbuf = (char *)malloc((size_t)(strlen(text) + 2));
	_cpdf_malloc_check((void *)txbuf);

	for(fsize = inifontsize; ; fsize -= fsdecrement) {
	    strcpy(txbuf, text);	/* need fresh copy */
	    pdf->font_size = fsize;
	    if( cpdf_rawTextBoxY(pdf, xl, yl, width, height, angle, linespace, &ydummy, &tbA, txbuf) == NULL)
		break;		/* We have a fit ! */
	}

	/* Now we have font size that will fit the text into the box, so do it. */
	strcpy(txbuf, text);	/* need fresh copy */
	if(pdf->useContentMemStream) {
	    sprintf(pdf->spbuf, "/%s %.3f Tf\n", pdf->fontInfos[pdf->currentFont].name, fsize);
	    cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
	}
	else
	    fprintf(pdf->fpcontent, "/%s %.3f Tf\n", pdf->fontInfos[pdf->currentFont].name, fsize);
	tbA.noMark = noMarkSave;	/* restore marking flag */
	cpdf_rawTextBoxY(pdf, xl, yl, width, height, angle, linespace, &ydummy, &tbA, txbuf);
	free(txbuf);
	return(fsize);
}


/* scaled domain coordinate system version */
char *cpdf_textBoxY(CPDFdoc *pdf, float xl, float yl, float width, float height, float angle,
			float linespace, float *yend, CPDFtboxAttr *tbattr, char *text)
{
float y2;
char *retptr = cpdf_rawTextBoxY(pdf, x_Domain2Points(pdf, xl), y_Domain2Points(pdf, yl),
			x_Domain2Points(pdf, xl+width) - x_Domain2Points(pdf, xl),
			y_Domain2Points(pdf, yl+height) - y_Domain2Points(pdf, yl),
			angle, linespace, &y2, tbattr, text);
    *yend = y_Points2Domain(pdf, y2);	/* convert back from points to domain unit */
    return(retptr);
}

char *cpdf_textBox(CPDFdoc *pdf, float xl, float yl, float width, float height, float angle,
			float linespace, CPDFtboxAttr *tbattr, char *text)
{
    return( cpdf_rawTextBox(pdf, x_Domain2Points(pdf, xl), y_Domain2Points(pdf, yl),
			x_Domain2Points(pdf, xl+width) - x_Domain2Points(pdf, xl),
			y_Domain2Points(pdf, yl+height) - y_Domain2Points(pdf, yl),
			angle, linespace, tbattr, text) );
}

float cpdf_textBoxFit(CPDFdoc *pdf, float xl, float yl, float width, float height, float angle,
		float inifontsize, float fsdecrement, float linespace, CPDFtboxAttr *tbattr, char *text)
{
    return( cpdf_rawTextBoxFit(pdf, x_Domain2Points(pdf, xl), y_Domain2Points(pdf, yl),
			x_Domain2Points(pdf, xl+width) - x_Domain2Points(pdf, xl),
			y_Domain2Points(pdf, yl+height) - y_Domain2Points(pdf, yl),
			angle, inifontsize, fsdecrement, linespace, tbattr, text) );
}

