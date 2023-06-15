/* cpdfRawPlot.c  -- plotting functions in raw domain of points
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

	Note: Better circle and arc (raw and domain) plot functions are in cpdfArc.c
1999-08-17 [IO] 
	## MULTI_THREADING_OK ##
*/

#include "config.h"
#include "version.h"

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cpdflib.h"		/* This must be included before all other local include files */
#include "cpdf_private.h"
#include "cglobals.h"



void cpdf_rawMoveto(CPDFdoc *pdf, float x, float y)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.3f %.3f m\n", x, y);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.3f %.3f m\n", x, y);
    pdf->xcurr = x; pdf->ycurr = y;
}

void cpdf_rawRmoveto(CPDFdoc *pdf, float x, float y)
{
    pdf->xcurr += x; pdf->ycurr += y;
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.3f %.3f m\n", pdf->xcurr, pdf->ycurr);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.3f %.3f m\n", pdf->xcurr, pdf->ycurr);
}


void cpdf_rawLineto(CPDFdoc *pdf, float x, float y)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.3f %.3f l\n", x, y);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.3f %.3f l\n", x, y);
    pdf->xcurr = x; pdf->ycurr = y;
}

void cpdf_rawRlineto(CPDFdoc *pdf, float x, float y)
{
    pdf->xcurr += x; pdf->ycurr += y;
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.3f %.3f l\n", pdf->xcurr, pdf->ycurr);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.3f %.3f l\n", pdf->xcurr, pdf->ycurr);
}


void cpdf_rawCurveto(CPDFdoc *pdf, float x1, float y1, float x2, float y2, float x3, float y3)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.3f %.3f %.3f %.3f %.3f %.3f c\n", x1, y1, x2, y2, x3, y3);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.3f %.3f %.3f %.3f %.3f %.3f c\n", x1, y1, x2, y2, x3, y3);
    pdf->xcurr = x3; pdf->ycurr = y3;
}

/* This is used only for circle marker in cpdfMark.c */
void cpdf_rawQuickCircle(CPDFdoc *pdf, float xc, float yc, float r)
{
float b = 1.3333333;
	cpdf_rawMoveto(pdf, xc-r, yc);
	cpdf_rawCurveto(pdf, xc-r, yc+b*r, xc+r, yc+b*r, xc+r, yc);
	cpdf_rawCurveto(pdf, xc+r, yc-b*r, xc-r, yc-b*r, xc-r, yc);
	cpdf_closepath(pdf);
	pdf->xcurr = xc; pdf->ycurr = yc;
}

void cpdf_rawRect(CPDFdoc *pdf, float x, float y, float w, float h)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.3f %.3f %.3f %.3f re\n", x, y, w, h);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.3f %.3f %.3f %.3f re\n", x, y, w, h);
    pdf->xcurr = x; pdf->ycurr = y;
}

void cpdf_rawRectRotated(CPDFdoc *pdf, float x, float y, float w, float h, float angle)
{
float xo, yo;
    cpdf_rawMoveto(pdf, x, y);
    rotate_xyCoordinate(w, 0.0, angle, &xo, &yo);
    cpdf_rawLineto(pdf, x + xo, y + yo);
    rotate_xyCoordinate(w, h, angle, &xo, &yo);
    cpdf_rawLineto(pdf, x + xo, y + yo);
    rotate_xyCoordinate(0.0, h, angle, &xo, &yo);
    cpdf_rawLineto(pdf, x + xo, y + yo);
    cpdf_closepath(pdf);
}

void cpdf_closepath(CPDFdoc *pdf)
{
    if(pdf->useContentMemStream)
	cpdf_memPuts("h\n", pdf->currentMemStream);
    else
	fprintf(pdf->fpcontent, "h\n");
}

void cpdf_stroke(CPDFdoc *pdf)
{
    if(pdf->useContentMemStream)
	cpdf_memPuts("S\n", pdf->currentMemStream);
    else
	fprintf(pdf->fpcontent, "S\n");
}

void cpdf_fill(CPDFdoc *pdf)
{
    if(pdf->useContentMemStream)
	cpdf_memPuts("f\n", pdf->currentMemStream);
    else
	fprintf(pdf->fpcontent, "f\n");
}

void cpdf_eofill(CPDFdoc *pdf)
{
    if(pdf->useContentMemStream)
	cpdf_memPuts("f*\n", pdf->currentMemStream);
    else
	fprintf(pdf->fpcontent, "f*\n");
}

void cpdf_fillAndStroke(CPDFdoc *pdf)
{
    if(pdf->useContentMemStream)
	cpdf_memPuts("B\n", pdf->currentMemStream);
    else
	fprintf(pdf->fpcontent, "B\n");
}

void cpdf_eofillAndStroke(CPDFdoc *pdf)
{
    if(pdf->useContentMemStream)
	cpdf_memPuts("B*\n", pdf->currentMemStream);
    else
	fprintf(pdf->fpcontent, "B*\n");
}

void cpdf_clip(CPDFdoc *pdf)
{
    if(pdf->useContentMemStream)
	cpdf_memPuts("W\n", pdf->currentMemStream);
    else
	fprintf(pdf->fpcontent, "W\n");
}

void cpdf_eoclip(CPDFdoc *pdf)
{
    if(pdf->useContentMemStream)
	cpdf_memPuts("W*\n", pdf->currentMemStream);
    else
	fprintf(pdf->fpcontent, "W*\n");
}

void cpdf_newpath(CPDFdoc *pdf)
{
    if(pdf->useContentMemStream)
	cpdf_memPuts("n\n", pdf->currentMemStream);
    else
	fprintf(pdf->fpcontent, "n\n");
}


