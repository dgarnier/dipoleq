/* cpdfGstate.c  -- graphics state 
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf.
 *  or http://www.fastio.com/LICENSE.pdf

1999-08-17 [IO] ## MULTI_THREADING_OK ##

1998-07-09 [IO]
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


void cpdf_gsave(CPDFdoc *pdf)
{
    if(pdf->useContentMemStream)
	cpdf_memPuts("\nq\n", pdf->currentMemStream);
    else
	fprintf(pdf->fpcontent, "\nq\n");
}

void cpdf_grestore(CPDFdoc *pdf)
{
    if(pdf->useContentMemStream)
	cpdf_memPuts("Q\n", pdf->currentMemStream);
    else
	fprintf(pdf->fpcontent, "Q\n");
}

/* E.g., cpdf_setdash("[6 1 2 1] 0"); for 6 ON, 1 OFF, 2 ON, 1 OFF, ... */
void cpdf_setdash(CPDFdoc *pdf, const char *dashspec)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%s d\n", dashspec);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%s d\n", dashspec);
}

void cpdf_nodash(CPDFdoc *pdf)
{
	cpdf_setdash(pdf, "[] 0");
}

void cpdf_setstrokeadjust(CPDFdoc *pdf, int flag)
{
/* static char *truefalse[] = {"false", "true", ""}; */
    if(pdf->pdfLevelMaj >= 1 && pdf->pdfLevelMin >= 2) {
	/* This can't be placed directly inside a content stream. Must use Ext graphics object. */
/*
	if(pdf->useContentMemStream) {
	    sprintf(pdf->spbuf, "%s SA\n", truefalse[flag&1]);
	    cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
	}
	else
	    fprintf(pdf->fpcontent, "%s SA\n", truefalse[flag&1]);
*/
    }
}

void cpdf_concat(CPDFdoc *pdf, float a, float b, float c, float d, float e, float f)
{
    cpdf_rawConcat(pdf, a, b, c, d, x_Domain2Points(pdf, e), y_Domain2Points(pdf, f));
}


/* concatenate CTM (coordinate transformation matrix) */
void cpdf_rawConcat(CPDFdoc *pdf, float a, float b, float c, float d, float e, float f)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.4f %.4f %.4f %.4f %.4f %.4f cm\n", a, b, c, d, e, f);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.4f %.4f %.4f %.4f %.4f %.4f cm\n", a, b, c, d, e, f);
}

void cpdf_rotate(CPDFdoc *pdf, float angle)
{
float vsin, vcos;
    angle *= PI/180.0;
    vsin = sin(angle);
    vcos = cos(angle);
    cpdf_rawConcat(pdf, vcos, vsin, -vsin, vcos, 0.0, 0.0);
}

void cpdf_translate(CPDFdoc *pdf, float xt, float yt)
{
    cpdf_rawTranslate(pdf, x_Domain2Points(pdf, xt), y_Domain2Points(pdf, yt));
}

void cpdf_rawTranslate(CPDFdoc *pdf, float xt, float yt)
{
    cpdf_rawConcat(pdf, 1.0, 0.0, 0.0, 1.0, xt, yt);
}

void cpdf_scale(CPDFdoc *pdf, float sx, float sy)
{
    cpdf_rawConcat(pdf, sx, 0.0, 0.0, sy, 0.0, 0.0);
}

void cpdf_setflat(CPDFdoc *pdf, int flatness)		/* flatness = 0 .. 100 */
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%d i\n", flatness);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%d i\n", flatness);
}

void cpdf_setlinejoin(CPDFdoc *pdf, int linejoin)	/* linejoin = 0(miter), 1(round), 2(bevel) */
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%d j\n", linejoin);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%d j\n", linejoin);
}

void cpdf_setlinecap(CPDFdoc *pdf, int linecap)	/* linecap = 0(butt end), 1(round), 2(projecting square) */
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%d J\n", linecap);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%d J\n", linecap);
}

void cpdf_setmiterlimit(CPDFdoc *pdf, float miterlimit)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.4f M\n", miterlimit);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.4f M\n", miterlimit);
}

void cpdf_setlinewidth(CPDFdoc *pdf, float width)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.4f w\n", width);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.4f w\n", width);
}


