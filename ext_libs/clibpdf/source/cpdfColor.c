/* cpdfColor.c  -- color functions
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf.
 *  or http://www.fastio.com/LICENSE.pdf

## MULTI_THREADING_OK ##

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



/* ## MULTI_THREADING_OK ## */
void cpdf_setgray(CPDFdoc *pdf, float gray)		/* set both fill and stroke grays */
{
	cpdf_setgrayFill(pdf, gray);
	cpdf_setgrayStroke(pdf, gray);
}


/* ## MULTI_THREADING_OK ## */
void cpdf_setrgbcolor(CPDFdoc *pdf, float r, float g, float b)	/* set both fill and stroke colors */
{
	cpdf_setrgbcolorFill(pdf, r, g, b);
	cpdf_setrgbcolorStroke(pdf, r, g, b);
}


/* ## MULTI_THREADING_OK ## */
void cpdf_setcmykcolor(CPDFdoc *pdf, float c, float m, float y, float k)
{
	cpdf_setcmykcolorFill(pdf, c, m, y, k);
	cpdf_setcmykcolorStroke(pdf, c, m, y, k);
}


/* ## MULTI_THREADING_OK ## */
void cpdf_setgrayFill(CPDFdoc *pdf, float gray)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.4f g\n", gray);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.4f g\n", gray);
}



/* ## MULTI_THREADING_OK ## */
void cpdf_setgrayStroke(CPDFdoc *pdf, float gray)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.4f G\n", gray);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.4f G\n", gray);
}



/* ## MULTI_THREADING_OK ## */
void cpdf_setrgbcolorFill(CPDFdoc *pdf, float r, float g, float b)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.4f %.4f %.4f rg\n", r, g, b);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.4f %.4f %.4f rg\n", r, g, b);
}



/* ## MULTI_THREADING_OK ## */
void cpdf_setrgbcolorStroke(CPDFdoc *pdf, float r, float g, float b)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.4f %.4f %.4f RG\n", r, g, b);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.4f %.4f %.4f RG\n", r, g, b);
}


/* ## MULTI_THREADING_OK ## */
void cpdf_setcmykcolorFill(CPDFdoc *pdf, float c, float m, float y, float k)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.4f %.4f %.4f %.4f k\n", c, m, y, k);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.4f %.4f %.4f %.4f k\n", c, m, y, k);
}


/* ## MULTI_THREADING_OK ## */
void cpdf_setcmykcolorStroke(CPDFdoc *pdf, float c, float m, float y, float k)
{
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "%.4f %.4f %.4f %.4f K\n", c, m, y, k);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
	fprintf(pdf->fpcontent, "%.4f %.4f %.4f %.4f K\n", c, m, y, k);
}

