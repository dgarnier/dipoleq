/* cpdfDomPlot.c  -- plotting functions in raw domain of points
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

void cpdf_moveto(CPDFdoc *pdf, float x, float y)
{
	cpdf_rawMoveto(pdf, x_Domain2Points(pdf, x), y_Domain2Points(pdf, y));
}

void cpdf_lineto(CPDFdoc *pdf, float x, float y)
{
	cpdf_rawLineto(pdf, x_Domain2Points(pdf, x), y_Domain2Points(pdf, y));
}

void cpdf_rmoveto(CPDFdoc *pdf, float x, float y)
{
	cpdf_rawRmoveto(pdf, x_Domain2Points(pdf, x), y_Domain2Points(pdf, y));
}

void cpdf_rlineto(CPDFdoc *pdf, float x, float y)
{
	cpdf_rawRlineto(pdf, x_Domain2Points(pdf, x), y_Domain2Points(pdf, y));
}

void cpdf_curveto(CPDFdoc *pdf, float x1, float y1, float x2, float y2, float x3, float y3)
{
	cpdf_rawCurveto(pdf, x_Domain2Points(pdf, x1), y_Domain2Points(pdf, y1),
			x_Domain2Points(pdf, x2), y_Domain2Points(pdf, y2),
			x_Domain2Points(pdf, x3), y_Domain2Points(pdf, y3));
}

void cpdf_quickCircle(CPDFdoc *pdf, float x, float y, float r)	/* center (x,y) and radius r */
{
float convr;
	convr = x_Domain2Points(pdf, x+r) - x_Domain2Points(pdf, x);
	cpdf_rawQuickCircle(pdf, x_Domain2Points(pdf, x), y_Domain2Points(pdf, y), convr);
}

void cpdf_rect(CPDFdoc *pdf, float x, float y, float w, float h)
{
float convw, convh;
	convw = x_Domain2Points(pdf, x+w) - x_Domain2Points(pdf, x);
	convh = y_Domain2Points(pdf, y+h) - y_Domain2Points(pdf, y);
	cpdf_rawRect(pdf, x_Domain2Points(pdf, x), y_Domain2Points(pdf, y), convw, convh);
}

void cpdf_rectRotated(CPDFdoc *pdf, float x, float y, float w, float h, float angle)
{
float convw, convh;
	convw = x_Domain2Points(pdf, x+w) - x_Domain2Points(pdf, x);
	convh = y_Domain2Points(pdf, y+h) - y_Domain2Points(pdf, y);
	cpdf_rawRectRotated(pdf, x_Domain2Points(pdf, x), y_Domain2Points(pdf, y), convw, convh, angle);
}

