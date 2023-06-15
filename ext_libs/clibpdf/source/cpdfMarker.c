/* cpdfMarker.c  -- plot markers 
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf


1999-08-17 [IO] 
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



void cpdf_marker(CPDFdoc *pdf, float x, float y, int markertype, float size)
{
	cpdf_rawMarker(pdf, x_Domain2Points(pdf, x), y_Domain2Points(pdf, y), markertype, size);
}

void cpdf_pointer(CPDFdoc *pdf, float x, float y, int direction, float size)
{
	cpdf_rawPointer(pdf, x_Domain2Points(pdf, x), y_Domain2Points(pdf, y), direction, size);
}

/* draws errorbar with caps, capsize in points */
void cpdf_errorbar(CPDFdoc *pdf, float x, float y1, float y2, float capsize)
{
	cpdf_rawErrorbar(pdf, x_Domain2Points(pdf, x), y_Domain2Points(pdf, y1), y_Domain2Points(pdf, y2), capsize);
}

/* primarily for daily stock prices: high, low and close */
void cpdf_highLowClose(CPDFdoc *pdf, float x, float vhigh, float vlow, float vclose, float ticklen)
{
	cpdf_rawHighLowClose(pdf, x_Domain2Points(pdf, x),
	     y_Domain2Points(pdf, vhigh), y_Domain2Points(pdf, vlow), y_Domain2Points(pdf, vclose), ticklen);
}

/* =============================================================================== */
/* Raw (point-based corrdinate) data point markers and other symbols */

/* draw a marker at a raw point-based coordinate */
void cpdf_rawMarker(CPDFdoc *pdf, float x, float y, int markertype, float size)
{
float halfsize = size*0.5;

    cpdf_gsave(pdf);
    switch(markertype) {
	case 0:		/* circle */
    		cpdf_setlinewidth(pdf, size*0.15);
		cpdf_rawQuickCircle(pdf, x, y, halfsize);
    		cpdf_fillAndStroke(pdf);
		break;
	case 1:		/* up triangle */
    		cpdf_setlinewidth(pdf, size*0.15);
		cpdf_rawMoveto(pdf, x, y+size*0.6667);
		cpdf_rawRlineto(pdf, -size/1.7321, -size);
		cpdf_rawRlineto(pdf, 1.1546*size, 0.0);
		cpdf_closepath(pdf);
    		cpdf_fillAndStroke(pdf);
		break;
	case 2:		/* diamond */
    		cpdf_setlinewidth(pdf, size*0.15);
		size *= 0.9;
		cpdf_rawMoveto(pdf, x, y+size/1.38);
		cpdf_rawRlineto(pdf, 0.546*size, -size/1.38);
		cpdf_rawRlineto(pdf, -0.546*size, -size/1.38);
		cpdf_rawRlineto(pdf, -0.546*size, size/1.38);
		cpdf_closepath(pdf);
    		cpdf_fillAndStroke(pdf);
		break;
	case 3:		/* square */
    		cpdf_setlinewidth(pdf, size*0.15);
		cpdf_rawRect(pdf, x - halfsize, y - halfsize, size, size);
    		cpdf_fillAndStroke(pdf);
		break;
	case 4:		/* star */
		size *= 1.2;
		halfsize = 0.5*size;
    		cpdf_setlinewidth(pdf, size*0.09);
		cpdf_rawMoveto(pdf, x, y+size*0.5);
		cpdf_rawLineto(pdf, x+0.112255*size, y+0.15451*size);
		cpdf_rawLineto(pdf, x+0.47552*size, y+0.15451*size);
		cpdf_rawLineto(pdf, x+0.181635*size, y-0.05902*size);
		cpdf_rawLineto(pdf, x+0.29389*size, y-0.40451*size);
		cpdf_rawLineto(pdf, x, y-0.19098*size);
		cpdf_rawLineto(pdf, x-0.29389*size, y-0.40451*size);
		cpdf_rawLineto(pdf, x-0.181635*size, y-0.05902*size);
		cpdf_rawLineto(pdf, x-0.47552*size, y+0.15451*size);
		cpdf_rawLineto(pdf, x-0.112255*size, y+0.15451*size);
		cpdf_closepath(pdf);
		cpdf_fillAndStroke(pdf);
		break;
	case 5:		/* pentagon */
    		cpdf_setlinewidth(pdf, size*0.15);
		cpdf_rawMoveto(pdf, x-0.5257*size, y+size*0.1708);
		cpdf_rawRlineto(pdf, 0.5257*size, 0.382*size);
		cpdf_rawRlineto(pdf, 0.5257*size, -0.382*size);
		cpdf_rawRlineto(pdf, -0.2008*size, -0.6181*size);
		cpdf_rawRlineto(pdf, -0.6499*size, 0.0);
		cpdf_closepath(pdf);
		cpdf_fillAndStroke(pdf);
		break;
	case 6:		/* down triangle */
    		cpdf_setlinewidth(pdf, size*0.15);
		cpdf_rawMoveto(pdf, x, y-size*0.6667);
		cpdf_rawRlineto(pdf, -size/1.7321, size);
		cpdf_rawRlineto(pdf, 1.1546*size, 0.0);
		cpdf_closepath(pdf);
    		cpdf_fillAndStroke(pdf);
		break;
	case 7:		/* horiz bow-tie */
    		cpdf_setlinewidth(pdf, size*0.13);
		cpdf_rawMoveto(pdf, x-0.5*size, y-0.5*size);
		cpdf_rawLineto(pdf, x+0.5*size, y+0.5*size);
		cpdf_rawLineto(pdf, x+0.5*size, y-0.5*size);
		cpdf_rawLineto(pdf, x-0.5*size, y+0.5*size);
		cpdf_closepath(pdf);
    		cpdf_fillAndStroke(pdf);
		break;
	case 8:		/* vert bow-tie */
    		cpdf_setlinewidth(pdf, size*0.13);
		cpdf_rawMoveto(pdf, x-0.5*size, y-0.5*size);
		cpdf_rawLineto(pdf, x+0.5*size, y+0.5*size);
		cpdf_rawLineto(pdf, x-0.5*size, y+0.5*size);
		cpdf_rawLineto(pdf, x+0.5*size, y-0.5*size);
		cpdf_closepath(pdf);
    		cpdf_fillAndStroke(pdf);
		break;
		break;
	case 9:		/* asterisk */
		size *= 1.05;
    		cpdf_setlinewidth(pdf, size*0.15);
		cpdf_rawMoveto(pdf, x, y+size*0.5);
		cpdf_rawRlineto(pdf, 0.0, -size);
		cpdf_rawMoveto(pdf, x+0.433*size, y+0.25*size);
		cpdf_rawLineto(pdf, x-0.433*size, y-0.25*size);
		cpdf_rawMoveto(pdf, x+0.433*size, y-0.25*size);
		cpdf_rawLineto(pdf, x-0.433*size, y+0.25*size);
    		cpdf_stroke(pdf);
		break;
	case 10:	/* sun */
    		cpdf_setlinewidth(pdf, size*0.15);
		cpdf_rawQuickCircle(pdf, x, y, size*0.25);
		cpdf_rawMoveto(pdf, x+size*0.5, y);
		cpdf_rawLineto(pdf, x+size*0.25, y);
		cpdf_rawMoveto(pdf, x-size*0.5, y);
		cpdf_rawLineto(pdf, x-size*0.25, y);
		cpdf_rawMoveto(pdf, x, y-size*0.5);
		cpdf_rawLineto(pdf, x, y-size*0.25);
		cpdf_rawMoveto(pdf, x, y+size*0.5);
		cpdf_rawLineto(pdf, x, y+size*0.25);
    		cpdf_fillAndStroke(pdf);
		break;

	default:
		break;
    } /* end of switch() */
    cpdf_grestore(pdf);
}


/* draw a marker at a raw point-based coordinate */
void cpdf_rawPointer(CPDFdoc *pdf, float x, float y, int direction, float size)
{
    cpdf_gsave(pdf);
    switch(direction & 3) {
	case 0:		/* right */
    		cpdf_setlinewidth(pdf, size*0.14);
		cpdf_rawMoveto(pdf, x, y);
		cpdf_rawLineto(pdf, x-size, y+size*0.3333);
		cpdf_rawLineto(pdf, x-size, y-size*0.3333);
		cpdf_closepath(pdf);
    		cpdf_fillAndStroke(pdf);
		break;
	case 1:		/* down */
    		cpdf_setlinewidth(pdf, size*0.14);
		cpdf_rawMoveto(pdf, x, y);
		cpdf_rawLineto(pdf, x+size*0.3333, y+size);
		cpdf_rawLineto(pdf, x-size*0.3333, y+size);
		cpdf_closepath(pdf);
    		cpdf_fillAndStroke(pdf);
		break;
	case 2:		/* left */
    		cpdf_setlinewidth(pdf, size*0.14);
		cpdf_rawMoveto(pdf, x, y);
		cpdf_rawLineto(pdf, x+size, y+size*0.3333);
		cpdf_rawLineto(pdf, x+size, y-size*0.3333);
		cpdf_closepath(pdf);
    		cpdf_fillAndStroke(pdf);
		break;
	case 3:		/* up */
    		cpdf_setlinewidth(pdf, size*0.14);
		cpdf_rawMoveto(pdf, x, y);
		cpdf_rawLineto(pdf, x+size*0.3333, y-size);
		cpdf_rawLineto(pdf, x-size*0.3333, y-size);
		cpdf_closepath(pdf);
    		cpdf_fillAndStroke(pdf);
		break;
	default:
		break;
    } /* end of switch() */
    cpdf_grestore(pdf);
}

void cpdf_rawErrorbar(CPDFdoc *pdf, float x, float y1, float y2, float capsize)
{
	cpdf_rawMoveto(pdf, x, y1);
	cpdf_rawLineto(pdf, x, y2);
	cpdf_rawMoveto(pdf, x-0.5*capsize, y1);
	cpdf_rawLineto(pdf, x+0.5*capsize, y1);
	cpdf_rawMoveto(pdf, x-0.5*capsize, y2);
	cpdf_rawLineto(pdf, x+0.5*capsize, y2);
	cpdf_stroke(pdf);
}

void cpdf_rawHighLowClose(CPDFdoc *pdf, float x, float vhigh, float vlow, float vclose, float ticklen)
{
	cpdf_rawMoveto(pdf, x, vlow);
	cpdf_rawLineto(pdf, x, vhigh);
	cpdf_rawMoveto(pdf, x, vclose);
	cpdf_rawLineto(pdf, x+ticklen, vclose);
	cpdf_stroke(pdf);

}


