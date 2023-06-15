/* cpdfArc.c
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf.
 *  or http://www.fastio.com/LICENSE.pdf

## MULTI_THREADING_OK  ##

	PDF lacks the "arc" operator, so here is an implementation in C which
	uses arcs implemented as Bezier curves of less than 90 degrees
	to produce arbitrary-angle arcs.  We compute all the control points
	in C and use PDF "curveto" operator multiple times as needed.

	We cannot duplicate all the functionality of the PostScript arc operator.
	Specifically, PostScript arc operator performs a test on whether
	there is a current point and does things differently depending on the
	result.  Since PDF lacks the programming language construct such as "if()",
	we cannot define an equivalent operator at the PDF level.  We can maintain
	a state for the presence of current point within C, but that seems overkill.
	We just let the programmer decide wether there should be a current
	point and let her pass that as the last argument "moveto0" in the
	arc functions here.

--
1998-09-05 [IO]
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


/* Draws a circle CCW (counter clock-wise) in the current domain coordinate.
   Circle is drawn by using 4 Bezier curves of 90 degrees each.
   It does a cosepath.
*/

void cpdf_circle(CPDFdoc *pdf, float x, float y, float r)
{
float convr;
	convr = x_Domain2Points(pdf, x+r) - x_Domain2Points(pdf, x);
	cpdf_rawCircle(pdf, x_Domain2Points(pdf, x), y_Domain2Points(pdf, y), convr);
}


/* Arc implementation where (x, y, r) are in the unit of current domain, and
       start and end angles (sangle, eangle) are in degrees.
   If eangle < sangle, arc is drawn CW (clock wise).
   If arc spans an angle > 90, it is drawn as multiple Bezier curves, each
	of which spans no more than 90 degrees.
   If moveto0 is non-zero, it will perform a moveto to the starting point of arc,
   Otherwise (moveto0 == 0), it assumes there is current point, and performs a
	lineto to the starting point of arc, and then draws the arc.
*/

void cpdf_arc(CPDFdoc *pdf, float x, float y, float r, float sangle, float eangle, int moveto0)
{
float convr;
	convr = x_Domain2Points(pdf, x+r) - x_Domain2Points(pdf, x);
	cpdf_rawArc(pdf, x_Domain2Points(pdf, x), y_Domain2Points(pdf, y), convr, sangle, eangle, moveto0);
}


/* This is a 4-arc version more accurate than the 2-arc version in cpdf_rawQuickCircle() */
/* Draws a circle CCW (counter clock-wise) in (xc, yc, r) are in the unit of points.
   Circle is drawn by using 4 Bezier curves of 90 degrees each.
   It does a cosepath.
*/

void cpdf_rawCircle(CPDFdoc *pdf, float xc, float yc, float r)
{
float b = 0.552;
	cpdf_rawMoveto(pdf, xc+r, yc);
	cpdf_rawCurveto(pdf, xc+r, yc+b*r, xc+b*r, yc+r, xc, yc+r);  /* Q1 */
	cpdf_rawCurveto(pdf, xc-b*r, yc+r, xc-r, yc+b*r, xc-r, yc);  /* Q2 */
	cpdf_rawCurveto(pdf, xc-r, yc-b*r, xc-b*r, yc-r, xc, yc-r);  /* Q3 */
	cpdf_rawCurveto(pdf, xc+b*r, yc-r, xc+r, yc-b*r, xc+r, yc);  /* Q4 */
	cpdf_closepath(pdf);
}

/* Arc implementation where (x, y, r) are in unit of points, and
       start and end angles (sangle, eangle) are in degrees.
   If eangle < sangle, arc is drawn CW (clock wise).
   If arc spans an angle > 90, it is drawn as multiple Bezier curves, each
	of which spans no more than 90 degrees.
   If moveto0 is non-zero, it will perform a moveto to the starting point of arc,
   Otherwise (moveto0 == 0), it assumes there is current point, and performs a
	lineto to the starting point of arc, and then draws the arc.
*/

void cpdf_rawArc(CPDFdoc *pdf, float x, float y, float r, float sangle, float eangle, int moveto0)
{
int i;
int Narc = 1;
int ccwcw = 1;
float angleBump = 0.0, haBump = 0.0;
float aspan = eangle - sangle;
float currAngle= 0.0;
    if(eangle < sangle) ccwcw = -1;
    while( fabs(aspan)/(float)Narc > 90.0) Narc++;
    angleBump = aspan/(float)Narc;
    haBump = 0.5*angleBump;
    currAngle = sangle + haBump;
    for(i=0; i<Narc; i++) {
        /* fprintf(stderr, "ccwcw=%d,  currAngle=%g,   haBump=%g\n", ccwcw, currAngle, haBump); */
	if(i==0) {
	    if(moveto0)
	        _cpdf_arc_small(pdf, x, y, r, currAngle, haBump, 1, ccwcw);	/* do moveto */
	    else
	        _cpdf_arc_small(pdf, x, y, r, currAngle, haBump, -1, ccwcw);	/* do lineto from assumed CP */
	}
	else
	    _cpdf_arc_small(pdf, x, y, r, currAngle, haBump, 0, ccwcw);
	currAngle += angleBump;
    }
}


/*
   Small (less than 90 deg) arc symmetrical about mid-point angle midtheta
       with angular span of +/- htheta (i.e., total span of 2*htheta).
   This is a generalization to an arbitrary center angle of the method described
	in the Jim Fitzsimmons TeX tutorial.  We just use standard coordinate
	rotation to do that.  Small arcs implemented here are used multiple times
	above to make arcs spanning more than 90 degrees.
   Set ccwcw for CCW = 1, CW = -1.
   Set mvlnto0 = 1 to do initial moveto, mvlnto0 = -1 to do initial lineto, and =0 to do nothing.
*/

void _cpdf_arc_small(CPDFdoc *pdf, float x, float y, float r, float midtheta, float hangle, int mvlnto0, int ccwcw)
{
float vcos, vsin;		/* cos, sin of htheta */
float x0=0.0, y0=0.0, x1, y1, x2, y2, x3, y3;
float fccwcw;
float htheta;
    if(ccwcw < 0) fccwcw = -1.0;
    else fccwcw = 1.0;
    htheta = fabs(hangle)*PI/180.0;
    vcos = cos(htheta);
    vsin = sin(htheta);

    /* OK now about trying to allow things that PostScript arc does, but we need to supply choices. */
    x0 = r*vcos;
    y0 = -fccwcw*r*vsin;
    rotate_xyCoordinate(x0, y0, midtheta, &x0, &y0);
    if(mvlnto0 == 1)
	cpdf_rawMoveto(pdf, x+x0, y+y0);
    else if(mvlnto0 == -1)
	cpdf_rawLineto(pdf, x+x0, y+y0);

    x1 = x2 = r*(4.0 - vcos)/3.0;
    y1 = r*fccwcw *(1.0 - vcos) * (vcos - 3.0) / (3.0*vsin);
    y2 = -y1;
    x3 = r*vcos;
    y3 = fccwcw*r*vsin;
    rotate_xyCoordinate(x1, y1, midtheta, &x1, &y1);	/* rotate conrol points */
    rotate_xyCoordinate(x2, y2, midtheta, &x2, &y2);
    rotate_xyCoordinate(x3, y3, midtheta, &x3, &y3);
    cpdf_rawCurveto(pdf, x+x1, y+y1, x+x2, y+y2, x+x3, y+y3);
#ifdef CPDF_DEBUG
    fprintf(stderr, "midTheta= %g,  halfAngle= %g\n",  midtheta, hangle);
    fprintf(stderr, "x1,y1=%g, %g,   x2,y2=%g, %g,   x3,y3=%g, %g\n",
			x1, y1, x2, y2, x3, y3);
#endif
}


/* 
    Rotate (x, y) by angle in degrees, and return output in (*xrot, *yrot).
    xrot, yrot may point to the same input variables x and y, but in that case,
    x, y must be plain float variables, and NOT expressions.
*/
void rotate_xyCoordinate(float x, float y, float angle, float *xrot, float *yrot)
{
float rcos, rsin;		/* cos, sine of midtheta */
float xr, yr;			/* temps for here */
    angle *= PI/180.0;
    rcos = cos(angle);
    rsin = sin(angle);
    xr = rcos*x - rsin*y;
    yr = rsin*x + rcos*y;
    *xrot = xr;
    *yrot = yr;
}


