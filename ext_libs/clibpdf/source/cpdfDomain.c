/* cpdfDomain.c  -- creation and management of plot domains 
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf.
 *  or http://www.fastio.com/LICENSE.pdf

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


static char *defaultMeshDash = "[2 3] 0";	/* This is OK for now.  Can't change this via API */

CPDFplotDomain *cpdf_createPlotDomain(float x, float y, float w, float h,
			float xL, float xH, float yL, float yH,
			int xtype, int ytype, int polar)
{
CPDFplotDomain *newDomain;
int i;
float sL, sH;
	newDomain = (CPDFplotDomain *) malloc((size_t)sizeof(CPDFplotDomain));
	_cpdf_malloc_check((void *)newDomain);
        newDomain->magic = (unsigned int)DOMAIN_MAGIC_NUMBER;
	newDomain->xloc = x;
	newDomain->yloc = y;
	newDomain->width = w;
	newDomain->height = h;
	if(xtype == LINEAR) {
	    cpdf_suggestLinearDomainParams(xL, xH, &sL, &sH,
		&newDomain->xvalFirstMeshLinMajor, &newDomain->xmeshIntervalLinMajor,
		&newDomain->xvalFirstMeshLinMinor, &newDomain->xmeshIntervalLinMinor);
	}
	newDomain->xvalL = xL;
	newDomain->xvalH = xH;

	if(ytype == LINEAR) {
	    cpdf_suggestLinearDomainParams(yL, yH, &sL, &sH,
		&newDomain->yvalFirstMeshLinMajor, &newDomain->ymeshIntervalLinMajor,
		&newDomain->yvalFirstMeshLinMinor, &newDomain->ymeshIntervalLinMinor);
	}
	newDomain->yvalL = yL;
	newDomain->yvalH = yH;

	newDomain->xtype = xtype;
	newDomain->ytype = ytype;
	newDomain->polar = polar;	/* reserved: not implemented */
	newDomain->enableMeshMajor = 1;
	newDomain->enableMeshMinor = 1;
	newDomain->meshDashMajor = defaultMeshDash;
	newDomain->meshLineWidthMajor = 0.15;
	newDomain->meshDashMinor = defaultMeshDash;
	newDomain->meshLineWidthMinor = 0.15;
	for(i=0; i<3; i++) {
	    newDomain->meshLineColorMajor[i] = 0.0;
	    newDomain->meshLineColorMinor[i] = 0.0;
	}

	return newDomain;
}


CPDFplotDomain *cpdf_createTimePlotDomain(float x, float y, float w, float h,
			struct tm *xTL, struct tm *xTH, float yL, float yH,
			int xtype, int ytype, int polar)
{
CPDFplotDomain *newDomain;
int i;
float sL, sH;
	newDomain = (CPDFplotDomain *) malloc((size_t)sizeof(CPDFplotDomain));
	_cpdf_malloc_check((void *)newDomain);
        newDomain->magic = (unsigned int)DOMAIN_MAGIC_NUMBER;
	newDomain->xloc = x;
	newDomain->yloc = y;
	newDomain->width = w;
	newDomain->height = h;
	cpdf_mktime(xTL); cpdf_mktime(xTH);	/* This sets all fields of struct tm correctly */
	memcpy(&newDomain->xvTL, xTL, sizeof(struct tm));
	memcpy(&newDomain->xvTH, xTH, sizeof(struct tm));
	newDomain->xvalL = 0.0;			/* keep also float representation for date */
	newDomain->xvalH = tm_to_NumDays(xTL, xTH);	 /* in float number of days */

	if(ytype == LINEAR) {
	    cpdf_suggestLinearDomainParams(yL, yH, &sL, &sH,
		&newDomain->yvalFirstMeshLinMajor, &newDomain->ymeshIntervalLinMajor,
		&newDomain->yvalFirstMeshLinMinor, &newDomain->ymeshIntervalLinMinor);
	}
	newDomain->yvalL = yL;
	newDomain->yvalH = yH;
	newDomain->xtype = xtype;
	newDomain->ytype = ytype;
	newDomain->polar = polar;	/* reserved: not implemented */
	newDomain->enableMeshMajor = 1;
	newDomain->enableMeshMinor = 1;
	newDomain->meshDashMajor = defaultMeshDash;
	newDomain->meshLineWidthMajor = 0.16;
	newDomain->meshDashMinor = defaultMeshDash;
	newDomain->meshLineWidthMinor = 0.1;
	for(i=0; i<3; i++) {
	    newDomain->meshLineColorMajor[i] = 0.0;
	    newDomain->meshLineColorMinor[i] = 0.0;
	}

	return newDomain;
}


void cpdf_freePlotDomain(CPDFplotDomain *aDomain)
{
	if(aDomain) free(aDomain);
}

/* fill plot domain with a gray */
void cpdf_fillDomainWithGray(CPDFplotDomain *aDomain, float gray)
{
CPDFdoc *pdf;
	pdf = aDomain->pdfdoc;
	cpdf_gsave(pdf);
	cpdf_rawRect(pdf, aDomain->xloc, aDomain->yloc, aDomain->width, aDomain->height);
	cpdf_setgrayFill(pdf, gray);
	cpdf_fill(pdf);
	cpdf_grestore(pdf);
}

/* fill plot domain with a gray */
void cpdf_fillDomainWithRGBcolor(CPDFplotDomain *aDomain, float r, float g, float b)
{
CPDFdoc *pdf;
	pdf = aDomain->pdfdoc;
	cpdf_gsave(pdf);
	cpdf_rawRect(pdf, aDomain->xloc, aDomain->yloc, aDomain->width, aDomain->height);
	cpdf_setrgbcolorFill(pdf, r, g, b);
	cpdf_fill(pdf);
	cpdf_grestore(pdf);
}


/* clip drawing at domain border. Bracket this call with gsave and grestore */
void cpdf_clipDomain(CPDFplotDomain *aDomain)
{
CPDFdoc *pdf;
	pdf = aDomain->pdfdoc;
	cpdf_newpath(pdf);
	cpdf_rawRect(pdf, aDomain->xloc, aDomain->yloc, aDomain->width, aDomain->height);
	cpdf_clip(pdf);
	cpdf_newpath(pdf);
}


/* This is the only Plot Domain related function that needs the PDF doc context argument.
   All other domain functions and Axis functions will get PDF doc context via what
   you specify here.
*/
/* ## MULTI_THREADING_OK ## */
CPDFplotDomain *cpdf_setPlotDomain(CPDFdoc *pdf, CPDFplotDomain *aDomain)
{
CPDFplotDomain *oldDomain;
	if(pdf == NULL || aDomain == NULL) return NULL;

	oldDomain = pdf->currentDomain;
	pdf->currentDomain = aDomain;
	aDomain->pdfdoc = pdf;		/* record the PDF document that this domain belongs to */
	pdf->x2points = (aDomain->width)/((aDomain->xvalH) - (aDomain->xvalL));
	pdf->y2points = (aDomain->height)/((aDomain->yvalH) - (aDomain->yvalL));
	if(aDomain->xtype == LOGARITHMIC) {
	    /* logarithmic */
	    pdf->xLlog = log10((double)(aDomain->xvalL));
	    pdf->xHlog = log10((double)(aDomain->xvalH));
	}
	if(aDomain->ytype == LOGARITHMIC) {
	    /* logarithmic */
	    pdf->yLlog = log10((double)(aDomain->yvalL));
	    pdf->yHlog = log10((double)(aDomain->yvalH));
	}
	return(oldDomain);
}

void cpdf_setMeshColor(CPDFplotDomain *aDomain, float meshMajorR, float meshMajorG, float meshMajorB,
						float meshMinorR, float meshMinorG, float meshMinorB)
{
	aDomain->meshLineColorMajor[0] = meshMajorR;
	aDomain->meshLineColorMajor[1] = meshMajorG;
	aDomain->meshLineColorMajor[2] = meshMajorB;
	aDomain->meshLineColorMinor[0] = meshMinorR;
	aDomain->meshLineColorMinor[1] = meshMinorG;
	aDomain->meshLineColorMinor[2] = meshMinorB;
}

void cpdf_drawMeshForDomain(CPDFplotDomain *aDomain)
{
CPDFdoc *pdf;
	pdf = aDomain->pdfdoc;
	cpdf_gsave(pdf);
	cpdf_setstrokeadjust(pdf, 1);	/* valid only in PDF-1.2 and later */
	_do_meshLines_X(aDomain);
	_do_meshLines_Y(aDomain);
	cpdf_grestore(pdf);
}

/* ToDo: consolidate X and Y versions of mesh functions below into one */

void _do_meshLines_X(CPDFplotDomain *aDomain)
{
float vL, vH, v, xp;
int m, ep;
int mL, expL, mH, expH;	/* mantissa MSD and exponents */
float fndays;
int minorBump = 1, majorBump = 1;
int minorBumpVar = MINUTE, majorBumpVar = HOUR;
struct tm vtm;
CPDFdoc *pdf;

    pdf = aDomain->pdfdoc;
    vL = aDomain->xvalL;
    vH = aDomain->xvalH;

    if(aDomain->xtype == LOGARITHMIC) {
	/* Log X domain */
	mL = (int)getMantissaExp(vL*1.0001, &expL);
	mH = (int)getMantissaExp(vH*1.0001, &expH);
	/* fprintf(stderr, "%dE%d -- %dE%d\n", mL, expL, mH, expH); */
	/* Do major and minor mesh lines as we go. */
	for(m=mL, ep=expL; (v = (float)m * pow(10.0, (double)ep)) <= vH*1.0001; ) {
	    /* fprintf(stderr, "tick %g\n", v); */
	    if(m == 1) {
		/* Major mesh lines at m == 1 */
	        cpdf_setlinewidth(pdf, aDomain->meshLineWidthMajor);
		cpdf_setdash(pdf, aDomain->meshDashMajor);
		cpdf_setrgbcolorStroke(pdf, aDomain->meshLineColorMajor[0],
			aDomain->meshLineColorMajor[1], aDomain->meshLineColorMajor[2]);
	    }
	    else {
		/* Minor mesh lines */
	        cpdf_setlinewidth(pdf, aDomain->meshLineWidthMinor);
		cpdf_setdash(pdf, aDomain->meshDashMinor);
		cpdf_setrgbcolorStroke(pdf, aDomain->meshLineColorMinor[0],
			aDomain->meshLineColorMinor[1], aDomain->meshLineColorMinor[2]);
	    }
	    xp = x_Domain2Points(pdf, v);
	    cpdf_rawMoveto(pdf, xp, aDomain->yloc);
	    cpdf_rawLineto(pdf, xp, aDomain->yloc + aDomain->height);
	    cpdf_stroke(pdf);
	    m++;
	    if(m >= 10) {	/* cycle m through 1 .. 9 */
		m = 1;
		ep++;
	    }
	} /* end for(m=mL,... ) */
    }
    else if(aDomain->xtype == LINEAR) {
	/* Linear X domain */
	/* do minor mesh lines */
	cpdf_setlinewidth(pdf, aDomain->meshLineWidthMinor);
	cpdf_setdash(pdf, aDomain->meshDashMinor);
	cpdf_setrgbcolorStroke(pdf, aDomain->meshLineColorMinor[0],
			aDomain->meshLineColorMinor[1], aDomain->meshLineColorMinor[2]);
	for(v = aDomain->xvalFirstMeshLinMinor; v <= vH*1.0001; v += aDomain->xmeshIntervalLinMinor) {
	    xp = x_Domain2Points(pdf, v);
	    cpdf_rawMoveto(pdf, xp, aDomain->yloc);
	    cpdf_rawLineto(pdf, xp, aDomain->yloc + aDomain->height);
	    cpdf_stroke(pdf);
	}
	/* do minor mesh lines */
	cpdf_setlinewidth(pdf, aDomain->meshLineWidthMajor);
	cpdf_setdash(pdf, aDomain->meshDashMajor);
	cpdf_setrgbcolorStroke(pdf, aDomain->meshLineColorMajor[0],
			aDomain->meshLineColorMajor[1], aDomain->meshLineColorMajor[2]);
	for(v = aDomain->xvalFirstMeshLinMajor; v <= vH*1.0001; v += aDomain->xmeshIntervalLinMajor) {
	    xp = x_Domain2Points(pdf, v);
	    cpdf_rawMoveto(pdf, xp, aDomain->yloc);
	    cpdf_rawLineto(pdf, xp, aDomain->yloc + aDomain->height);
	    cpdf_stroke(pdf);
	}
    }
    else if(aDomain->xtype == TIME) {
	/* Linear TIME X domain */
	/* do minor mesh lines */
	fndays = tm_to_NumDays(&aDomain->xvTL, &aDomain->xvTH);
	_setDefaultTimeBumpVar(fndays, &minorBumpVar, &majorBumpVar, &minorBump, &majorBump);
	/* do minor ticks */
	cpdf_setlinewidth(pdf, aDomain->meshLineWidthMinor);
	cpdf_setdash(pdf, aDomain->meshDashMinor);
	cpdf_setrgbcolorStroke(pdf, aDomain->meshLineColorMinor[0],
			aDomain->meshLineColorMinor[1], aDomain->meshLineColorMinor[2]);
	memcpy(&vtm, &aDomain->xvTL, sizeof(struct tm));
	for(v = tm_to_NumDays(&aDomain->xvTL, &vtm) ; v <= fndays*1.0001 ;
			v = _bump_tm_Time(&aDomain->xvTL, &vtm, minorBumpVar, minorBump)) {
	    xp = x_Domain2Points(pdf, v);
	    cpdf_rawMoveto(pdf, xp, aDomain->yloc);
	    cpdf_rawLineto(pdf, xp, aDomain->yloc + aDomain->height);
	    cpdf_stroke(pdf);
	}
	/* do minor mesh lines */
	cpdf_setlinewidth(pdf, aDomain->meshLineWidthMajor);
	cpdf_setdash(pdf, aDomain->meshDashMajor);
	cpdf_setrgbcolorStroke(pdf, aDomain->meshLineColorMajor[0],
			aDomain->meshLineColorMajor[1], aDomain->meshLineColorMajor[2]);
	memcpy(&vtm, &aDomain->xvTL, sizeof(struct tm));
	for(v = tm_to_NumDays(&aDomain->xvTL, &vtm) ; v <= fndays*1.0001 ;
			v = _bump_tm_Time(&aDomain->xvTL, &vtm, majorBumpVar, majorBump)) {
	    xp = x_Domain2Points(pdf, v);
	    cpdf_rawMoveto(pdf, xp, aDomain->yloc);
	    cpdf_rawLineto(pdf, xp, aDomain->yloc + aDomain->height);
	    cpdf_stroke(pdf);
	}
    } /* end the else part of: if(aDomain-> xtype) {} else {} */
}

void _do_meshLines_Y(CPDFplotDomain *aDomain)
{
float vL, vH, v, yp;
int m, ep;
int mL, expL, mH, expH;	/* mantissa MSD and exponents */
CPDFdoc *pdf;

    pdf = aDomain->pdfdoc;
    vL = aDomain->yvalL;
    vH = aDomain->yvalH;

    if(aDomain->ytype == LOGARITHMIC) {
	/* Log Y domain */
	mL = (int)getMantissaExp(vL*1.0001, &expL);
	mH = (int)getMantissaExp(vH*1.0001, &expH);
	/* fprintf(stderr, "%dE%d -- %dE%d\n", mL, expL, mH, expH); */
	/* Do major and minor mesh lines as we go. */
	for(m=mL, ep=expL; (v = (float)m * pow(10.0, (double)ep)) <= vH*1.0001; ) {
	    /* fprintf(stderr, "tick %g\n", v); */
	    if(m == 1) {
		/* Major mesh lines at m == 1 */
	        cpdf_setlinewidth(pdf, aDomain->meshLineWidthMajor);
		cpdf_setdash(pdf, aDomain->meshDashMajor);
		cpdf_setrgbcolorStroke(pdf, aDomain->meshLineColorMajor[0],
			aDomain->meshLineColorMajor[1], aDomain->meshLineColorMajor[2]);
	    }
	    else {
		/* Minor mesh lines */
	        cpdf_setlinewidth(pdf, aDomain->meshLineWidthMinor);
		cpdf_setdash(pdf, aDomain->meshDashMinor);
		cpdf_setrgbcolorStroke(pdf, aDomain->meshLineColorMinor[0],
			aDomain->meshLineColorMinor[1], aDomain->meshLineColorMinor[2]);
	    }
	    yp = y_Domain2Points(pdf, v);
	    cpdf_rawMoveto(pdf, aDomain->xloc, yp);
	    cpdf_rawLineto(pdf, aDomain->xloc + aDomain->width, yp);
	    cpdf_stroke(pdf);
	    m++;
	    if(m >= 10) {	/* cycle m through 1 .. 9 */
		m = 1;
		ep++;
	    }
	} /* end for(m=mL,... ) */
    }
    else {
	/* Linear Y domain */
	/* do minor mesh lines */
	cpdf_setlinewidth(pdf, aDomain->meshLineWidthMinor);
	cpdf_setdash(pdf, aDomain->meshDashMinor);
	cpdf_setrgbcolorStroke(pdf, aDomain->meshLineColorMinor[0],
			aDomain->meshLineColorMinor[1], aDomain->meshLineColorMinor[2]);
	for(v = aDomain->yvalFirstMeshLinMinor; v <= vH*1.0001; v += aDomain->ymeshIntervalLinMinor) {
	    yp = y_Domain2Points(pdf, v);
	    cpdf_rawMoveto(pdf, aDomain->xloc, yp);
	    cpdf_rawLineto(pdf, aDomain->xloc + aDomain->width, yp);
	    cpdf_stroke(pdf);
	}
	/* do major mesh lines */
	cpdf_setlinewidth(pdf, aDomain->meshLineWidthMajor);
	cpdf_setdash(pdf, aDomain->meshDashMajor);
	cpdf_setrgbcolorStroke(pdf, aDomain->meshLineColorMajor[0],
			aDomain->meshLineColorMajor[1], aDomain->meshLineColorMajor[2]);
	for(v = aDomain->yvalFirstMeshLinMajor; v <= vH*1.0001; v += aDomain->ymeshIntervalLinMajor) {
	    yp = y_Domain2Points(pdf, v);
	    cpdf_rawMoveto(pdf, aDomain->xloc, yp);
	    cpdf_rawLineto(pdf, aDomain->xloc + aDomain->width, yp);
	    cpdf_stroke(pdf);
	}
    } /* end the else part of: if(aDomain-> xtype) {} else {} */
}


void cpdf_setLinearMeshParams(CPDFplotDomain *aDomain, int xy, float mesh1ValMajor, float intervalMajor,
					      float mesh1ValMinor, float intervalMinor)
{
    if(xy == X_MESH) {
	aDomain->xvalFirstMeshLinMajor = mesh1ValMajor;	/* value of first major mesh line */
	aDomain->xvalFirstMeshLinMinor = mesh1ValMinor;	/* value of first minor mesh line */
	aDomain->xmeshIntervalLinMajor = intervalMajor;	/* mesh interval for linear axis */
	aDomain->xmeshIntervalLinMinor = intervalMinor;
    }
    else {
	aDomain->yvalFirstMeshLinMajor = mesh1ValMajor;	/* value of first major mesh line */
	aDomain->yvalFirstMeshLinMinor = mesh1ValMinor;	/* value of first minor mesh line */
	aDomain->ymeshIntervalLinMajor = intervalMajor;	/* mesh interval for linear axis */
	aDomain->ymeshIntervalLinMinor = intervalMinor;
    }
}


/* Thse must be called after cpdf_setPlotDomain() has been called */

float x_Domain2Points(CPDFdoc *pdf, float x)
{
float xrval = 0.0;
double xvlog = 0.0, fraction=0.0;
    if(!(pdf->currentDomain)) {
	cpdf_Error(pdf, 2, "ClibPDF", "Plot domain has not been set. See cpdf_setPlotDomain()");
	return(0.0);
    }

    xrval = pdf->currentDomain->xloc;	/* min point value for X */
    if(pdf->currentDomain->xtype == LOGARITHMIC) {
	/* logarithmic */
	if(x > 0.0) {
	    xvlog = log10((double)x);
	    fraction = (xvlog - pdf->xLlog)/(pdf->xHlog - pdf->xLlog);
		/* xHlog, xLlog precomputed in cpdf_setPlotDomain() */
	    xrval += (pdf->currentDomain->width) * fraction;
	}
    }
    else {
	/* linear */
	xrval += (x - (pdf->currentDomain->xvalL)) * pdf->x2points;
    }
    return xrval;
}


float y_Domain2Points(CPDFdoc *pdf, float y)
{
float yrval = 0.0;
double yvlog = 0.0, fraction=0.0;
    if(!(pdf->currentDomain)) {
	cpdf_Error(pdf, 1, "ClibPDF", "Plot domain has not been set. See cpdf_setPlotDomain()");
	return(0.0);
    }

    yrval = pdf->currentDomain->yloc;	/* min point value for Y */
    if(pdf->currentDomain->ytype == LOGARITHMIC){
	/* logarithmic */
	if(y > 0.0) {
	    yvlog = log10((double)y);
	    fraction = (yvlog - pdf->yLlog)/(pdf->yHlog - pdf->yLlog);
		/* yHlog, yLlog precomputed in cpdf_setPlotDomain() */
	    yrval += (pdf->currentDomain->height) * fraction;
	}
    }
    else {
	/* linear */
	yrval += (y - (pdf->currentDomain->yvalL)) * pdf->y2points;
    }
    return yrval;
}

/* Convert point based Y location back to domain unit */
float y_Points2Domain(CPDFdoc *pdf, float ypt)
{
float yrval = 0.0;	/* return value */

    if(!(pdf->currentDomain)) {
	cpdf_Error(pdf, 1, "ClibPDF", "Plot domain has not been set. See cpdf_setPlotDomain()");
	return(0.0);
    }

    if(pdf->currentDomain->ytype == LOGARITHMIC){
	/* logarithmic */
	double yvlog, fraction;
	fraction = (double)(ypt - pdf->currentDomain->yloc) / (double)(pdf->currentDomain->height);
	yvlog = (pdf->yHlog - pdf->yLlog) * fraction + pdf->yLlog;
	/* yHlog, yLlog precomputed in cpdf_setPlotDomain() */
        yrval = pow((double)10.0, yvlog);
    }
    else {
	/* linear domain */
	yrval = (ypt - pdf->currentDomain->yloc) / pdf->y2points + pdf->currentDomain->yvalL;
    }
    return yrval;
}

/* Convert point based X location back to domain unit */
float x_Points2Domain(CPDFdoc *pdf, float xpt)
{
float xrval = 0.0;	/* return value */

    if(!(pdf->currentDomain)) {
	cpdf_Error(pdf, 1, "ClibPDF", "Plot domain has not been set. See cpdf_setPlotDomain()");
	return(0.0);
    }

    if(pdf->currentDomain->xtype == LOGARITHMIC){
	/* logarithmic */
	double xvlog, fraction;
	fraction = (double)(xpt - pdf->currentDomain->xloc) / (double)(pdf->currentDomain->width);
	xvlog = (pdf->xHlog - pdf->xLlog) * fraction + pdf->xLlog;
	/* xHlog, xLlog precomputed in cpdf_setPlotDomain() */
        xrval = pow((double)10.0, xvlog);
    }
    else {
	/* linear domain */
	xrval = (xpt - pdf->currentDomain->xloc) / pdf->x2points + pdf->currentDomain->xvalL;
    }
    return xrval;
}

