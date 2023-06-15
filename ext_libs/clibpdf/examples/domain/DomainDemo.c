/* DomainDemo.c -- Demo program for LogLog and LinearLinear plot domains.
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

cc -Wall -o DomainDemo DomainDemo.c -lcpdf

1999-08-22 [io] for v2.00
1998-09-09 [io]
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "cpdflib.h"


/* Some plot data. In a real program, these come from external sources. */
/* --------------------------------------------------------------------------- */
static float ymax = 100.0;
static float ymaxL = 80.0;
static float ymin = 0.01;
static float  xmax = 3.0;
static float  xmin = 0.03;
/* Fit parameter definition for Gaussian:  y = A*exp( - ((x - X0)/Sigma)^2)
 */
static float A = 70.5;
static float X0 = 1.2;
static float Sigma = 0.4;

/* computed Gaussian data */
static float x[120], y[120];

/* --------------------------------------------------------------------------- */
extern void do_LogLog(CPDFdoc *pdf);
extern void do_LinearLinear(CPDFdoc *pdf);
extern void plot_Curve(CPDFdoc *pdf);
extern void generateGaussianData(void);
void do_myTitle(CPDFdoc *pdf);
/* --------------------------------------------------------------------------- */


int main(int argc, char *argv[])
{
CPDFdoc *pdf;
char *bufPDF;
int length= 0;
	pdf = cpdf_open(0, NULL);
	cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
	cpdf_init(pdf);
	cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */

	cpdf_setCreator(pdf, "ClibPDF Domain Demo: DomainDemo.c");

	/* Your drawing code in PDF below */
	generateGaussianData();
	do_LinearLinear(pdf);
	do_LogLog(pdf);
	do_myTitle(pdf);

	cpdf_finalizeAll(pdf);			/* PDF file is actually written here */

	/* if you wish to use memory-resident PDF stream do as follows:
	   bufPDF points to the beginning of buffer containing PDF, and length is a byte count.
	*/
	bufPDF = cpdf_getBufferForPDF(pdf, &length);
	printf("PDF buffer length= %d\n", length);

	/* Or you can save generated PDF to file */
	cpdf_savePDFmemoryStreamToFile(pdf, "domain.pdf");

	cpdf_launchPreview(pdf);	/* launch Acrobat/PDF viewer on the output file */
	cpdf_close(pdf);		/* shut down the library resources */
	return(0);
}

void generateGaussianData(void)
{
int i;
float xp = 0.03;
float xbump = 0.03;
	for(i=0; i<100; i++) {
	    /* A Gaussian curve */
	    x[i] = xp;
	    y[i] = A * exp( - ((xp - X0) / Sigma) * ((xp - X0) / Sigma));
	    xp += xbump;
	}
}

/* This will be called from both do_LinearLinear(), and do_LogLog().
   Respective domains will take care of transformation and scaling.
   So, we just use the raw data values.
*/

void plot_Curve(CPDFdoc *pdf)
{
int i;
	cpdf_comments(pdf, "\n%% plot the curve (linear-linear).\n");
	cpdf_setlinewidth(pdf, 2.0);
	cpdf_setrgbcolorStroke(pdf, 0.0, 0.0, 1.0);	/* blue */
	for(i=0; i<100; i++) {
	    if(i==0) cpdf_moveto(pdf, x[i], y[i]);
	    else     cpdf_lineto(pdf, x[i], y[i]);
	}
	cpdf_stroke(pdf);
}


void do_LinearLinear(CPDFdoc *pdf)
{
CPDFplotDomain *myDomain, *oldDomain;
CPDFaxis *xAxis, *yAxis;

	myDomain = cpdf_createPlotDomain( 1.5*inch, 1.5*inch, 4.5*inch, 3.0*inch,
				0.0, xmax, 0.0, ymaxL, LINEAR, LINEAR, 0);
	oldDomain = cpdf_setPlotDomain(pdf, myDomain);	/* save oldDomain for later restore */
	cpdf_fillDomainWithRGBcolor(myDomain, 1.0, 0.9, 0.9);	/* light pink */
	cpdf_drawMeshForDomain(myDomain);

	cpdf_setgray(pdf, 0.0);
	/* X-Axis --------------------------------------------------------------------------------- */
	xAxis = cpdf_createAxis( 0.0, 4.5*inch, LINEAR, 0.0, xmax);
	cpdf_attachAxisToDomain(xAxis, myDomain, 0.0, -0.2*inch);
	cpdf_setAxisNumberFormat(xAxis, "%g", "Helvetica", 16.0);
	cpdf_setAxisLabel(xAxis, "X axis (linear)", "Times-Roman", "MacRomanEncoding", 20.0);
	cpdf_drawAxis(xAxis);
	cpdf_freeAxis(xAxis);

	/* Y-Axis --------------------------------------------------------------------------------- */
	yAxis = cpdf_createAxis(90.0, 3.0*inch, LINEAR, 0.0, ymaxL);
	cpdf_attachAxisToDomain(yAxis, myDomain, -0.2*inch, 0.0);
	cpdf_setAxisNumberFormat(yAxis, "%2.f", "Helvetica", 16.0);
	cpdf_setAxisLabel(yAxis, "Y axis (linear)", "Times-Roman", "MacRomanEncoding", 20.0);
	cpdf_drawAxis(yAxis);
	cpdf_freeAxis(yAxis);

	/* Plot it */
	plot_Curve(pdf);

	cpdf_setPlotDomain(pdf, oldDomain);		/* restore previous plot domain */
	cpdf_freePlotDomain(myDomain);		/* deallocate the plot domain */
}


void do_LogLog(CPDFdoc *pdf)
{
CPDFplotDomain *myDomain, *oldDomain;
CPDFaxis *xAxis, *yAxis;

	myDomain = cpdf_createPlotDomain( 1.5*inch, 6.0*inch, 4.5*inch, 3.0*inch,
				xmin, xmax, ymin, ymax, LOGARITHMIC, LOGARITHMIC, 0);
	oldDomain = cpdf_setPlotDomain(pdf, myDomain);	/* save oldDomain for later restore */
	cpdf_fillDomainWithRGBcolor(myDomain, 1.0, 0.9, 0.9);	/* light pink */
	cpdf_drawMeshForDomain(myDomain);

	cpdf_setgray(pdf, 0.0);
	/* X-Axis --------------------------------------------------------------------------------- */
	xAxis = cpdf_createAxis( 0.0, 4.5*inch, LOGARITHMIC, xmin, xmax);
	cpdf_setLogAxisNumberSelector(xAxis, LOGAXSEL_13);	/* Nums 1 and 3s */
	cpdf_attachAxisToDomain(xAxis, myDomain, 0.0, -0.2*inch);
	cpdf_setAxisNumberFormat(xAxis, "%g", "Helvetica", 16.0);
	cpdf_setAxisLabel(xAxis, "X axis (logarithmic)", "Times-Roman", "MacRomanEncoding", 20.0);
	cpdf_drawAxis(xAxis);
	cpdf_freeAxis(xAxis);

	/* Y-Axis --------------------------------------------------------------------------------- */
	yAxis = cpdf_createAxis(90.0, 3.0*inch, LOGARITHMIC, ymin, ymax);
	cpdf_attachAxisToDomain(yAxis, myDomain, -0.2*inch, 0.0);
	cpdf_setAxisNumberFormat(yAxis, "%g", "Helvetica", 16.0);
	cpdf_setAxisLabel(yAxis, "Y axis (logarithmic)", "Times-Roman", "MacRomanEncoding", 20.0);
	cpdf_drawAxis(yAxis);
	cpdf_freeAxis(yAxis);

	/* We have to do clipping of domain, because log plot can under-shoot */
	cpdf_gsave(pdf);				/* to later undo clipping */
	cpdf_clipDomain(myDomain);		/* clip to myDomain: TRY it without this line */
	plot_Curve(pdf);				/* Plot it */
	cpdf_grestore(pdf);			/* remove clipping */

	cpdf_setPlotDomain(pdf, oldDomain);		/* restore previous plot domain */
	cpdf_freePlotDomain(myDomain);		/* deallocate the plot domain */
}


void do_myTitle(CPDFdoc *pdf)
{
	cpdf_setrgbcolorFill(pdf, 0.0, 0.0, 1.0);		/* Blue */
	cpdf_beginText(pdf, 0);
	cpdf_setFont(pdf, "Times-Italic", "MacRomanEncoding", 18.0);
	cpdf_text(pdf, 1.5, 9.7, 0.0, "Y = A * exp( - ((X-Xo)/S)^2)");
	cpdf_text(pdf, 1.7, 9.4, 0.0, "where  A=70,  Xo=1.2,  S=0.4");
	cpdf_endText(pdf);
}

