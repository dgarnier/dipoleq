/*
**
**    plot_contour.c
**
**    12/12/98
**
**    make a contour plot using the ClibPDF routines and
**    contour.c from the tokamac code.
**
**    (c) D. Garnier -- Columbia University
*/

#include <stdio.h>
#include "contour.h"
#include "cpdflib.h"
#include "plot_contour.h"

void do_contour_pt_open(double x, double y, double z, int ccode);
void do_contour_pt_closed(double x, double y, double z, int ccode);

static CPDFdoc *curPDF;
static int doFill=0;

void do_contour_pt_open(double x, double y, double dummy, int ccode)
{
	switch (ccode) {
		case CONTOUR_START:
/*			cpdf_newpath(curPDF);
*/			cpdf_moveto(curPDF,(float) x, (float) y);
			break;
		case CONTOUR_TRACE:
			cpdf_lineto(curPDF,(float) x, (float) y);
			break;
		case CONTOUR_STOP:
			cpdf_moveto(curPDF,(float) x, (float) y);
/*			if (doFill==FILL) {
				cpdf_fillAndStroke(curPDF);
			} else
				cpdf_stroke(curPDF);
*/			break;
		}
}
void do_contour_pt_closed(double x, double y, double dummy, int ccode)
{
	switch (ccode) {
		case CONTOUR_START:
/*			cpdf_newpath(curPDF);
*/			cpdf_moveto(curPDF,(float) x, (float) y);
			break;
		case CONTOUR_TRACE:
			cpdf_lineto(curPDF,(float) x, (float) y);
			break;
		case CONTOUR_STOP:
/*		    cpdf_closepath(curPDF);
			if (doFill==FILL) {
				cpdf_fillAndStroke(curPDF);
			} else
				cpdf_stroke(curPDF);
*/
			break;
		}
}


void  plot_contour(CPDFdoc *pdf, int flags, double *xarg, double *yarg, double **zarg, double level,
	  int nx1, int nx2, int ny1, int ny2)
{
	curPDF = pdf;
	doFill = flags & FILL;

	cpdf_newpath(pdf);

	contour(xarg, yarg, zarg, nx1, nx2, ny1, ny2, level,
			CONTOUR_ONLY_OPEN,!(flags & MIDPT), do_contour_pt_open);
	contour(xarg, yarg, zarg, nx1, nx2, ny1, ny2, level,
			CONTOUR_ONLY_CLOSED,!(flags & MIDPT), do_contour_pt_closed);

	cpdf_closepath(curPDF);
	if (doFill==FILL) {
		cpdf_fillAndStroke(curPDF);
	} else
		cpdf_stroke(curPDF);

}
