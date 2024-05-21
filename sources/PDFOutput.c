/*
**
**    PDFOutput.c
**
**    12/12/98
**
**    make a contour plot using the ClibPDF routines and
**    contour.c from the tokamac code.
**
**    (c) D. Garnier -- Columbia University
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <cpdflib.h>

#include "tokamak.h"
#include "contour.h"
#include "plot_contour.h"
#include "PDFOutput.h"
#include "multitask.h"
#include "nrutil.h"

#include <unistd.h>

#define APPNAME  "DipolEq v0.9"
void    Eq_Plot(CPDFdoc *pdf,TOKAMAK *td);
void	showGeometry(CPDFdoc *pdf,TOKAMAK *td);
void 	contour_Psi(CPDFdoc *pdf,TOKAMAK *td);
void 	contour_ModB(CPDFdoc *pdf,TOKAMAK *td);
void    contour_IsPlasma(CPDFdoc *pdf,TOKAMAK *td);


static const int    nB2lev = 24;
static const double B2_Levels[] = {1.0e4,   2.5e3,   4.0e2,
							       1.0e2,   2.5e1,   4.0e0,
							       1.0e0,   2.5e-1,  4.0e-2,
							       1.0e-2,  2.5e-3,  4.0e-4,
							       1.0e-4,  2.5e-5,  4.0e-4,
							       1.0e-6,  2.5e-7,  4.0e-8,
							       1.0e-8,  2.5e-9,  4.0e-10,
							       1.0e-10, 2.5e-11, 4.0e-12};

static const char *B2Dash[] = { "[] 0",
								"[4 4] 0"
								};

static const int B2DI[] = {0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1,
                           0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1};

static CPDFviewerPrefs theVP = { PM_NONE, NO, NO, NO, NO, NO, PL_SINGLE, PM_NONE };

#define BSIZE 8191
void PDFOutput(TOKAMAK *td)
{
    CPDFdoc *pdf;

    char  buf[BSIZE+1];
    pdf = cpdf_open(0, NULL);

    snprintf(buf, BSIZE, "%s.pdf", td->Oname);
    cpdf_setOutputFilename(pdf, buf);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);

    /* This will have outline (book marks) visible upon document open */
    cpdf_setViewerPreferences(pdf,&theVP);

    /* PDF Info object */
	snprintf(buf, BSIZE, "%s: built %s, %s", APPNAME, __DATE__, __TIME__);
	cpdf_setCreator(pdf, buf);
	snprintf(buf, BSIZE, "%s Equilibrium Output", td->Name);
	cpdf_setTitle(pdf, buf);
    cpdf_setSubject(pdf,td->Info);
	snprintf(buf, BSIZE, "Equilibrium, Grad-Shafranov,  %s, %s, %s", APPNAME,
			 td->Name, td->Info);
	cpdf_setKeywords(pdf, buf);

    cpdf_pageInit(pdf,1, PORTRAIT, LETTER, LETTER);		/* page orientation */
    Eq_Plot(pdf,td);

    cpdf_finalizeAll(pdf);			/* PDF file is actually written here */
	if (isatty(STDOUT_FILENO)) {	/* only launch preview if output is interactive terminal */
	    cpdf_launchPreview(pdf);	/* launch Acrobat/PDF viewer on the output file */
	}
    cpdf_close(pdf);		/* shut down the library resources */

}


void	showGeometry(CPDFdoc *pdf,TOKAMAK *td)
{
	LIMITER *l;
	COIL *c;
	SUBCOIL *sc;
	MEAS *m;
	int i,j;

	cpdf_gsave(pdf);

	/* do limiters */

	cpdf_setrgbcolorStroke(pdf,1.0, 0.0, 0.0);	/* red */


	for (i=0;i<td->NumLimiters;i++) {
		l = td->Limiters[i];
		if (l->Enabled != Limiter_Off) {
			cpdf_newpath(pdf);
			cpdf_moveto(pdf,l->X1,l->Z1);
			cpdf_lineto(pdf,l->X2,l->Z2);
			cpdf_stroke(pdf);
		}
	}

	/* do coils .... */
	cpdf_setgray(pdf,1.0);
	cpdf_setrgbcolorFill(pdf,0.0, 1.0, 0.0);	/* green */
	cpdf_setrgbcolorStroke(pdf,0.2, 0.2, 0.2);	/* gray */

	for (i=0;i<td->NumCoils;i++) {
		c = td->Coils[i];
		if (c->Enabled >= Coil_On) {

			if (c->dX > 0) {
				/* make box */
				cpdf_newpath(pdf);
				cpdf_moveto(pdf,c->X-c->dX/2,c->Z-c->dZ/2);
				cpdf_lineto(pdf,c->X-c->dX/2,c->Z+c->dZ/2);
				cpdf_lineto(pdf,c->X+c->dX/2,c->Z+c->dZ/2);
				cpdf_lineto(pdf,c->X+c->dX/2,c->Z-c->dZ/2);
				cpdf_lineto(pdf,c->X-c->dX/2,c->Z-c->dZ/2);
				cpdf_fillAndStroke(pdf);
				/* put X in it */
/*				cpdf_newpath(pdf);
				cpdf_moveto(pdf,c->X-c->dX/2,c->Z-c->dZ/2);
				cpdf_lineto(pdf,c->X+c->dX/2,c->Z+c->dZ/2);
				cpdf_moveto(pdf,c->X+c->dX/2,c->Z-c->dZ/2);
				cpdf_lineto(pdf,c->X-c->dX/2,c->Z+c->dZ/2);
				cpdf_stroke(pdf);
*/
				/* mini grid points... */
//				for (j=0; j<c->NumSubCoils; j++) {
//					sc = c->SubCoils[j];
//					cpdf_marker(pdf,sc->X,sc->Z,0,.5);
//				}
			} else {
				for (j=0; j<c->NumSubCoils; j++) {
					sc = c->SubCoils[j];
					cpdf_marker(pdf,sc->X,sc->Z,3,2);
				}
			}
		}
	}

	/* do measurements */
	cpdf_setrgbcolorStroke(pdf,1.0, 0.0, 1.0);	/* magenta */
	cpdf_setrgbcolorFill(pdf,0.5, 0.5, 0.5);	/* gray */

#define ARROW 8.0
#define DEG2RAD 0.017453

	for (i=0; i<td->NumMeasures;i++) {
		m = td->Measures[i];
		switch (m->mType) {
			case meas_bp:
			case meas_bpangle:
			case meas_bangle:  /* do something with an angle */
				cpdf_newpath(pdf);
				cpdf_moveto(pdf,m->X,m->Z);
				cpdf_rawRlineto(pdf,ARROW*sin(m->parm.bp.Angle*DEG2RAD),ARROW*cos(m->parm.bp.Angle*DEG2RAD));
				cpdf_stroke(pdf);
				break;
			case meas_flux:
				cpdf_marker(pdf,m->X,m->Z,10,5);
				break;
		}
	}
	cpdf_grestore(pdf);
}

void contour_Psi(CPDFdoc *pdf,TOKAMAK *td)
{
	PSIGRID   *pg;
	PLASMA    *pl;
	int		nmax,npsi;

	int			i;
	double	 	dPsi,PsiAxis,Psi;

	pg = td->PsiGrid;
	pl = td->Plasma;
	nmax = pg->Nsize;


    PsiAxis = pg->PsiAxis;
    npsi = pl->NumPsiPts/3;
    dPsi = pl->PsiXmax*(pg->PsiLim-PsiAxis)/npsi;

	for (i=npsi;i>=0;i--) { // plot every third psi point....
    MULTI;

	    Psi = i*dPsi+PsiAxis;
		plot_contour(pdf, NOFILL, pg->X, pg->Z, pg->Psi, Psi, 0, nmax, 0, nmax);

	}

	cpdf_setdash(pdf,"[1 2] 0");
	plot_contour(pdf, NOFILL, pg->X, pg->Z, pg->Psi, pg->PsiLim, 0, nmax, 0, nmax);
	cpdf_nodash(pdf);
}

void contour_ModB(CPDFdoc *pdf,TOKAMAK *td)
{
	PSIGRID   *pg;
	PLASMA    *pl;
	double   **B2;
	int			ix, iz, i, nmax;
	double      BB, maxB2, minB2;

	pg = td->PsiGrid;
	pl = td->Plasma;
	nmax = pg->Nsize;
	B2 = pl->B2;

	maxB2 = minB2 = B2[nmax/2][nmax/2];
	for (ix=1; ix<nmax-1; ix++)
		for (iz=1; iz<nmax-1 ; iz++) {
			BB = B2[ix][iz];
			if (BB < minB2) minB2 = BB;
			if (BB > maxB2) maxB2 = BB;
		}

	for (i = 0; i < nB2lev ; i++) {
MULTI;
		if ((BB = B2_Levels[i]) < minB2) break;
		/* contour B2, but ignore boundary edges */
		if (BB <= maxB2) {
			cpdf_setdash(pdf,B2Dash[B2DI[i]]);
			plot_contour(pdf, NOFILL, pg->X, pg->Z, B2, BB, 1, nmax-1, 1, nmax-1);
//			plot_contour(pg->X, pg->Z, B2, BB, 0, nmax, 0, nmax);
		}

	}
}


void contour_IsPlasma(CPDFdoc *pdf,TOKAMAK *td) {

	PSIGRID *pg;
	int		nmax,ix,iz;
	double	**ipd;

	pg = td->PsiGrid;

	nmax = pg->Nsize;
	ipd = dmatrix(0,nmax,0,nmax);
	for (ix=0;ix<=nmax;ix++)
		for (iz=0;iz<=nmax;iz++)
			ipd[ix][iz] = (double) pg->IsPlasma[ix][iz];
	plot_contour(pdf, FILL | MIDPT , pg->X, pg->Z, ipd, 0.05, 0, nmax, 0, nmax);
	free_dmatrix(ipd,0,nmax,0,nmax);
}


void    Eq_Plot(CPDFdoc *pdf,TOKAMAK *td)
{

	PSIGRID   		*pg;

	CPDFplotDomain 	*myDomain, *oldDomain;
	CPDFaxis 		*xAxis, *yAxis;

	float 			aspect, xmin, xmax, ymin, ymax, plx, ply;

	pg = td->PsiGrid;
	xmin = (float) pg->Xmin;
	xmax = (float) pg->Xmax;
	ymin = (float) pg->Zmin;
	ymax = (float) pg->Zmax;
	aspect = (float) pg->dz/pg->dx;
	plx = 6.5*inch;
	ply = plx*aspect;
	if (ply > 9.0*inch) {
		ply = 9.0*inch;
		plx = ply/aspect;
	}

	myDomain = cpdf_createPlotDomain(.5*inch+(8.0*inch-plx)/2, (11*inch - ply)/2, plx, ply,
				xmin,xmax,ymin,ymax, LINEAR, LINEAR, 0);
	oldDomain = cpdf_setPlotDomain(pdf,myDomain);	/* save oldDomain for later restore */
	cpdf_fillDomainWithRGBcolor(myDomain, 1.0, 1.0, 0.9);	/* light yellow */
//	cpdf_drawMeshForDomain(myDomain);

	cpdf_setgray(pdf,0.0);
	/* X-Axis --------------------------------------------------------------------------------- */
	xAxis = cpdf_createAxis(0.0, plx, LINEAR, xmin, xmax);
	cpdf_attachAxisToDomain(xAxis, myDomain, 0.0, 0.0);
	cpdf_setAxisNumberFormat(xAxis, "%g", "Helvetica", 14.0);
	cpdf_setAxisLabel(xAxis, "R (m)", "Times-Roman", "MacRomanEncoding", 18.0);
	cpdf_drawAxis(xAxis);
	cpdf_freeAxis(xAxis);

	/* Y-Axis --------------------------------------------------------------------------------- */
	yAxis = cpdf_createAxis(90.0, ply, LINEAR, ymin, ymax);
	cpdf_attachAxisToDomain(yAxis, myDomain, 0.0, 0.0);
	cpdf_setAxisNumberFormat(yAxis, "%g", "Helvetica", 14.0);
	cpdf_setAxisLabel(yAxis, "Z (m)", "Times-Roman", "MacRomanEncoding", 18.0);
	cpdf_drawAxis(yAxis);
	cpdf_freeAxis(yAxis);

    /* title on top */
    cpdf_setgrayFill(pdf,0.0); /* Black */
    cpdf_beginText(pdf,0);
    cpdf_setFont(pdf,"Times-Roman", "MacRomanEncoding", 24.0);
    cpdf_textAligned(pdf,(xmin+xmax)/2., ymax+(ymax-ymin)*.035, 0.0, TEXTPOS_LM,
    	"Psi and |B| Contours");
    cpdf_endText(pdf);

	cpdf_gsave(pdf); /* to later undo clipping */
	cpdf_clipDomain(myDomain);

	/* first show roughly where the code thinks there is plasma */

	cpdf_setrgbcolorFill(pdf,1.0, 0.7, .7);	/* pink plasma */
	cpdf_setrgbcolorStroke(pdf,1.0, 0.7, .7);	/* pink plasma */
	cpdf_setlinewidth(pdf,.1);

	contour_IsPlasma(pdf,td);



	cpdf_setlinewidth(pdf,.5);
    showGeometry(pdf,td);

	/* Plot it */



	cpdf_setlinewidth(pdf,.1);
	cpdf_setrgbcolorStroke(pdf,0.0, 0.0, 1.0);	/* blue */

	contour_Psi(pdf,td);


	cpdf_setlinewidth(pdf,.1);
	cpdf_setrgbcolorStroke(pdf,0.0, 0.5, 0.5);
	cpdf_setdash(pdf,"[4 4] 0");

	contour_ModB(pdf,td);

	printf("\n");

	cpdf_grestore(pdf); /* remove clipping */


	cpdf_setPlotDomain(pdf,oldDomain);		/* restore previous plot domain */
	cpdf_freePlotDomain(myDomain);		/* deallocate the plot domain */
}
