/*
** TokaMac v2.0
**
** GetFluxParameters.c
**
**
**
** File:		GetFluxParameters.c
** Date:		April 8, 1993
**
** Modifications:
**
**		May 3, 1993		Fixed typo in qStar formula
**		May 4, 1993		Bilinear interpolation in q(0) formula
**		May 20, 1993	Spline fit to flux surface for integrals
**		June 2, 1993	Added Grimm's mapper formula, Eq. 9
**		Jun 10, 1993	Removed possibility of atan2(0,0)
**		Sept 16, 1993	Added virial integrals (Lao, Nuc. Fus., 1985)
**
** Routine list:
**
**		NearPlasma(ip, ix, iz);			-- returns 1 if near plasma
**		Trace_Integrand(x, z, p, flag);	-- for tracing a psi surface
**		GetFluxParameters(td);			-- the main routine in this file
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdio.h>
#include <stdlib.h>
#include "VAX_Alloc.h"
#include <math.h>
#include "nrutil.h"
#include "contour.h"
#include "nrRomberg.h"
#include "nrSpline.h"
#include "interpolate.h"
#include "fpoly.h"
#include "psigrid.h"
#include "plasma.h"
#include "tokamak.h"
#include "GetFluxMoments.h"
#include "GetFluxParameters.h"
#include "multitask.h"
#include "CPlasmaModel.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

#define P_EDGE		0.0

#define SQUARE(x)	((x) * (x))
#define ISMIN4(x1,x2,x3,x4) (((x1) <= (x2)) && ((x1) <= (x2)) && ((x1) <= (x3)) && ((x1) <= (x4)))
#define BILIN(y1,y2,y3,y4) ((1.0-hx)*(1.0-hz)*(y1) + hx*(1.0-hz)*(y2) + hx*hz*(y3) + (1.0-hx)*hz*(y4))

#define COMPUTE_INT(pg, PsiX)   quick_Int(pg, PsiX, 0)
#define COMPUTE_INT1(pg, PsiX)  quick_Int(pg, PsiX, 1)

#define FIND_MIN(pg, PsiX)  find_Min(pg, PsiX)

//#define COMPUTE_INT(pg, PsiX)   compute_Int(pg, PsiX)
//#define COMPUTE_INT1(pg, PsiX)  compute_Int1(pg, PsiX)



/*
**	G L O B A L   V A R I A B L E S
*/

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif

EXTERN FILE  *LogFile;

EXTERN double gX0, gZ0;			/* the center of the flux surface used to define theta */
EXTERN double *gX, *gZ;			/* arrays containing the boundary coordinates */
EXTERN double *gTheta;
EXTERN double *gXsplines, *gZsplines;
EXTERN int    gCount;			/* the number of points representing the boundary */

double      **gIntegrand;		/* global integrand for flux integrals */
double       *gInt;				/* array containing integrand along flux surface */
double       *gIntsplines;

double        gXmin, gZmin;		/* for interpolations of the integrands */
double        gdx, gdz;
double		  gIntegral, gPathInt;

/*
**  Private Prototypes
*/

int     NearPlasma(int **ip, int ix, int iz);
void    Trace_Integrand(double x, double z, double p, int flag);
void    splint_dx(double xa[], double ya[], double y2a[], int n, double x, double *dydx);
double  Int_Theta(double theta);
double  Int_Theta1(double theta);
double  compute_Int(PSIGRID * pg, double PsiX);
double  compute_Int1(PSIGRID * pg, double PsiX);
void   	Fill_q_integrand(PSIGRID * pg, PLASMA * pl);
void	quick_Int_Step(double x, double z, double psi, int flag);
double  quick_Int(PSIGRID * pg, double PsiX, int average);

void	find_Min_Step(double x, double z, double psi, int flag);
double  find_Min(PSIGRID * pg, double PsiX);



void		quick_Int_Step(double x, double z, double dummy, int flag)
{
	static double 	Xlast, Zlast, IntLast;
	double 			Integrand, hx, hz, dS;
	int				ix, iz;

	ix = (int) floor((x - gXmin) / gdx);
	iz = (int) floor((z - gZmin) / gdz);
	hx = (x - gXmin) / gdx - ix;
	hz = (z - gZmin) / gdz - iz;

	Integrand = ((1.0 - hx) * (1.0 - hz) * gIntegrand[ix][iz]
				 + hx * (1.0 - hz) * gIntegrand[ix + 1][iz]
				 + hx * hz * gIntegrand[ix + 1][iz + 1]
				 + (1.0 - hx) * hz * gIntegrand[ix][iz + 1]);


	switch (flag) {
	  case CONTOUR_START:
	  	  gIntegral = gPathInt = 0;
		  break;
	  case CONTOUR_TRACE:
	  case CONTOUR_STOP:
	  	  dS = sqrt(SQUARE(x-Xlast)+SQUARE(z-Zlast));
	  	  gIntegral += dS*(Integrand+IntLast)/2.0;
	  	  gPathInt += dS;
		  break;
	}


	IntLast = Integrand;
	Xlast = x;
	Zlast = z;
}

static double gIntMin, gXimin, gZimin;

void		find_Min_Step(double x, double z, double dummy, int flag)
{
	double 			Integrand, hx, hz;
	int				ix, iz;

	ix = (int) floor((x - gXmin) / gdx);
	iz = (int) floor((z - gZmin) / gdz);
	hx = (x - gXmin) / gdx - ix;
	hz = (z - gZmin) / gdz - iz;

	Integrand = ((1.0 - hx) * (1.0 - hz) * gIntegrand[ix][iz]
				 + hx * (1.0 - hz) * gIntegrand[ix + 1][iz]
				 + hx * hz * gIntegrand[ix + 1][iz + 1]
				 + (1.0 - hx) * hz * gIntegrand[ix][iz + 1]);


	switch (flag) {
	  case CONTOUR_START:
	  	  gIntMin = Integrand;
	  	  		  break;
	  case CONTOUR_TRACE:
	  case CONTOUR_STOP:

		  break;
	}

	if (Integrand <= gIntMin) {
 		gIntMin = Integrand;
 		gXimin = x;
 		gZimin = z;
	}
}

double        quick_Int(PSIGRID * pg, double PsiX, int average)
{
	int           nmax;
	double      **Psi;
	double       *X, *Z;
	double        PsiBnd;

    MULTI;
	nmax = pg->Nsize;
	Psi = pg->Psi;
	X = pg->X;
	Z = pg->Z;
	PsiBnd = pg->PsiAxis + PsiX * pg->DelPsi;

	contour(X, Z, Psi, 0, nmax, 0, nmax, PsiBnd, CONTOUR_ONLY_CLOSED,
			CONTOUR_MIDPOINT, quick_Int_Step);

	if (average == 0) {
		return gIntegral;
	} else {
		return gIntegral/gPathInt;
	}
}


double        find_Min(PSIGRID * pg, double PsiX)
{
	int           nmax;
	double      **Psi;
	double       *X, *Z;
	double        PsiBnd;

    MULTI;
	nmax = pg->Nsize;
	Psi = pg->Psi;
	X = pg->X;
	Z = pg->Z;
	PsiBnd = pg->PsiAxis + PsiX * pg->DelPsi;

	contour(X, Z, Psi, 0, nmax, 0, nmax, PsiBnd, CONTOUR_ONLY_CLOSED,
			CONTOUR_MIDPOINT, find_Min_Step);

	return gIntMin;

}

/*
**	N E A R P L A S M A
*/
int           NearPlasma(int **ip, int ix, int iz)
{
	return (ip[ix][iz] || ip[ix + 1][iz] || ip[ix][iz + 1] || ip[ix + 1][iz + 1]
			|| ip[ix - 1][iz] || ip[ix][iz - 1] || ip[ix - 1][iz - 1]
			|| ip[ix - 1][iz + 1] || ip[ix + 1][iz - 1]);
}

/*
**	T R A C E _ I N T E G R A N D
*/
#define CHKSET(chk, op, set)   if (chk op set) set = chk

void          Trace_Integrand(double x, double z, double dummy, int flag)
{
	int           ix, iz;
	double        hx, hz, Integrand;
	static double        Xmin, Xmax, Zmin, Zmax;

	switch (flag) {
	  case CONTOUR_START:
	  	  Xmax = Xmin = x;
		  Zmax = Zmin = z;
		  gCount = 0;
		  break;
	  case CONTOUR_TRACE:
		  gCount++;
		  CHKSET(x,>,Xmax);
		  CHKSET(x,<,Xmin);
		  CHKSET(z,>,Zmax);
		  CHKSET(z,<,Zmin);
		  break;
	  case CONTOUR_STOP:
          gX0 = (Xmin+Xmax)/2;
          gZ0 = (Zmin+Zmax)/2;
		  break;
	}

	ix = (int) floor((x - gXmin) / gdx);
	iz = (int) floor((z - gZmin) / gdz);
	hx = (x - gXmin) / gdx - ix;
	hz = (z - gZmin) / gdz - iz;

	Integrand = ((1.0 - hx) * (1.0 - hz) * gIntegrand[ix][iz]
				 + hx * (1.0 - hz) * gIntegrand[ix + 1][iz]
				 + hx * hz * gIntegrand[ix + 1][iz + 1]
				 + (1.0 - hx) * hz * gIntegrand[ix][iz + 1]);

	gX[gCount] = x;
	gZ[gCount] = z;
	gInt[gCount] = Integrand;
}

/*
**	s p l i n t _ d x
**
**	This returns the first derivative of y(x) at x using the previously
** 	computed array of second derivatives y2a[].
*/
void          splint_dx(double xa[], double ya[], double y2a[], int n, double x, double *dydx)
{
	int           klo, khi, k;
	double        h, b, a;

	klo = 1;
	khi = n;
	while (khi - klo > 1) {
		k = (khi + klo) >> 1;
		if (xa[k] > x)
			khi = k;
		else
			klo = k;
	}
	h = xa[khi] - xa[klo];
	if (h == 0.0)
		nrerror("Bad xa input to routine splint_dx");
	a = (xa[khi] - x) / h;
	b = (x - xa[klo]) / h;
	*dydx = (ya[khi] - ya[klo]) / h - (3.0 * a * a - 1.0) * h * y2a[klo] / 6.0 +
		(3.0 * b * b - 1.0) * h * y2a[khi] / 6.0;
}

/*
** 	I N T _ T H E T A
**
**	The integrand I(theta) * ds/dtheta
*/
double        Int_Theta(double theta)
{
	double        Int, dsdth, dxdth, dzdth;

	splint(gTheta - 1, gInt - 1, gIntsplines - 1, gCount + 1, theta, &Int);

	splint_dx(gTheta - 1, gX - 1, gXsplines - 1, gCount + 1, theta, &dxdth);
	splint_dx(gTheta - 1, gZ - 1, gZsplines - 1, gCount + 1, theta, &dzdth);
	dsdth = sqrt(dxdth * dxdth + dzdth * dzdth);

	return Int * dsdth;
}

/*
** 	I N T _ T H E T A 1
**
**	The integrand I(theta)
*/
double        Int_Theta1(double theta)
{
	double        Integrand;

	splint(gTheta - 1, gInt - 1, gIntsplines - 1, gCount + 1, theta, &Integrand);

	return Integrand;
}

/*
**	c o m p u t e _ I n t
**
**	The function performs all of the required calculations to find
**	q(PsiX) and other flux integrals assuming that
**	gIntegrand has been properly set up.
*/
double        compute_Int(PSIGRID * pg, double PsiX)
{
	int           nmax;
	int           count;
	double      **Psi;
	double       *X, *Z;
	double        dXdt, dZdt, dIdt;
	double        PsiBnd;
	double        q;

    MULTI;

	nmax = pg->Nsize;
	Psi = pg->Psi;
	X = pg->X;
	Z = pg->Z;
	PsiBnd = pg->PsiAxis + PsiX * pg->DelPsi;

	contour(X, Z, Psi, 0, nmax, 0, nmax, PsiBnd, CONTOUR_ONLY_CLOSED,
			CONTOUR_NO_MIDPOINT, Trace_Count);

	count = gCount;				/* the number of points to represent this flux surface */

	/*	A L L O C A T E   M E M O R Y   F O R   S P L I N E S */
	gX = dvector(0, count);
	gZ = dvector(0, count);
	gTheta = dvector(0, count);
	gInt = dvector(0, count);
	gXsplines = dvector(0, count);
	gZsplines = dvector(0, count);
	gIntsplines = dvector(0, count);

	/*	T R A C E _ I N T E G R A N D */
	contour(X, Z, Psi, 0, nmax, 0, nmax, PsiBnd, CONTOUR_ONLY_CLOSED,
			CONTOUR_NO_MIDPOINT, Trace_Integrand);

// Trace_Integrand sets gX0 and gZ0
//	gX0 = pg->XMagAxis;
//	gZ0 = pg->ZMagAxis;
	FindTheta();

	/*  F I N D   S P L I N E S */
	dXdt = 0.5 * ((gX[count] - gX[count - 1]) / (gTheta[count] - gTheta[count - 1]) +
				  (gX[1] - gX[0]) / (gTheta[1] - gTheta[0]));
	spline(gTheta - 1, gX - 1, count + 1, dXdt, dXdt, gXsplines - 1);

	dZdt = 0.5 * ((gZ[count] - gZ[count - 1]) / (gTheta[count] - gTheta[count - 1]) +
				  (gZ[1] - gZ[0]) / (gTheta[1] - gTheta[0]));
	spline(gTheta - 1, gZ - 1, count + 1, dZdt, dZdt, gZsplines - 1);

	dIdt = 0.5 * ((gInt[count] - gInt[count - 1]) / (gTheta[count] - gTheta[count - 1]) +
				  (gInt[1] - gInt[0]) / (gTheta[1] - gTheta[0]));
	spline(gTheta - 1, gInt - 1, count + 1, dIdt, dIdt, gIntsplines - 1);

	/* I N T E G R A T E   O V E R   F L U X   S U R F A C E */
	q = qromb(Int_Theta, gTheta[0], gTheta[count]);

	/*  F R E E    M E M O R Y */
	free_dvector(gIntsplines, 0, count);
	free_dvector(gZsplines, 0, count);
	free_dvector(gXsplines, 0, count);
	free_dvector(gInt, 0, count);
	free_dvector(gTheta, 0, count);
	free_dvector(gZ, 0, count);
	free_dvector(gX, 0, count);

	return q;
}

/*
**	c o m p u t e _ I n t 1
**
**	The function performs all of the required calculations to find
**	q(PsiX) and other flux integrals assuming that
**	gIntegrand has been properly set up.
*/
double        compute_Int1(PSIGRID * pg, double PsiX)
{
	int           nmax;
	int           count;
	double      **Psi;
	double       *X, *Z;
	double        dXdt, dZdt, dIdt;
	double        PsiBnd;
	double        q;
MULTI;
	nmax = pg->Nsize;
	Psi = pg->Psi;
	X = pg->X;
	Z = pg->Z;
	PsiBnd = pg->PsiAxis + PsiX * pg->DelPsi;

	contour(X, Z, Psi, 0, nmax, 0, nmax, PsiBnd, CONTOUR_ONLY_CLOSED,
			CONTOUR_NO_MIDPOINT, Trace_Count);

	count = gCount;				/* the number of points to represent this flux surface */

	/*	A L L O C A T E   M E M O R Y   F O R   S P L I N E S */
	gX = dvector(0, count);
	gZ = dvector(0, count);
	gTheta = dvector(0, count);
	gInt = dvector(0, count);
	gXsplines = dvector(0, count);
	gZsplines = dvector(0, count);
	gIntsplines = dvector(0, count);

	/*	T R A C E _ I N T E G R A N D */
	contour(X, Z, Psi, 0, nmax, 0, nmax, PsiBnd, CONTOUR_ONLY_CLOSED,
			CONTOUR_NO_MIDPOINT, Trace_Integrand);

//	gX0 = pg->XMagAxis;
//	gZ0 = pg->ZMagAxis;
	FindTheta();

	/*  F I N D   S P L I N E S */
	dXdt = 0.5 * ((gX[count] - gX[count - 1]) / (gTheta[count] - gTheta[count - 1]) +
				  (gX[1] - gX[0]) / (gTheta[1] - gTheta[0]));
	spline(gTheta - 1, gX - 1, count + 1, dXdt, dXdt, gXsplines - 1);

	dZdt = 0.5 * ((gZ[count] - gZ[count - 1]) / (gTheta[count] - gTheta[count - 1]) +
				  (gZ[1] - gZ[0]) / (gTheta[1] - gTheta[0]));
	spline(gTheta - 1, gZ - 1, count + 1, dZdt, dZdt, gZsplines - 1);

	dIdt = 0.5 * ((gInt[count] - gInt[count - 1]) / (gTheta[count] - gTheta[count - 1]) +
				  (gInt[1] - gInt[0]) / (gTheta[1] - gTheta[0]));
	spline(gTheta - 1, gInt - 1, count + 1, dIdt, dIdt, gIntsplines - 1);

	/* I N T E G R A T E   O V E R   F L U X   S U R F A C E */
	q = qromb(Int_Theta1, gTheta[0], gTheta[count]);

	/*  F R E E    M E M O R Y */
	free_dvector(gIntsplines, 0, count);
	free_dvector(gZsplines, 0, count);
	free_dvector(gXsplines, 0, count);
	free_dvector(gInt, 0, count);
	free_dvector(gTheta, 0, count);
	free_dvector(gZ, 0, count);
	free_dvector(gX, 0, count);

	return q;
}

/*
**	F I L L _ Q _ I N T E G R A N D
**
**
*/
void          Fill_q_integrand(PSIGRID * pg, PLASMA * pl)
{
	int           ix, iz, nmax;
	double        r, theta;
	double        xa, za;
	double        dx, dz;
	double       *X, *Z, **dPsiX, **dPsiZ;

	xa = pg->XMagAxis;
	za = pg->ZMagAxis;

	nmax = pg->Nsize;
	dPsiX = pl->GradPsiX;
	dPsiZ = pl->GradPsiZ;
	X = pg->X;
	Z = pg->Z;

	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (NearPlasma(pg->IsPlasma, ix, iz)) {
				dx = X[ix] - xa;
				dz = Z[iz] - za;
				r = sqrt(dx * dx + dz * dz);
				if (r == 0.0)	/* then we're at the magnetic axis & atan2 is not defined */
					gIntegrand[ix][iz] = 0.5 * (gIntegrand[ix - 1][iz] + gIntegrand[ix][iz - 1]);
				else {
					theta = atan2(dz, dx);
					gIntegrand[ix][iz] =
						pl->Bt[ix][iz] * r / (dPsiX[ix][iz] * cos(theta) + dPsiZ[ix][iz] * sin(theta));
				}
			} else
				gIntegrand[ix][iz] = 0.0;
}


/*
**  ComputeFluxFunctions
*/
void		ComputeFluxFunctions(TOKAMAK *td)
{
	PLASMA       *pl;
	int           i, npts;
	double        PsiXmax, PsiX, Psi, DelPsi, P, G, Pp, G2p;
	double		  *PV, *GV, *PpV, *G2V, *PsiV, *PsiXV;

	/* F L U X   P R O F I L E S */
	pl = td->Plasma;
	npts = pl->NumPsiPts;
	PsiXmax = pl->PsiXmax;
	DelPsi = pl->PsiLim - pl->PsiAxis;

	/*  A L L O C A T E   M E M O R Y */
	PsiV  = pl->Psi_pr  = dvector(0, npts-1);
	PsiXV = pl->PsiX_pr = dvector(0, npts-1);
	PV    = pl->P_pr    = dvector(0, npts-1);
	GV    = pl->G_pr    = dvector(0, npts-1);
	PpV   = pl->Pp_pr   = dvector(0, npts-1);
	G2V   = pl->G2p_pr  = dvector(0, npts-1);

	for (i = 0; i < npts; i++) {
		PsiX = PsiXV[i] = i * PsiXmax / (npts - 1);
		Psi = PsiV[i] = pl->PsiAxis + PsiX * DelPsi;
		Pp = P = 0.0;
		switch (pl->ModelType) {
			case Plasma_Std:
				P = -DelPsi * pl->Pp[1] * pow(1.0 - PsiX, pl->StndP) / pl->StndP;
				G = 1.0 - DelPsi * pl->G2p[1] * pow(1.0 - PsiX, pl->StndG) / pl->StndG;
				Pp = pl->Pp[1] * pow(1.0 - PsiX, pl->StndP - 1.0);
				G2p = pl->G2p[1] * pow(1.0 - PsiX, pl->StndG - 1.0);
			break;

			case  Plasma_IsoNoFlow:
				P = fpoly_int(pl->Pp, PsiX, pl->PpTerms, DelPsi, P_EDGE);
				G = fpoly_int(pl->G2p, PsiX, pl->G2pTerms, DelPsi, 1.0);
				Pp = fpoly(pl->Pp, PsiX, pl->PpTerms);
				G2p = fpoly(pl->G2p, PsiX, pl->G2pTerms);
			break;

			default :
			    if (pl->Model) {
			    	P   = pl->Model->P(Psi);
			    	Pp  = pl->Model->Pp(Psi);
			    	G   = pl->Model->G2(Psi);
			    	G2p = pl->Model->G2p(Psi);
			    }

		}
		GV[i]  = G  = sqrt(G);
		PV[i]  = P  = P / MU0;
		PpV[i] = Pp = Pp / MU0;
		G2V[i] = G2p;
	}

}


/*
**	G E T F L U X P A R A M E T E R S
**
**	This routine performs contour integrals around the flux functions
** 	in order to compute surface quantities such as the safety factor, q.
**
**	Our main reference is Chapter 6 of Jeff Freidberg's textbook
**	Ideal Magnetohydrodynamics.
**
**	In the following, all are "flux functions" and variables ending
**	with "p" (or "prime") denote derivatives of flux functions
**	with respect to flux (Psi).
**
**	We compute:
**		 	q(x)		safety factor
**			S(x)		global shear = 2(Vol/Volp)(qp/q)
**			Volp(x)		flux derivative of volume
**			Well(x)		Magnetic well
**			...and a few others.
**
**	NOTE:
**
**	An important part of this routine is the calculation of q(psi).
**	This has taken some work.  We first find a arrays representing the
**	grid locations of a flux surface (x,z,theta), where theta is an angle
**	parameter defined from a reference.  We also compute the integrand (I)
**	at each of the grid points.  Then we spline fit x(theta), z(theta), and
**	I(theta).  The flux integral is then calculated by finding...
**
**		Int(ds I(s)) = Int(dth I(th)*sqrt((dx/dth)^2 + (dz/dth)^2))
**
**	where th = theta, and ds = dth*sqrt((dx/dth)^2 + (dz/dth)^2)).
**
*/
void          GetFluxParameters(TOKAMAK * td)
{
	PLASMA       *pl;
	PSIGRID      *pg;
	int           npts;			/* number of psi profile points */
	int           nmax, ix, iz, ixa, iza;
	int           i;
	double      **Psi;
	double       *X, *Z, dx, dz, hx, hz;
	double        Bt[4], dPsiX2[4], dPsiZ2[4];
	double        t1, t2, t3;
	double        PsiX;			/* normalized Psi */
	double        DelPsi;
	double       *Pp;			/* �<P>/�Psi */
	double        R0, Wnorm;

	ComputeFluxFunctions(td);

	pg = td->PsiGrid;
	pl = td->Plasma;
	nmax = pg->Nsize;
	Psi = pg->Psi;
	DelPsi = pg->DelPsi;
	X = pg->X;
	Z = pg->Z;

	npts = pl->NumPsiPts;

	gXmin = pg->Xmin;
	gZmin = pg->Zmin;
	gdx = dx = pg->dx;
	gdz = dz = pg->dz;


#ifdef USE_ROLL_DOWNHILL_POINTS

        xpsi = dvector(0,npts-1);
        zpsi = dvector(0,npts-1);
        vpsi = dvector(0,npts-1);
        ppsi = dvector(0,npts-1);

        for (i = 0; i < npts; i++) {
            ppsi[i] = i * pl->PsiXMax / (npts - 1.0);
            GetRandV(td, pl->XMagAxis, pl->ZMagAxis, npts, ppsi, xpsi, zpsi, vpsi);
#endif

	/* A L L O C A T E   S O M E   M E M O R Y */

	gIntegrand = dmatrix(0, nmax, 0, nmax);


	/* S A F E T Y   F A C T O R   O N   A X I S */

	/* find closest grid point to magnetic axis */
	/* use bilinear interpolation near axis */
	ixa = (int) floor((pl->XMagAxis - pg->Xmin) / dx);
	iza = (int) floor((pl->ZMagAxis - pg->Zmin) / dz);
	hx = (pl->XMagAxis - X[ixa]) / dx;
	hz = (pl->ZMagAxis - Z[iza]) / dz;
	/* toroidal field at axis */
	Bt[0] = pl->Bt[ixa][iza];
	Bt[1] = pl->Bt[ixa + 1][iza];
	Bt[2] = pl->Bt[ixa + 1][iza + 1];
	Bt[3] = pl->Bt[ixa][iza + 1];
	t1 = BILIN(Bt[0], Bt[1], Bt[2], Bt[3]);
	/* d2 Psi/ dx2 and d2 Psi/ dz2 */
	dPsiX2[0] = (Psi[ixa + 1][iza] - 2.0 * Psi[ixa][iza] + Psi[ixa - 1][iza]);
	dPsiX2[1] = (Psi[ixa + 2][iza] - 2.0 * Psi[ixa + 1][iza] + Psi[ixa][iza]);
	dPsiX2[2] = (Psi[ixa + 2][iza + 1] - 2.0 * Psi[ixa + 1][iza + 1] + Psi[ixa][iza + 1]);
	dPsiX2[3] = (Psi[ixa + 1][iza + 1] - 2.0 * Psi[ixa][iza + 1] + Psi[ixa - 1][iza + 1]);
	t2 = BILIN(dPsiX2[0], dPsiX2[1], dPsiX2[2], dPsiX2[3]) / dx / dx;
	dPsiZ2[0] = (Psi[ixa][iza + 1] - 2.0 * Psi[ixa][iza] + Psi[ixa][iza - 1]);
	dPsiZ2[1] = (Psi[ixa + 1][iza + 1] - 2.0 * Psi[ixa + 1][iza] + Psi[ixa + 1][iza - 1]);
	dPsiZ2[2] = (Psi[ixa + 1][iza + 2] - 2.0 * Psi[ixa + 1][iza + 1] + Psi[ixa + 1][iza]);
	dPsiZ2[3] = (Psi[ixa][iza + 2] - 2.0 * Psi[ixa][iza + 1] + Psi[ixa][iza]);
	t3 = BILIN(dPsiZ2[0], dPsiZ2[1], dPsiZ2[2], dPsiZ2[3]) / dz / dz;
	if (t2 && t3)
		pl->q0 = TWOPI * t1 / sqrt(t2 * t3);
	else
		pl->q0 = 1.0e3;			/* a big number */
	pl->qCircular = 5.0e6 * DSQR(pl->HalfWidth) *
		(pl->B0R0 / pl->RCenter / pl->RCenter) / pl->Ip;
	pl->qStar = pl->qCircular * (1.0 + DSQR(pl->Elongation)) / 2.0;


	/* S A F E T Y   F A C T O R   P R O F I L E */

	/*
	**	From Grimm, et al, in "Methods of Computational Physics"
	**	(Killeen, ed.) Vol 16., p. 253, Academic Press, New York, 1976.
	**
	**	Eq. 9 used with the definition of q, gives the forumla
	**
	**		q = Int( Bt r dtheta / ( dPsi/dx cos(theta) + dPsi/dz sin(theta)) )
	**
	**	where the integral is around the flux surface from 0 to 2pi.
	**
	**	We use the magnetic axis as the reference point for the evaluation of
	**	(r, theta).  Note, we use a right-handed definition of theta.
	**
	*/

	pl->q_pr = dvector(0, npts - 1);
	pl->q_pr[0] = pl->q0;

	Fill_q_integrand(pg, pl);

	for (i = 1; i < npts; i++) {
		PsiX = i * pl->PsiXmax / (npts - 1.0);
		pl->q_pr[i] = COMPUTE_INT1(pg, PsiX);
	}

	/* d V O L U M E  / d P S I */

	pl->Volp_pr = dvector(0, npts - 1);
#ifndef DIPOLE
	pl->Volp_pr[0] = 2.0 * TWOPI / pg->Current[ixa][iza];
#endif
	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (NearPlasma(pg->IsPlasma, ix, iz)) {
				t1 = sqrt(pl->GradPsi2[ix][iz]);
				gIntegrand[ix][iz] = TWOPI * X[ix] / t1;
			}
#ifndef DIPOLE
	for (i = 1; i < npts; i++) {
#else
    for (i = 0; i < npts; i++) {
#endif
		PsiX = i * pl->PsiXmax / (npts - 1.0);
		pl->Volp_pr[i] = COMPUTE_INT(pg, PsiX);
	}

	/* V O L U M E    P R O F I L E */

	pl->Vol_pr = dvector(0, npts - 1);
	pl->Vol_pr[0] = 0.0;

	t1 = DelPsi * pl->PsiXmax / (npts - 1.0);	/* del Psi between profile points */
	for (i = 1; i < npts; i++) {
		pl->Vol_pr[i] = pl->Vol_pr[i - 1] + 0.5 * t1 * (pl->Volp_pr[i - 1] + pl->Volp_pr[i]);
	}

	/* G L O B A L   S H E A R */

	pl->S_pr = dvector(0, npts - 1);
	pl->S_pr[0] = 0.0;

	t1 = pl->PsiXmax / (npts - 1.0);	/* del PsiX between profile points */
	for (i = 1; i < npts - 1; i++) {
		pl->S_pr[i] = (pl->Vol_pr[i] / pl->Vol_pr[i])
			* (pl->q_pr[i + 1] - pl->q_pr[i - 1]) / t1 / pl->q_pr[i];
	}
	/* edge point */
	pl->S_pr[npts - 1] = 2.0 * (pl->Vol_pr[npts - 1] / pl->Vol_pr[npts - 1])
		* (pl->q_pr[npts - 1] - pl->q_pr[npts - 2]) / t1 / pl->q_pr[npts - 1];

	/* F L U X   S U R F A C E   A V G   O F   B 2 */

	pl->B2_pr = dvector(0, npts - 1);
	pl->B2_pr[0] = pl->B2[ixa][iza];	/* toroidal field at axis */

	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (NearPlasma(pg->IsPlasma, ix, iz)) {
				t1 = sqrt(pl->GradPsi2[ix][iz]);
				gIntegrand[ix][iz] = pl->B2[ix][iz] * TWOPI * X[ix] / t1;
			}
#ifndef DIPOLE
	for (i = 1; i < npts; i++) {
#else
    for (i = 0; i < npts; i++) {
#endif
		PsiX = i * pl->PsiXmax / (npts - 1.0);
		pl->B2_pr[i] = COMPUTE_INT(pg, PsiX) / pl->Volp_pr[i];
	}

	/* F L U X   T U B E   A V G   O F   B E T A */


	pl->Beta_pr = dvector(0, npts - 1);

	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (NearPlasma(pg->IsPlasma, ix, iz)) {
				t1 = sqrt(pl->GradPsi2[ix][iz]);
				gIntegrand[ix][iz] = 2*pl->Piso[ix][iz]/pl->B2[ix][iz] * TWOPI * X[ix] / t1;
			}

#ifndef DIPOLE
	pl->Beta_pr[0] = 2*pl->Piso[ixa][iza]/pl->B2[ixa][iza];	/* beta on axis */

	for (i = 1; i < npts; i++) {
#else
    for (i = 0; i < npts; i++) {
#endif
		PsiX = i * pl->PsiXmax / (npts - 1.0);
		pl->Beta_pr[i] = COMPUTE_INT(pg, PsiX) / pl->Volp_pr[i];

	}

    /* F I N D   M A X   B E T A   P O S I T I O N */

    pl->BetaMax_pr  = dvector(0, npts - 1);
	pl->XBetaMax_pr = dvector(0, npts - 1);
	pl->ZBetaMax_pr = dvector(0, npts - 1);
	pl->BBetaMax_pr = dvector(0, npts - 1);

	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (NearPlasma(pg->IsPlasma, ix, iz)) {
				gIntegrand[ix][iz] = pl->B2[ix][iz];
			}

#ifndef DIPOLE
    pl->BetaMax_pr[0] = pl->Beta_pr[0];
    pl->XBetaMax_pr[0] = pg->X[ixa];
    pl->ZBetaMax_pr[0] = pg->Z[iza];
    pl->BBetaMax_pr[0] = sqrt(pg->B2[ixa][iza]);

	for (i = 1; i < npts; i++) {
#else
    for (i = 0; i < npts; i++) {
#endif
		PsiX = i * pl->PsiXmax / (npts - 1.0);
		FIND_MIN(pg, PsiX);
		pl->BetaMax_pr[i]  = 2*PlasmaP(pl, pg->PsiAxis + PsiX * pg->DelPsi)/gIntMin;
		pl->XBetaMax_pr[i] = gXimin;
		pl->ZBetaMax_pr[i] = gZimin;
		pl->BBetaMax_pr[i] = sqrt(gIntMin);
	}

    /* F I N D   M A X   F I E L D   P O S I T I O N */

        pl->BMax_pr = dvector(0, npts - 1);
	pl->XBMax_pr = dvector(0, npts - 1);
	pl->ZBMax_pr = dvector(0, npts - 1);

	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (NearPlasma(pg->IsPlasma, ix, iz)) {
				gIntegrand[ix][iz] = 1.0/pl->B2[ix][iz];
			}

#ifndef DIPOLE
    pl->XBMax_pr[0] = pg->X[ixa];
    pl->ZBMax_pr[0] = pg->Z[iza];
    pl->BMax_pr[0] = sqrt(pg->B2[ixa][iza]);

	for (i = 1; i < npts; i++) {
#else
    for (i = 0; i < npts; i++) {
#endif
		PsiX = i * pl->PsiXmax / (npts - 1.0);
		FIND_MIN(pg, PsiX);
		pl->XBMax_pr[i] = gXimin;
		pl->ZBMax_pr[i] = gZimin;
		pl->BMax_pr[i] = sqrt(1.0/gIntMin);
	}



	/* F L U X   S U R F A C E   A V G   O F   C U R R E N T */

	pl->J_pr = dvector(0, npts - 1);
	pl->J_pr[0] = pg->Current[ixa][iza] / MU0;	/* current density at axis */

	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (NearPlasma(pg->IsPlasma, ix, iz)) {
				t1 = sqrt(pl->GradPsi2[ix][iz]);
				gIntegrand[ix][iz] = pg->Current[ix][iz] * TWOPI * X[ix] / t1;
			}
#ifndef DIPOLE
	for (i = 1; i < npts; i++) {
#else
    for (i = 0; i < npts; i++) {
#endif
		PsiX = i * pl->PsiXmax / (npts - 1.0);
		pl->J_pr[i] = COMPUTE_INT(pg, PsiX) / pl->Volp_pr[i] / MU0;
	}

	/* M A G N E T I C   W E L L */

	Pp = dvector(0, npts - 1);
	Pp[0] = pl->Piso[ixa][iza];	/* pressure at magnetic axis */

	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (NearPlasma(pg->IsPlasma, ix, iz)) {
				t1 = sqrt(pl->GradPsi2[ix][iz]);
				gIntegrand[ix][iz] = pl->Piso[ix][iz] * TWOPI * X[ix] / t1;
			}
#ifndef DIPOLE
	for (i = 1; i < npts; i++) {
#else
    for (i = 0; i < npts; i++) {
#endif
		PsiX = i * pl->PsiXmax / (npts - 1.0);
		Pp[i] = COMPUTE_INT(pg, PsiX) / pl->Volp_pr[i];
	}

	pl->Well_pr = dvector(0, npts - 1);
	for (i = 0; i < npts; i++)
		pl->Well_pr[i] = 0.0;

	t1 = pl->PsiXmax / (npts - 1.0);	/* del PsiX between profile points */
	for (i = 1; i < npts; i++)
		switch (pl->ModelType) {
		  case Plasma_Std:
		  case Plasma_IsoNoFlow:
		  case Plasma_DipoleStd:
		  case Plasma_DipoleIntStable:
		  case Plasma_DipoleStablePsiN:
			  t2 = MU0 * (Pp[i + 1] - Pp[i - 1]) / t1;	/* Pressure gradient */
			  pl->Well_pr[i] = (pl->Vol_pr[i] / pl->Vol_pr[i])
				  * ((pl->B2_pr[i + 1] - pl->B2_pr[i - 1]) / 2.0 / t1 + t2) / pl->B2_pr[i];
			  break;
		  case Plasma_IsoFlow:
			  break;
		  case Plasma_AnisoNoFlow:
			  break;
		  case Plasma_AnisoFlow:
			  break;
		}
	/* edge point */
	switch (pl->ModelType) {
	  case Plasma_Std:
	  case Plasma_IsoNoFlow:
	  case Plasma_DipoleStd:
	  case Plasma_DipoleIntStable:
	  case Plasma_DipoleStablePsiN:
		  t2 = 2.0 * MU0 * (Pp[npts - 1] - Pp[npts - 2]) / t1;	/* Pressure gradient */
		  pl->Well_pr[npts - 1] = (pl->Vol_pr[npts - 1] / pl->Vol_pr[npts - 1])
			  * ((pl->B2_pr[npts - 1] - pl->B2_pr[npts - 2]) / 4.0 / t1 + t2) / pl->B2_pr[npts - 1];
		  break;
	  case Plasma_IsoFlow:
		  break;
	  case Plasma_AnisoNoFlow:
		  break;
	  case Plasma_AnisoFlow:
		  break;
	}

	/*  F R E E   S O M E   M E M O R Y  */
	free_dvector(Pp, 0, npts - 1);

	/* V I R I A L   I N T E G R A L S   */
	R0 = pl->RCenter;
	Wnorm = 0.25 * MU0 * pl->RStar * pl->Ip * pl->Ip;

	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (NearPlasma(pg->IsPlasma, ix, iz)) {
				t1 = sqrt(pl->GradPsi2[ix][iz]);
				t2 = (X[ix] - R0) * pl->GradPsiX[ix][iz] + Z[iz] * pl->GradPsiZ[ix][iz];
				gIntegrand[ix][iz] = t1 * t2 / (TWOPI * X[ix] * 2.0 * MU0);
			}
	pl->S1_vr = COMPUTE_INT(pg, pl->PsiXmax) / Wnorm;

	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (NearPlasma(pg->IsPlasma, ix, iz)) {
				t1 = sqrt(pl->GradPsi2[ix][iz]);
				t2 = R0 * pl->GradPsiX[ix][iz];
				gIntegrand[ix][iz] = t1 * t2 / (TWOPI * X[ix] * 2.0 * MU0);
			}
	pl->S2_vr = COMPUTE_INT(pg, pl->PsiXmax) / Wnorm;

	for (ix = 1; ix < nmax; ix++)
		for (iz = 1; iz < nmax; iz++)
			if (NearPlasma(pg->IsPlasma, ix, iz)) {
				t1 = sqrt(pl->GradPsi2[ix][iz]);
				t2 = Z[iz] * pl->GradPsiZ[ix][iz];
				gIntegrand[ix][iz] = t1 * t2 / (TWOPI * X[ix] * 2.0 * MU0);
			}
	pl->S3_vr = COMPUTE_INT(pg, pl->PsiXmax) / Wnorm;

	/*  F R E E   S O M E   M E M O R Y  */
	free_dmatrix(gIntegrand, 0, nmax, 0, nmax);
}
