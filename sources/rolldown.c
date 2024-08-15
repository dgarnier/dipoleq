/*
**
**   rolldown.c
**
**   given a starting point on an equilibrium, roll down the hill and
**   output the points as a function of psi
**
**   2/5/99 -- Darren T. Garnier, Columbia University
**
**
**
**
**
*/

#include <math.h>
#include "nrutil.h"
#include "tokamak.h"
#include "psigrid.h"
#include "interpolate.h"
#include "rolldown.h"
#include "contour.h"

#define PI          3.14159265358979323

/* a simple runge-kutta routine */
void rk4(double *y, double *dydx, int n, double x, double h, double *yout,
         void (*derivs)(double, double *, double *));

void RollDerivs(double , double *p, double *dpdz);

/* a simple runge-kutta routine */
void rk4(double *y, double *dydx, int n, double x, double h, double *yout,
         void (*derivs)(double, double *, double *))
{
	int i;
	double xh,hh,h6,*dym,*dyt,*yt;

	dym=dvector(1,n);
	dyt=dvector(1,n);
	yt=dvector(1,n);
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];
	(*derivs)(xh,yt,dyt);
	for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
	(*derivs)(xh,yt,dym);
	for (i=1;i<=n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	(*derivs)(x+h,yt,dyt);
	for (i=1;i<=n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	free_dvector(yt,1,n);
	free_dvector(dyt,1,n);
	free_dvector(dym,1,n);
}

static TOKAMAK *theTok;

void RollDerivs(double dummy, double *p, double *dpdz)
{
	double gpx, gpz, gp2;

	gpx = interpolate(theTok->PsiGrid, theTok->Plasma->GradPsiX, p[1], p[2]);
	gpz = interpolate(theTok->PsiGrid, theTok->Plasma->GradPsiZ, p[1], p[2]);
	gp2 = interpolate(theTok->PsiGrid, theTok->Plasma->GradPsi2, p[1], p[2]);

	dpdz[1] = gpx/gp2;
	dpdz[2] = gpz/gp2;
}

void RollDownHill(TOKAMAK *td, double r0, double z0,
                  int n, double *psi, double *rr, double *zz)
{
	double y_in[3], y_out[3], dydpsi[3];
	int i;

	theTok = td;

	*rr = y_in[1] = r0;
	*zz = y_in[2] = z0;
	*psi = GetPsi(td->PsiGrid, r0, z0);

	for (i=0; i<n-1 ; i++) {
		RollDerivs(psi[i],y_in,dydpsi);
		rk4(y_in, dydpsi, 2, psi[i], (psi[i+1]-psi[i]), y_out, RollDerivs);
		rr[i+1] = y_in[1] = y_out[1];
		zz[i+1] = y_in[2] = y_out[2];
		psi[i+1] = GetPsi(td->PsiGrid, rr[i+1], zz[i+1]);
	}
}

#define PT_FLOOR(x, xmin, dx)	(int) floor( ((x)-(xmin))/(dx) )

/***************************************************************
* invert
*
* This function interpolates/extrapolates an array A in order to
* return its value at the requested point, (x, z).
*
* We use a simple bilinear interpolation.
*
***************************************************************/
double        binterpolate(PSIGRID * pg, double **A, double x, double z)
{
    int           ix, iz;
    double        xm, zm;
    double        t, u;

    xm = DMAX(x, pg->Xmin);
    xm = DMIN(x, pg->Xmax);
    zm = DMIN(z, pg->Zmax);
    zm = DMAX(z, pg->Zmin);

    ix = PT_FLOOR(xm, pg->Xmin, pg->dx);
    iz = PT_FLOOR(zm, pg->Zmin, pg->dz);
    t = (xm - pg->X[ix]) / pg->dx;
    u = (zm - pg->Z[iz]) / pg->dz;

    return ((1.0 - t) * (1.0 - u) * A[ix][iz] + t * (1.0 - u) * A[ix + 1][iz] +
            t * u * A[ix + 1][iz + 1] + (1.0 - t) * u * A[ix][iz + 1]);
}

#define GETPSI(x,z) binterpolate(pg, pg->Psi, x, z)

#if 0
void RollOutToPsi(TOKAMAK *td, double *x, double *z, double psiNew)
{
    // this rolls out to major radius until it reaches new psi value

	int    i, ix, iz;
	int    nmax;
    double psiOld;
    double tx = *x;
    double tz = *z;
    PSIGRID *pg = td->PsiGrid;

    psiOld = GETPSI(tx, tz);

    ix = PT_FLOOR(*x, pg->Xmin, pg->dx);
    iz = PT_FLOOR(*z, pg->Zmin, pg->dz);

	nmax = pg->Nsize;

    // now lets set up the next contour
    for (tx=*x; tx<pg->Xmax, <nmax-1; i++) {
        if (pg->Psi[i+1] > psiNew) break;
    }
    if (i == nmax-1) nrerror("Couldn't determine midplane.\n");

        if ((pg->X[j]-rr[i-1])*(pg->X[j+1]-rr[i-1]) <= 0) ix = j;
        if ((pg->Z[j]-zz[i-1])*(pg->Z[j+1]-zz[i-1]) <= 0) iz = j;
    }
    while (pg->Psi[ix+1][iz] <= psi[i]) {
        if (ix >= nmax-1) nrerror("GetRandV: Couldn't follow flux contours.\n");
        ix++;
    }

}

#endif


static double gV, gRmax, gZrmax;
void	RmaxVStep(double x, double z, double y, int flag);

void	RmaxVStep(double x, double z, double dummy, int flag)
{
	static double 	Xlast, Zlast, IntLast;
	double 			Integrand, B, dS;

	B = sqrt(interpolate(theTok->PsiGrid,theTok->Plasma->GradPsi2,x,z));
	Integrand = 2*PI * x / B;

	switch (flag) {
	  case CONTOUR_START:
	      gRmax = x;
	      gZrmax = z;
	  	  gV = 0;
		  break;
	  case CONTOUR_TRACE:
	  case CONTOUR_STOP:
	      if (x >= gRmax) {
	      	gRmax = x;
	      	gZrmax = z;
	      }
	  	  dS = hypot((x-Xlast),(z-Zlast));
	  	  gV += dS*(Integrand+IntLast)/2.0;
		  break;
	}
	IntLast = Integrand;
	Xlast = x;
	Zlast = z;
}

void GetRandVfromPsi(TOKAMAK *td, int n, 
	double *psi, double *rr, double *zz, double *vv)
{
	PSIGRID *pg = td->PsiGrid;
	double r0, z0;
	theTok = td;  // static global for RmaxVStep because contour doesn't allow for user data

	// find the peak location via psi.. assume its set in Psi[0]
	contour(pg->X, pg->Z, pg->Psi, 1, pg->Nsize, 1, pg->Nsize,
		psi[0], CONTOUR_ONLY_CLOSED, CONTOUR_MIDPOINT, RmaxVStep);

	r0 = gRmax;
	z0 = gZrmax;

	// use the previous routine to get the rest
	GetRandV(td, r0, z0, n, psi, rr, zz, vv);
}

void GetRandV(TOKAMAK *td, double r0, double z0,
                  int n, double *psi, double *rr, double *zz, double *vv)
{
	PSIGRID *pg;

	int i,j,nmax,ix,iz;

	theTok = td;

	pg = td->PsiGrid;

	nmax = pg->Nsize;

	*rr = r0;
	*zz = z0;


	i=0;
	contour_point(pg->X, pg->Z, pg->Psi, 1, nmax, 1, nmax,
		   rr[i], zz[i], CONTOUR_MIDPOINT, RmaxVStep);
		vv[i] = gV;
		rr[i] = gRmax;
		zz[i] = gZrmax;
		psi[i] = GetPsi(td->PsiGrid, rr[i], zz[i]);


	for (i=1; i<n ; i++) {

		// now lets set up the next contour
		for (j=0; j<nmax; j++) {
			if ((pg->X[j]-rr[i-1])*(pg->X[j+1]-rr[i-1]) <= 0) ix = j;
			if ((pg->Z[j]-zz[i-1])*(pg->Z[j+1]-zz[i-1]) <= 0) iz = j;
		}
		while (pg->Psi[ix+1][iz] <= psi[i]) {
			if (ix >= nmax-1) nrerror("GetRandV: Couldn't follow flux contours.\n");
			ix++;
		}
		rr[i] =   pg->X[ix] +
		        ( psi[i]            - pg->Psi[ix][iz] ) *
		        ( pg->X[ix+1]       - pg->X[ix]       ) /
		        ( pg->Psi[ix+1][iz] - pg->Psi[ix][iz] ) ;
		zz[i] = pg->Z[iz];

		// now do the contour
		contour_point(pg->X, pg->Z, pg->Psi, 1, nmax, 1, nmax,
		   rr[i], zz[i], CONTOUR_MIDPOINT, RmaxVStep);
		vv[i] = gV;
		rr[i] = gRmax;
		zz[i] = gZrmax;
		psi[i] = GetPsi(td->PsiGrid, rr[i], zz[i]);
	}
}

#if 0
void GetRmidPsi(TOKAMAK *td, int npts, double *psi, double *xx, double *zz)
{
    // determine rr and zz of midplane, starting at center and working way from
    // inner to out using contouring technique
    // if you are a dipole, don't start from center but from LCFS */
    PSIGRID *pg;
    int i,j,nmax,ix,iz;
    double lastX, lastZ;

    pg = td->PsiGrid;
    theTok = td;

    nmax = pg->Nsize;

    lastX = pg->XMagAxis;
    lastZ = pg->ZMagAxis;

#ifndef DIPOLE
    *xx = lastX;
    *zz = lastZ;
    *psi = pg->PsiAxis;

    for(i=1; i<npts; i++)
#else
    for(i=0; i<npts; i++)
#endif
    {
        psi[i] = i*(pg->PsiLim-pg->PsiAxis)/(npts-1.d)+pg->PsiAxis;

        newX = 1;

        contour_point(pg->X, pg->Z, pg->Psi, 1, nmax, 1, nmax,
                  rr[i], zz[i], CONTOUR_MIDPOINT, RmaxVStep);
    vv[i] = gV;
    rr[i] = gRmax;
    zz[i] = gZrmax;
    psi[i] = GetPsi(td->PsiGrid, rr[i], zz[i]);


    for (i=1; i<n ; i++) {

        // now lets set up the next contour
        for (j=0; j<nmax; j++) {
            if ((pg->X[j]-rr[i-1])*(pg->X[j+1]-rr[i-1]) <= 0) ix = j;
            if ((pg->Z[j]-zz[i-1])*(pg->Z[j+1]-zz[i-1]) <= 0) iz = j;
        }
        while (pg->Psi[ix+1][iz] <= psi[i]) {
            if (ix >= nmax-1) nrerror("GetRandV: Couldn't follow flux contours.\n");
            ix++;
        }
        rr[i] =   pg->X[ix] +
        interpolate(	gpx = interpolate(pg, td->Plasma->GradPsiX, p[1], p[2]);
        ( psi[i]            - pg->Psi[ix][iz] ) *
        ( pg->X[ix+1]       - pg->X[ix]       ) /
        ( pg->Psi[ix+1][iz] - pg->Psi[ix][iz] ) ;
        zz[i] = pg->Z[iz];

        // now do the contour
        contour_point(pg->X, pg->Z, pg->Psi, 1, nmax, 1, nmax,
                      rr[i], zz[i], CONTOUR_MIDPOINT, RmaxVStep);
        vv[i] = gV;
        rr[i] = gRmax;
        zz[i] = gZrmax;
        psi[i] = GetPsi(td->PsiGrid, rr[i], zz[i]);
    }



}
#endif
