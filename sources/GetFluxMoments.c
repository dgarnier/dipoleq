/*
** TokaMac v2.0
**
** GetFluxMoments.c
**
**
**
** File:		GetFluxMoments.c
** Date:		April 10, 1993
**
** Modifications:
**
**		May 3, 1993				Allow iteration for (X0,Z0) in FindCenterMoment
**		May 11, 1993			Added spline interpolation for flux boundary
**		May 11, 1993			Added GetFluxMoments with new argument list
**		May 24, 1993			Whoops. Repaired symmetric zm moments.
**
** Routine list:
**
**		GetFluxMoments			Find average position of outer flux surface.
**								Takes arguments for arrays for Xm and Zm and the
**								number of moments between [0..m].
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdio.h>
#include <stdlib.h>
#include "VAX_Alloc.h"
#include <math.h>
#include "nrutil.h"
#include "contour.h"
#include "interpolate.h"
#include "nrRomberg.h"
#include "nrSpline.h"
#include "psigrid.h"
#include "GetFluxMoments.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

#define BNDZERONUM	20
#define BNDDIFFTHR	0.0001

void          Trace_Boundary(double x, double z, double p, int flag);
double        CosNorm(double theta);
double        SinNorm(double theta);
double        XCosm(double theta);
double        ZSinm(double theta);


/*
**	G L O B A L   V A R I A B L E S
*/

double        gX0, gZ0;			/* the center of the flux surface used to define theta */

double       *gX, *gZ;			/* arrays containing the boundary coordinates */
double       *gTheta;
double       *gXsplines, *gZsplines;
int           gCount;			/* the number of points representing the boundary */
int	      gLen;
int           gm;				/* the moment being integrated */

/*
**	Trace_Count
**
**	This routine is used to trace the number of points used to describe
**	the flux surface.
**
**	The count begins at 0, and ends at N. (N + 1 entries.)
*/

void          Trace_Count(double x, double z, double dummy, int flag)
{
#pragma unused(x,z,dummy)
    static int count = 0;
    switch (flag) {
        case CONTOUR_START:
            count = 0;
            break;
        case CONTOUR_TRACE:
            count++;
            break;
        case CONTOUR_STOP:
            gCount = count;
            return;
    }
 //   printf("I = % 5d, X = %f, Z = %f\n",count,x,z);
}

/*
**	Trace_Boundary
*/
void          Trace_Boundary(double x, double z, double dummy, int flag)
{
#pragma unused(dummy)
    static int count = 0;

    switch (flag) {
        case CONTOUR_START:
            count = 0;
            break;
        case CONTOUR_TRACE:
            count++;
            break;
        case CONTOUR_STOP:
            return;
    }

    if (count > gLen || count < 0) {
        return;
    }
//    printf("I = % 5d, X = %f, Z = %f\n",count,x,z);
    gX[count] = x;
    gZ[count] = z;
}

/*
** CosNorm
**
*/
double        CosNorm(double theta)
{
	double        t;
	if (gm == 0)
		t = 1.0;
	else
		t = cos(gm * theta);
	return t * t;
}

/*
** SinNorm
**
*/
double        SinNorm(double theta)
{
	double        t;
	if (gm == 0)
		t = 1.0;
	else
		t = sin(gm * theta);
	return t * t;
}

/*
** XCosm
**
*/
double        XCosm(double theta)
{
	double        t, xx;

	splint(gTheta - 1, gX - 1, gXsplines - 1, gCount + 1, theta, &xx);

	if (gm == 0)
		t = xx;
	else
		t = xx * cos(gm * theta);
	return t;
}

/*
** ZSinm
**
*/
double        ZSinm(double theta)
{
	double        t, zz;

	splint(gTheta - 1, gZ - 1, gZsplines - 1, gCount + 1, theta, &zz);

	if (gm == 0)
		t = zz;
	else
		t = zz * sin(gm * theta);
	return t;
}

/*
**	FindTheta
**
**	Assumes that countour(..) traces in the positive theta direction.
**  This works for PsiAxis < PsiLim.
*/
void          FindTheta(void)
{
	int           i;
	double        trange = 0.0;

	for (i = 0; i <= gCount; i++) {
		gTheta[i] = atan2(gZ[i] - gZ0, gX[i] - gX0);
		if (gTheta[i] < 0.0)
			gTheta[i] += TWOPI;
		gTheta[i] += trange;	/* insures that theta always increases */
		if ((i > 0) && ((gTheta[i - 1] - gTheta[i]) > PI))
			trange += TWOPI;
	}
}

/*
**	GetFluxContour
**
**   THIS ROUTINE HAS A BUG... IT CAN GET MORE THAN ONE CONTOUR!
*/

void	GetFluxContour(PSIGRID * pg, double PsiX, double **x, double **z, int *len)
{
        int nmax;
        double **Psi, PsiBnd, *X, *Z;

    	nmax = pg->Nsize;
	Psi = pg->Psi;
	X = pg->X;
	Z = pg->Z;

	PsiBnd = pg->PsiAxis + PsiX * pg->DelPsi;

        gCount=0;

	contour(X, Z, Psi, 0, nmax, 0, nmax, PsiBnd, CONTOUR_ONLY_CLOSED,
			CONTOUR_MIDPOINT, Trace_Count);

   //     printf("Count on contour is %ld\n",gCount);

        if (gCount == 0) {
            nrinfo("Couldn't get boundary flux contour in GetFluxContour.\n");
 //           nrinfo("Boundary flux value = %f\n",PsiBnd);
            *x = NULL;
            *z = NULL;
            return;
        }


	*len = gCount;		/* the number of points to represent this flux surface */

        gLen = *len;

	*x = gX = dvector(0, *len);
	*z = gZ = dvector(0, *len);

	contour(X, Z, Psi, 0, nmax, 0, nmax, PsiBnd, CONTOUR_ONLY_CLOSED,
			CONTOUR_MIDPOINT, Trace_Boundary);
}

void          GetFluxMoments(PSIGRID * pg, double PsiX, double *Xm, double *Zm, int m)
{
	int           nmax;
	double      **Psi;
	double       *X, *Z, Zsym;
        double	     *bX, *bZ;
        int	     blen;
	double        PsiBnd;
	int           count, i;
	double        dXdt, dZdt, del;
	double        xnorm, znorm;

	nmax = pg->Nsize;
	Psi = pg->Psi;
	X = pg->X;
	Z = pg->Z;
	Zsym = Z[nmax / 2];

        GetFluxContour(pg, PsiX, &bX, &bZ, &blen);

        gTheta = dvector(0, blen);
	gXsplines = dvector(0, blen);
	gZsplines = dvector(0, blen);

	/* iterate flux integral in order to find center, X0, Z0 */
	gm = 0;
	Xm[0] = pg->XMagAxis;
	Zm[0] = pg->ZMagAxis;
	if (pg->Symmetric)
		Zm[0] = Zsym;
	for (i = 1; i < BNDZERONUM; i++) {
		gX0 = Xm[0];
		gZ0 = Zm[0];

		FindTheta();

		/* Fit spline to x[theta] and z[theta] */
		dXdt = 0.5 * ((gX[blen] - gX[blen - 1]) / (gTheta[blen] - gTheta[blen - 1]) +
					  (gX[1] - gX[0]) / (gTheta[1] - gTheta[0]));
		spline(gTheta - 1, gX - 1, blen + 1, dXdt, dXdt, gXsplines - 1);

		dZdt = 0.5 * ((gZ[blen] - gZ[blen - 1]) / (gTheta[blen] - gTheta[blen - 1]) +
					  (gZ[1] - gZ[0]) / (gTheta[1] - gTheta[0]));
		spline(gTheta - 1, gZ - 1, blen + 1, dZdt, dZdt, gZsplines - 1);

		xnorm = qromb(CosNorm, gTheta[0], gTheta[blen]);
		Xm[0] = qromb(XCosm, gTheta[0], gTheta[blen]) / xnorm;
		if (pg->Symmetric)
			Zm[0] = Zsym;
		else {
			znorm = qromb(SinNorm, gTheta[0], gTheta[blen]);
			Zm[0] = qromb(ZSinm, gTheta[0], gTheta[blen]) / znorm;
		}

		del = DMAX(fabs(Xm[0] - gX0), fabs(Zm[0] - gZ0));
		if (del < BNDDIFFTHR)
			break;
	}
	gX0 = Xm[0];
	gZ0 = Zm[0];

	FindTheta();

	/* Fit spline to x[theta] and z[theta] */
	dXdt = 0.5 * ((gX[blen] - gX[blen - 1]) / (gTheta[blen] - gTheta[blen - 1]) +
				  (gX[1] - gX[0]) / (gTheta[1] - gTheta[0]));
	spline(gTheta - 1, gX - 1, blen + 1, dXdt, dXdt, gXsplines - 1);

	dZdt = 0.5 * ((gZ[blen] - gZ[blen - 1]) / (gTheta[blen] - gTheta[blen - 1]) +
				  (gZ[1] - gZ[0]) / (gTheta[1] - gTheta[0]));
	spline(gTheta - 1, gZ - 1, blen + 1, dZdt, dZdt, gZsplines - 1);

	/* Evaluate other moments using X0 & Z0 */
	for (i = 1; i <= m; i++) {
		gm = i;
		xnorm = qromb(CosNorm, gTheta[0], gTheta[blen]);
		Xm[i] = qromb(XCosm, gTheta[0], gTheta[blen]) / xnorm;
		znorm = qromb(SinNorm, gTheta[0], gTheta[blen]);
		Zm[i] = qromb(ZSinm, gTheta[0], gTheta[blen]) / znorm;
	}

	free_dvector(gZsplines, 0, blen);
	free_dvector(gXsplines, 0, blen);
	free_dvector(gTheta, 0, blen);
	free_dvector(gZ, 0, blen);
	free_dvector(gX, 0, blen);
}
