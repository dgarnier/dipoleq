/*
** TokaMac v2.0
**
** green.c
**
** Utility functions for the poloidal field Green's Function
** in cylindrical coordinates.
**
**
** File:		green.c
** Date:		April 20, 1992
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <math.h>
#include "green.h"

/*Constants for EllipK...*/
const double  ka0 = 1.38629436112;
const double  ka1 = 0.09666344259;
const double  ka2 = 0.03590092383;
const double  ka3 = 0.03742563713;
const double  ka4 = 0.01451196212;
const double  kb0 = 0.5;
const double  kb1 = 0.12498593597;
const double  kb2 = 0.06880248576;
const double  kb3 = 0.03328355346;
const double  kb4 = 0.00441787012;

/* Constants for EllipE... */
const double  ea1 = 0.44325141463;
const double  ea2 = 0.06260601220;
const double  ea3 = 0.04757383546;
const double  ea4 = 0.01736506451;
const double  eb1 = 0.24998368310;
const double  eb2 = 0.09200180037;
const double  eb3 = 0.04069697526;
const double  eb4 = 0.00526449639;

/*
** The limit ROFF is due to a limit of how close a grid point can
** be to a current element.  k2(max) < 1 - 0.25*(d/x)^2 wher d is the
** closest distance of separation.
**
** This limit seems to be important.  How big or small should it be?
**
** For the following,
**
**		(dZ/X) > 2 SQRT( ROFF ) ~ 1.0e-3
**
** For X = 1 m, this gives a 1 mm minimum separation.
**
*/
const double  ROFF = 0.0100 / (64 * 64) / 3.0;

/*
** This function is the cylindrical Green's function
** ** It does not have a Two Pi **
** This uses the definition of elliptic functions in Abromowitz and Segun.
*/

/*
** DTG 1/12/98 -- added fix for zero x
*/

double        Green(double x, double z, double xc, double zc)
{
	double        fGreen;
	double        k2, m1, EllipK, EllipE, dl, denom;

	denom =  ((x + xc) * (x + xc) + (z - zc) * (z - zc));
	k2 = 4.0 * x * xc / denom;

	if (k2 > (1.0 - ROFF))
		k2 = 1.0 - ROFF;
	m1 = 1.0 - k2;
	dl = log(1.0 / m1);
	EllipK = ka0 + m1 * (ka1 + m1 * (ka2 + m1 * (ka3 + m1 * ka4)))
		+ (kb0 + m1 * (kb1 + m1 * (kb2 + m1 * (kb3 + m1 * kb4)))) * dl;
	EllipE = 1.0 + m1 * (ea1 + m1 * (ea2 + m1 * (ea3 + m1 * ea4)))
		+ (m1 * (eb1 + m1 * (eb2 + m1 * (eb3 + m1 * eb4)))) * dl;

	fGreen = -sqrt(denom/4.0) * ((2.0 - k2) * EllipK - 2.0 * EllipE);

	return fGreen;
}

/*
** This function "GetdGreen" returns the gradient of the greens functions.
*/
void          GetdGreen(double *GG, double *dGx, double *dGz, double x, double z, double xc, double zc)
{
	double        k2, m1, EllipK, EllipE, dl, dKE, denom;

	denom =  ((x + xc) * (x + xc) + (z - zc) * (z - zc));
	k2 = 4.0 * x * xc / denom;
	if (k2 > (1.0 - ROFF))
		k2 = 1.0 - ROFF;
	m1 = 1.0 - k2;
	dl = log(1.0 / m1);
	EllipK = ka0 + m1 * (ka1 + m1 * (ka2 + m1 * (ka3 + m1 * ka4)))
		+ (kb0 + m1 * (kb1 + m1 * (kb2 + m1 * (kb3 + m1 * kb4)))) * dl;
	EllipE = 1.0 + m1 * (ea1 + m1 * (ea2 + m1 * (ea3 + m1 * ea4)))
		+ (m1 * (eb1 + m1 * (eb2 + m1 * (eb3 + m1 * eb4)))) * dl;
	*GG = -sqrt(denom/4.0) * ((2.0 - k2) * EllipK - 2.0 * EllipE);
	dKE = EllipE / m1 - EllipK;
	*dGx = 0.25 * (*GG * 4.0 * (x + xc) / denom  - sqrt(4.0/denom) * (2.0 * xc - k2 * (x + xc)) * dKE);
	*dGz = ((z - zc) / denom) * (*GG + sqrt(k2 * x * xc) * dKE);
}
