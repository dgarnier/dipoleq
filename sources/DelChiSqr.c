/*
**	DelChiSqr.c
**
**	This routine returns the value of Delta-ChiSqr
**	corresponding to a confidence level of "p" and
**	a degree of freedom, "nu".
**
**	From Numerical Recipies in C, 2nd Ed., p. 697.
**
**	M. E. Mauel
**	Columbia University
**
**	October 18, 1993
**
*/

#include <math.h>
#include "nrutil.h"
#include "DelChiSqr.h"

/*
**	Root Bisection
**
**
**	Find the root of (*func) between (x1, x2) to +-xacc.
*/

#ifdef __cplusplus
extern "C" {
#endif

double        rtbis(double (*func) (double), double x1, double x2, double xacc);
double        gammln(double xx);
void          gser(double *gamser, double a, double x, double *gln);
void          gcf(double *gammcf, double a, double x, double *gln);
double        FDEL(double del);

#ifdef __cplusplus
}
#endif

#define JMAX 40

double        rtbis(double (*func) (double), double x1, double x2, double xacc)
{
	int           j;
	double        dx, f, fmid, xmid, rtb;

	f = (*func) (x1);
	fmid = (*func) (x2);
	if (f * fmid >= 0.0)
		nrerror("Root must be bracketed for bisection in rtbis");
	rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
	for (j = 1; j <= JMAX; j++) {
		fmid = (*func) (xmid = rtb + (dx *= 0.5));
		if (fmid <= 0.0)
			rtb = xmid;
		if (fabs(dx) < xacc || fmid == 0.0)
			return rtb;
	}
	nrerror("Too many bisections in rtbis");
	return 0.0;
}

#undef JMAX

/*
**	Gamma Function
**
*/

#define 	ITMAX 		100
#define 	EPS 		3.0e-7
#define 	FPMIN 		1.0e-30

double        gammln(double xx)
{
	double        x, y, tmp, ser;
	static double cof[6] =
	{76.18009172947146, -86.50532032941677,
	 24.01409824083091, -1.231739572450155,
	 0.1208650973866179e-2, -0.5395239384953e-5};
	int           j;

	y = x = xx;
	tmp = x + 5.5;
	tmp -= (x + 0.5) * log(tmp);
	ser = 1.000000000190015;
	for (j = 0; j <= 5; j++)
		ser += cof[j] / ++y;
	return -tmp + log(2.5066282746310005 * ser / x);
}

void          gser(double *gamser, double a, double x, double *gln)
{
	int           n;
	double        sum, del, ap;

	*gln = gammln(a);
	if (x <= 0.0) {
		if (x < 0.0)
			nrerror("x less than 0 in routine gser");
		*gamser = 0.0;
		return;
	} else {
		ap = a;
		del = sum = 1.0 / a;
		for (n = 1; n <= ITMAX; n++) {
			++ap;
			del *= x / ap;
			sum += del;
			if (fabs(del) < fabs(sum) * EPS) {
				*gamser = sum * exp(-x + a * log(x) - (*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}

void          gcf(double *gammcf, double a, double x, double *gln)
{
	int           i;
	double        an, b, c, d, del, h;

	*gln = gammln(a);
	b = x + 1.0 - a;
	c = 1.0 / FPMIN;
	d = 1.0 / b;
	h = d;
	for (i = 1; i <= ITMAX; i++) {
		an = -i * (i - a);
		b += 2.0;
		d = an * d + b;
		if (fabs(d) < FPMIN)
			d = FPMIN;
		c = b + an / c;
		if (fabs(c) < FPMIN)
			c = FPMIN;
		d = 1.0 / d;
		del = d * c;
		h *= del;
		if (fabs(del - 1.0) < EPS)
			break;
	}
	if (i > ITMAX)
		nrerror("a too large, ITMAX too small in gcf");
	*gammcf = exp(-x + a * log(x) - (*gln)) * h;
}

#undef ITMAX
#undef EPS
#undef FPMIN

double        gammq(double a, double x)
{
	double        gamser, gammcf, gln;

	if (x < 0.0 || a <= 0.0)
		nrerror("Invalid arguments in routine gammq");
	if (x < (a + 1.0)) {
		gser(&gamser, a, x, &gln);
		return 1.0 - gamser;
	} else {
		gcf(&gammcf, a, x, &gln);
		return gammcf;
	}
}

/*
**	D E L C H I S Q R
**
**
**	The value of DelChiSqr for nu degrees of freedom corresponding
**	to p confidence level is the root of the equation
**
**		FDEL(Del) == gammq(nu/2,Del/2) + p - 1.0 = 0.0
**
**
*/

double        G_p;				/* global value of p */
double        G_nu;				/* nu/2.0 */

double        FDEL(double del)
{
	return gammq(G_nu, del / 2.0) + G_p - 1.0;
}

#define DACC	1.0e-4

double        DelChiSqr(double p, int nu)
{
	G_p = p;
	G_nu = nu / 2.0;

	return rtbis(FDEL, 1.0, 100.0, DACC);
}

#undef DACC
