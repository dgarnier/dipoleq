/*
**	nrRomberg.c
**
**	Romberg integration from Numerical Recipies.
**	May 4, 1993
**
**	Modifications:
**		May 12, 1993		Added absolute accuracy
**
**	Main routine:
**
**	integral = qromb(function, a, b);
**	integral = qsimp(function, a, b);
**
**	where
**			function is defined as
**				double function(double x);
**
**			a, b are the limits of integration
**
**	Modified for double precision math by M. E. Mauel
**	Dept. of Applied Physics
**
*/

#include <math.h>
#define NRANSI
#include "nrutil.h"
#include "nrRomberg.h"

#define FUNC(x) 	((*func)(x))

/*
**	EPS		==		Relative accuracy
**	ABS_EPS	==		Absolute accuracy
*/

#define EPS 		1.0e-5
#define ABS_EPS		1.0e-6
#define JMAX 		20
#define JMAXP 		(JMAX+1)
#define K 			5


#ifdef __cplusplus
extern "C" {
#endif

void          polint(double xa[], double ya[], int n, double x, double *y, double *dy);
double        trapzd(double (*func) (double), double a, double b, int n);


#ifdef __cplusplus
}
#endif

/*
**	p o l i n t
**
*/
void          polint(double *xa, double *ya, int n, double x, double *y, double *dy)
{
	int           i, m, ns = 1;
	double        den, dif, dift, ho, hp, w;
	double       *ca, *da;

	dif = fabs(x - xa[1]);
	ca = dvector(1, n);
	da = dvector(1, n);
	for (i = 1; i <= n; i++) {
		if ((dift = fabs(x - xa[i])) < dif) {
			ns = i;
			dif = dift;
		}
		ca[i] = ya[i];
		da[i] = ya[i];
	}
	*y = ya[ns--];
	for (m = 1; m < n; m++) {
		for (i = 1; i <= n - m; i++) {
			ho = xa[i] - x;
			hp = xa[i + m] - x;
			w = ca[i + 1] - da[i];
			den = xa[i] - xa[i + m];
			if (den == 0.0)	nrerror("Error in routine polint");
			den = w / den;
			da[i] = hp * den;
			ca[i] = ho * den;
		}
		if ( 2 * ns < (n - m) ) {
			*y += *dy = ca[ns + 1] ;
		} else {
			*y += *dy = da[ns--] ;
		}
	}

	free_dvector(da, 1, n);
	free_dvector(ca, 1, n);
}

/*
**	t r a p z d
**
*/
double        trapzd(double (*func) (double), double a, double b, int n)
{
	double        x, tnm, sum, del;
	static double s;
	int           it, j;

	if (n == 1) {
		return (s = 0.5 * (b - a) * (FUNC(a) + FUNC(b)));
	} else {
		for (it = 1, j = 1; j < n - 1; j++)
			it <<= 1;
		tnm = it;
		del = (b - a) / tnm;
		x = a + 0.5 * del;
		for (sum = 0.0, j = 1; j <= it; j++, x += del)
			sum += FUNC(x);
		s = 0.5 * (s + (b - a) * sum / tnm);
		return s;
	}
}

/*
** q s i m p
**
*/
double        qsimp(double (*func) (double), double a, double b)
{
	int           j;
	double        s, st, ost, os;

	ost = os = -1.0e30;
	for (j = 1; j <= JMAX; j++) {
		st = trapzd(func, a, b, j);
		s = (4.0 * st - ost) / 3.0;
		if ((fabs(s - os) < EPS * fabs(os)) || (fabs(s - os) < ABS_EPS))
			return s;
		os = s;
		ost = st;
	}
	nrinfo("Too many steps in routine qsimp");
	return s;
}

/*
** q r o m b
**
*/
double        qromb(double (*func) (double), double a, double b)
{
	double        ss, dss;
	double        s[JMAXP + 1], h[JMAXP + 1];
	int           j;

	h[1] = 1.0;
	for (j = 1; j <= JMAX; j++) {
		s[j] = trapzd(func, a, b, j);
		if (j >= K) {
			polint(&h[j - K], &s[j - K], K, 0.0, &ss, &dss);
			if ((fabs(dss) < EPS * fabs(ss)) || (fabs(dss) < ABS_EPS))
				return ss;
		}
		s[j + 1] = s[j];
		h[j + 1] = 0.25 * h[j];
	}
	nrinfo("Too many steps in routine qromb");
	return ss;
}

#undef EPS
#undef ABS_EPS
#undef JMAX
#undef JMAXP
#undef K
#undef NRANSI
#undef FUNC
