/*
**	nrSpline.c
**
**	Cubic Spline Interpolation from Numerical Recipies.
**
**	Modified for double precision math by M. E. Mauel
**	Dept. of Applied Physics
**	May 4, 1993
**
**	Modifications:
**
**		Sept. 24, 1993		Added 2D bi-cubic splines with derivatives.
**
*/

#define NRANSI
#include "nrutil.h"
#include "nrSpline.h"

/*
**	s p l i n e
**
*/
void  spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
{
	int           i, k;
	double        p, qn, sig, un, *u;

	u = dvector(1, n - 1);
	if (yp1 > 0.99e30)
		y2[1] = u[1] = 0.0;
	else {
		y2[1] = -0.5;
		u[1] = (3.0 / (x[2] - x[1])) * ((y[2] - y[1]) / (x[2] - x[1]) - yp1);
	}
	for (i = 2; i <= n - 1; i++) {
		sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
		p = sig * y2[i - 1] + 2.0;
		y2[i] = (sig - 1.0) / p;
		u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
		u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
	}
	if (ypn > 0.99e30)
		qn = un = 0.0;
	else {
		qn = 0.5;
		un = (3.0 / (x[n] - x[n - 1])) * (ypn - (y[n] - y[n - 1]) / (x[n] - x[n - 1]));
	}
	y2[n] = (un - qn * u[n - 1]) / (qn * y2[n - 1] + 1.0);
	for (k = n - 1; k >= 1; k--)
		y2[k] = y2[k] * y2[k + 1] + u[k];
	free_dvector(u, 1, n - 1);
}

/*
**	s p l i n t
**
*/
void          splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
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
		nrerror("Bad xa input to routine splint");
	a = (xa[khi] - x) / h;
	b = (x - xa[klo]) / h;
	*y = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
}

/*
**	s p l i n t _ d e r v s
**
**  Similar to splint excepting it returns the gradient of y at x as well
**	as y(x).
*/
void          splint_dervs(double xa[], double ya[], double y2a[], int n, double x, double *y, double *dy)
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
		nrerror("Bad xa input to routine splint");
	a = (xa[khi] - x) / h;
	b = (x - xa[klo]) / h;
	*y = a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
	*dy = (ya[khi] - ya[klo]) / h - (3.0 * a * a - 1.0) * h * y2a[klo] / 6.0 + (3.0 * b * b - 1.0) * h * y2a[khi] / 6.0;
}

/*
**	s p l i e 2
**
*/
void          splie2(double x1a[], double x2a[], double **ya, int m, int n, double **y2a)
{
#pragma unused( x1a )
	int           j;

	for (j = 1; j <= m; j++)
		spline(x2a, ya[j], n, 1.0e30, 1.0e30, y2a[j]);
}

/*
**	s p l i e 2 _ a l t
**
**	This version of splie2 computes the splines for the alternate direction
**	of the array **ya.  NOTE: The output array must have the OPPOSITE dimensions
**	as **ya!  i.e.  y2a[1..n][1..m].
**
*/
void          splie2_alt(double x1a[], double x2a[], double **ya, int m, int n, double **y2a)
{
#pragma unused( x2a )
	int           i, j;
	double       *yv;

	yv = dvector(1, m);

	for (i = 1; i <= n; i++) {
		for (j = 1; j <= m; j++)
			yv[j] = ya[j][i];
		spline(x1a, yv, m, 1.0e30, 1.0e30, y2a[i]);
	}

	free_dvector(yv, 1, m);
}

/*
**	s p l i n 2
**
*/
void          splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n,
					 double x1, double x2, double *y)
{
	int           j;
	double       *ytmp, *yytmp;

	ytmp = dvector(1, n);
	yytmp = dvector(1, n);

	for (j = 1; j <= m; j++)
		splint(x2a, ya[j], y2a[j], n, x2, &yytmp[j]);
	spline(x1a, yytmp, m, 1.0e30, 1.0e30, ytmp);
	splint(x1a, yytmp, ytmp, m, x1, y);

	free_dvector(yytmp, 1, n);
	free_dvector(ytmp, 1, n);
}

/*
**	s p l i n 2 _ d e r v s
**
**	Modified to return both the value of y(x1,x2) and the two derivatives.
**  It must do a spline in two dimensions in order to get both derivatives.
**
**	Notice that y2aa must have interchanged dimensions.
**
*/
void          splin2_dervs(double x1a[], double x2a[], double **ya, double **y2a, double **y2aa,
 int m, int n, double x1, double x2, double *y, double *dy1, double *dy2)
{
	int           i, j, mx;
	double       *ytmp, *yytmp;
	double        y1, y2;

	mx = IMAX(n, m);
	ytmp = dvector(1, mx);
	yytmp = dvector(1, mx);

	/* first, the x1 direction... */
	for (j = 1; j <= m; j++)
		splint(x2a, ya[j], y2a[j], n, x2, &yytmp[j]);
	spline(x1a, yytmp, m, 1.0e30, 1.0e30, ytmp);
	splint_dervs(x1a, yytmp, ytmp, m, x1, &y1, dy1);

	/* then, the x2 direction... */
	for (j = 1; j <= n; j++) {
		for (i = 1; i <= m; i++)
			ytmp[i] = ya[i][j];
		splint(x1a, ytmp, y2aa[j], m, x1, &yytmp[j]);
	}
	spline(x2a, yytmp, n, 1.0e30, 1.0e30, ytmp);
	splint_dervs(x1a, yytmp, ytmp, n, x2, &y2, dy2);

	*y = 0.5 * (y1 + y2);		/* average */

	free_dvector(yytmp, 1, mx);
	free_dvector(ytmp, 1, mx);
}

#undef NRANSI
