/*
**	TokaMac v2.0
**
**	SVDFit.c
**
**	Routines from Numerical Recipes
**	Used to compute the LeastSquares best fit
**  by applying Singular Value Decomposition.
**
**
**	File:		SVDFit.c
**	Date:		March 21, 1993
**
**	Modifications:
**
**		Oct 3, 1993		Added svdvar to calculate the covariance matrix.
**
**
**	(c) L. Bai and M. Mauel -- Columbia University
*/

#include <math.h>
#include "nrutil.h"
#include "SVDFit.h"

#ifdef __cplusplus
extern "C" {
#endif

double pythag(double, double);
#ifdef __cplusplus
}
#endif

#define TOL 		1.0e-5

/*
**	pythag
**
**	computes sqrt of (a*a+b*b) without destructive overflow or underflow
**
*/
double        pythag(double a, double b)
{
	double        at, bt, ct;
	double        p;

	at = fabs(a);
	bt = fabs(b);
	if (at > bt) {
		ct = bt / at;
		p = at * sqrt(1.0 + ct * ct);
	} else {
		if (bt != 0.0) {
			ct = at / bt;
			p = bt * sqrt(1.0 + ct * ct);
		} else
			p = 0.0;
	}

	return p;
}

/*
**
**	svdcmp
** 	Singular Value Decomposition
**
** 	From NUMERICAL RECIPES IN C, Ch. 2.
**
** 	Given a matrix a[1...m][1...n], this routine computes its singular value decomposition,
**
** 	a ==> u.w.vt
**
** 	The matrix a is replaced on output by the matrix u.
** 	The diagonal matrix w is output as a vector w[1...n].
** 	The matrix v (not the transpose, vt) is output as v[1...n][1...n];
** 	if it is smaller,then a should be filled up to square with zero row.
**
**
*/
void          svdcmp(double **a, int m, int n, double *w, double **v)
{
	int           flag, i, its, j, jj, k, l, nm;
	double        c, f, h, s, x, y, z;
	double        anorm = 0.0, g = 0.0, scale = 0.0;
	double       *rv1;

	if (m < n)
		nrerror("SVDCMP: You must augment A with extra zero rows");
	rv1 = dvector(1, n);
	for (i = 1; i <= n; i++) {
		l = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m) {
			for (k = i; k <= m; k++)
				scale += fabs(a[k][i]);
			if (scale != 0.0) {
				for (k = i; k <= m; k++) {
					a[k][i] /= scale;
					s += a[k][i] * a[k][i];
				}
				f = a[i][i];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][i] = f - g;
				if (i != n) {
					for (j = l; j <= n; j++) {
						for (s = 0.0, k = i; k <= m; k++)
							s += a[k][i] * a[k][j];
						f = s / h;
						for (k = i; k <= m; k++)
							a[k][j] += f * a[k][i];
					}
				}
				for (k = i; k <= m; k++)
					a[k][i] *= scale;
			}
		}
		w[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m && i != n) {
			for (k = l; k <= n; k++)
				scale += fabs(a[i][k]);
			if (scale != 0.0) {
				for (k = l; k <= n; k++) {
					a[i][k] /= scale;
					s += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][l] = f - g;
				for (k = l; k <= n; k++)
					rv1[k] = a[i][k] / h;
				if (i != m) {
					for (j = l; j <= m; j++) {
						for (s = 0.0, k = l; k <= n; k++)
							s += a[j][k] * a[i][k];
						for (k = l; k <= n; k++)
							a[j][k] += s * rv1[k];
					}
				}
				for (k = l; k <= n; k++)
					a[i][k] *= scale;
			}
		}
		anorm = DMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
	}
	for (i = n; i >= 1; i--) {
		if (i < n) {
			if (g != 0.0) {
				for (j = l; j <= n; j++)
					v[j][i] = (a[i][j] / a[i][l]) / g;
				for (j = l; j <= n; j++) {
					for (s = 0.0, k = l; k <= n; k++)
						s += a[i][k] * v[k][j];
					for (k = l; k <= n; k++)
						v[k][j] += s * v[k][i];
				}
			}
			for (j = l; j <= n; j++)
				v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}
	for (i = n; i >= 1; i--) {
		l = i + 1;
		g = w[i];
		if (i < n)
			for (j = l; j <= n; j++)
				a[i][j] = 0.0;
		if (g != 0.0) {
			g = 1.0 / g;
			if (i != n) {
				for (j = l; j <= n; j++) {
					for (s = 0.0, k = l; k <= m; k++)
						s += a[k][i] * a[k][j];
					f = (s / a[i][i]) * g;
					for (k = i; k <= m; k++)
						a[k][j] += f * a[k][i];
				}
			}
			for (j = i; j <= m; j++)
				a[j][i] *= g;
		} else {
			for (j = i; j <= m; j++)
				a[j][i] = 0.0;
		}
		a[i][i] += 1.0;
	}
	for (k = n; k >= 1; k--) {
		for (its = 1; its <= 30; its++) {
			flag = 1;
			for (l = k; l >= 1; l--) {
				nm = l - 1;
				if (fabs(rv1[l]) + anorm == anorm) {
					flag = 0;
					break;
				}
				if (fabs(w[nm]) + anorm == anorm)
					break;
			}
			/* the following error should never happen M.E.M. 4/6/93 */
			if (flag && (nm < 1))
				nrerror("SVDCMP:	Diagonalization error");
			if (flag) {
				c = 0.0;
				s = 1.0;
				for (i = l; i <= k; i++) {
					f = s * rv1[i];
					if (fabs(f) + anorm != anorm) {
						g = w[i];
						h = pythag(f, g);
						w[i] = h;
						h = 1.0 / h;
						c = g * h;
						s = (-f * h);
						for (j = 1; j <= m; j++) {
							y = a[j][nm];
							z = a[j][i];
							a[j][nm] = y * c + z * s;
							a[j][i] = z * c - y * s;
						}
					}
				}
			}
			z = w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j = 1; j <= n; j++)
						v[j][k] = (-v[j][k]);
				}
				break;
			}
			if (its == 30)
				nrerror("SVDCMP:	No convergence in 30 iterations");
			x = w[l];
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
			g = pythag(f, 1.0);
			f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
			c = s = 1.0;
			for (j = l; j <= nm; j++) {
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y = y * c;
				for (jj = 1; jj <= n; jj++) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x * c + z * s;
					v[jj][i] = z * c - x * s;
				}
				z = pythag(f, h);
				w[j] = z;
				if (z != 0.0) {
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = (c * g) + (s * y);
				x = (c * y) - (s * g);
				for (jj = 1; jj <= m; jj++) {
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y * c + z * s;
					a[jj][i] = z * c - y * s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
	free_dvector(rv1, 1, n);
}

/*
**
**	svbksb
** 	From Ch. 2 of NUMERICAL RECIPES.
**
** 	Solves the equation
**
** 	a.x = b
**
** 	for a vector, x.
** 	Where a is specified by the arrays u,w,v as returned by SVDCMP.
** 	B is the input right-hand-side. It has dimensions of NumMeas = m.
** 	X is the output solution vector. It has dimensions of NumUnkn = n.
** 	No input quantities are destroyed.
**
*/
void          svbksb(double **u, double *w, double **v, int m, int n,
					 double *b, double *x)
{
	int           jj, j, i;
	double        s, *tmp;

	tmp = dvector(1, n);
	for (j = 1; j <= n; j++) {
		s = 0.0;
		if (w[j] != 0.0) {
			for (i = 1; i <= m; i++)
				s += u[i][j] * b[i];
			s /= w[j];
		}
		tmp[j] = s;
	}
	for (j = 1; j <= n; j++) {
		s = 0.0;
		for (jj = 1; jj <= n; jj++)
			s += v[j][jj] * tmp[jj];
		x[j] = s;
	}
	free_dvector(tmp, 1, n);
}

/*
**
**	svdfit
** 	Modified from Ch. 14 of NUMERICAL RECIPES.
**
** 	This finds the least square fit using SVDCMP.
** 	This solves the minimization
**
** 	chisq = |L.y - b|**2
**
** 	where L is input,
** 	b contains the "measurements" (and is of dimensions NumMeas).
** 	This procedure first completes a SVDCMP.  This replaces u,v,w with the
** 	decomposition matricies.
**
** 	The final output is contained in the vector y...which contains the best fit.
** 	On output u,v,w can be used to calculate the covariance, etc.
*/
void          svdfit(double **L, double **u, double *w, double **v, double *y,
					 double *b, int NumMeas, int NumUnkn, double *chisq)
{
	int           j, i;
	double        wmax;
	double        thresh;
	double        sum;

	for (i = 1; i <= NumMeas; i++)
		for (j = 1; j <= NumUnkn; j++)
			u[i][j] = L[i][j];

	svdcmp(u, NumMeas, NumUnkn, w, v);
	wmax = 0.0;
	for (j = 1; j <= NumUnkn; j++)
		wmax = DMAX(w[j], wmax);
	thresh = TOL * wmax;
	for (j = 1; j <= NumUnkn; j++)
		if (w[j] < thresh)
			w[j] = 0.0;
	svbksb(u, w, v, NumMeas, NumUnkn, b, y);
	*chisq = 0.0;
	for (i = 1; i <= NumMeas; i++) {
		sum = 0.0;
		for (j = 1; j <= NumUnkn; j++)
			sum += L[i][j] * y[j];
		*chisq += DSQR(b[i] - sum);
	}
}

void          svdvar(double **v, int ma, double *w, double **cvm)
{
	int           k, j, i;
	double        sum, *wti;

	wti = dvector(1, ma);
	for (i = 1; i <= ma; i++) {
		wti[i] = 0.0;
		if (w[i])
			wti[i] = 1.0 / (w[i] * w[i]);
	}
	for (i = 1; i <= ma; i++) {
		for (j = 1; j <= i; j++) {
			for (sum = 0.0, k = 1; k <= ma; k++)
				sum += v[i][k] * v[j][k] * wti[k];
			cvm[j][i] = cvm[i][j] = sum;
		}
	}
	free_dvector(wti, 1, ma);
}

#undef TOL
