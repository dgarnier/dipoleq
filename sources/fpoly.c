/*
** TokaMac v2.0
**
** fpoly.c
**
**
**
** File:		fpoly.c
** Date:		March 22, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdlib.h>
#include <math.h>
#include "fpoly.h"

/*
**	fpoly
**
** 	we expand  the flux function F() in terms of this parameter X
**
**    X = (Psi-PsiAxis)/(PsiLim-PsiAxis)
**
**	F(X) = F0 + ( F1 + X F2 + x^2 F3...) - X^n (F1 + F2 +...+ Fn)
**
** 	F0 is equal to the value of this flux function at the edge where X=1
**
*/
double        fpoly(double *f, double x, int nTerms)
{
	int           i;
	double        t = 0.0;
	double        s;

	s = f[0];
	for (i = 1; i <= nTerms; i++) {
//		t += f[i];
		s += f[i] * pow(x, i - 1);
	}

	return (s - t * pow(x, nTerms));
}

/*
**	fpoly_int
**
** 	The integral of fpoly.
**
**	We need to know DelPsi and fpoly_int(X=1).
**
**	fpoly_int = - Int(F' dx, x, 1.0)
**			  = F(1) - DelPsi*[F'(edge) (1 - X) + Sum of terms]
**
*/
double        fpoly_int(double *f, double x, int nTerms,
						double DelPsi, double f1)
{
	int           i;
	double        t = 0.0;
	double        s;

	s = f[0] * (1.0 - x);
	for (i = 1; i <= nTerms; i++) {
//		t += f[i];
		s += f[i] * (1.0 - pow(x, i)) / i;
	}

	s = s - t * (1.0 - pow(x, nTerms + 1)) / (nTerms + 1);
	return (f1 - DelPsi * s);
}
