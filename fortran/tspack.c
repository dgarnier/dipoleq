/* tspack.f -- translated by f2c (version 19971204).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/*
#include "f2c.h"
*/

#include "tspack.h"

/* Common Block Declarations */

Extern struct {
    doublereal y;
} stcom_;

#define stcom_1 stcom_

/* Table of constant values */

static doublereal c_b82 = 1.;
static doublereal c_b192 = 0.;


/*                        TSPT1 */
/*                       10/05/98 */

/*   This program provides a quick easily-verified test of */
/* all TSPACK modules except the high-level interface */
/* routines TSPxx and TSVALx.  The primary test uses data */
/* values taken from the quadratic function f(x) = x**2 for */
/* which knot-derivative estimation and interpolation with no */
/* tension is exact.  Since this test function is positive, */
/* strictly increasing, and convex, zero tension is suffici- */
/* ent to preserve these properties in the data.  The */
/* abscissae are taken to be a set of N points uniformly */
/* distributed in [0,1] */

/*   The maximum absolute error in the derivative estimates */
/* is printed for each of the following. */

/*      YPC1 (local approximation of knot-derivatives) */
/*      YPC2 with end conditions from YPC1 estimates */
/*      YPC2 with specified endpoint first derivatives */
/*      YPC2 with specified endpoint second derivatives */
/*      YPC2 with end conditions computed by ENDSLP */

/*   The maximum error in a tension factor is printed for */
/* each of the following.  These should all be zero since no */
/* tension is required in any of the cases. */

/*      SIGBI (minimum tension factor required to satisfy */
/*             constraints defined by array B) */
/*      SIGS (minimum tension factor required to preserve */
/*            monotonicity and convexity on each interval) */
/*      SIG0 with a lower bound of Y(I) on (X(I),X(I+1)) */
/*      SIG1 with a lower bound of YP(I) on (X(I),X(I+1)) */
/*      SIG2 (minimum tension factor required to preserve */
/*            convexity on each interval) */

/*   The maximum absolute error on a set of 3*(N-1)+1 points */
/* uniformly distributed in [0,1] is printed for each of the */
/* following evaluation routines. */

/*      HVAL (function values) */
/*      HPVAL (first derivative values) */
/*      HPPVAL (second derivative values) */
/*      TSINTL (integral from 0 to T for each evaluation */
/*              point T) */

/*   The following four routines are used to construct para- */
/* metric tension spline fits to N points uniformly */
/* distributed on the unit circle.  In the case of SMCRV, */
/* periodic end conditions are used to fit cos(A) and natural */
/* end conditions are used for sin(A), where A is the parame- */
/* ter value (angle).  Thus, both B2TRI and B2TRIP are exer- */
/* cized.  The weights W are based on a standard deviation of */
/* EPS, for machine precision EPS, so that the smoothing */
/* curves are close to the interpolatory curves.  The maximum */
/* distance from the circle to a set of 3*(N-1)+1 evaluation */
/* points (angles uniformly distributed in [0,2*PI]) is */
/* printed for each routine. */

/*      YPC1P (local approximations to knot-derivatives with */
/*             periodic end conditions) */
/*      YPC2P (global approximations to knot-derivatives with */
/*             periodic end conditions) */
/*      SIGBP (minimum tension factor required to satisfy */
/*             bounds on distance between line segments and */
/*             planar curve segments) */
/*      SMCRV (global approximations to both function values */
/*             and first derivatives at the knots) */

/*   The final test consists of computing a set of parameter */
/* values associated with the N points on the unit circle. */
/* These are computed by both ARCL2D and ARCL3D with constant */
/* Z values, and the maximum difference is printed. */

/*      ARCL2D (cumulative arc lengths for a planar curve) */
/*      ARCL3D (cumulative arc lengths for a space curve) */

/* Subroutine */ int arcl2d_(integer *n, doublereal *x, doublereal *y, 
	doublereal *t, integer *ier)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal ds;
    static integer nn;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/10/92 */

/*   Given an ordered sequence of points (X,Y) defining a */
/* polygonal curve in the plane, this subroutine computes the */
/* sequence T of cumulative arc lengths along the curve: */
/* T(1) = 0 and, for 2 .LE. K .LE. N, T(K) is the sum of */
/* Euclidean distances between (X(I-1),Y(I-1)) and (X(I),Y(I)) */
/* for I = 2,...,K.  A closed curve corresponds to X(1) = */
/* X(N) and Y(1) = Y(N), and more generally, duplicate points */
/* are permitted but must not be adjacent.  Thus, T contains */
/* a strictly increasing sequence of values which may be used */
/* as parameters for fitting a smooth curve to the sequence */
/* of points. */

/* On input: */

/*       N = Number of points defining the curve.  N .GE. 2. */

/*       X,Y = Arrays of length N containing the coordinates */
/*             of the points. */

/* The above parameters are not altered by this routine. */

/*       T = Array of length at least N. */

/* On output: */

/*       T = Array containing cumulative arc lengths defined */
/*           above unless IER > 0. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if N < 2. */
/*             IER = I if X(I-1) = X(I) and Y(I-1) = Y(I) for */
/*                     some I in the range 2,...,N. */

/* Modules required by ARCL2D:  None */

/* Intrinsic function called by ARCL2D:  SQRT */

/* *********************************************************** */


    /* Parameter adjustments */
    --t;
    --y;
    --x;

    /* Function Body */
    nn = *n;
    if (nn < 2) {
	goto L2;
    }
    t[1] = 0.;
    i__1 = nn;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = x[i__] - x[i__ - 1];
/* Computing 2nd power */
	d__2 = y[i__] - y[i__ - 1];
	ds = d__1 * d__1 + d__2 * d__2;
	if (ds == 0.) {
	    goto L3;
	}
	t[i__] = t[i__ - 1] + sqrt(ds);
/* L1: */
    }
    *ier = 0;
    return 0;

/* N is outside its valid range. */

L2:
    *ier = 1;
    return 0;

/* Points I-1 and I coincide. */

L3:
    *ier = i__;
    return 0;
} /* arcl2d_ */

/* Subroutine */ int arcl3d_(integer *n, doublereal *x, doublereal *y, 
	doublereal *z__, doublereal *t, integer *ier)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal ds;
    static integer nn;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/10/92 */

/*   Given an ordered sequence of points (X,Y,Z) defining a */
/* polygonal curve in 3-space, this subroutine computes the */
/* sequence T of cumulative arc lengths along the curve: */
/* T(1) = 0 and, for 2 .LE. K .LE. N, T(K) is the sum of */
/* Euclidean distances between (X(I-1),Y(I-1),Z(I-1)) and */
/* (X(I),Y(I),Z(I)) for I = 2,...,K.  A closed curve corre- */
/* sponds to X(1) = X(N), Y(1) = Y(N), and Z(1) = Z(N).  More */
/* generally, duplicate points are permitted but must not be */
/* adjacent.  Thus, T contains a strictly increasing sequence */
/* of values which may be used as parameters for fitting a */
/* smooth curve to the sequence of points. */

/* On input: */

/*       N = Number of points defining the curve.  N .GE. 2. */

/*       X,Y,Z = Arrays of length N containing the coordi- */
/*               nates of the points. */

/* The above parameters are not altered by this routine. */

/*       T = Array of length at least N. */

/* On output: */

/*       T = Array containing cumulative arc lengths defined */
/*           above unless IER > 0. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if N < 2. */
/*             IER = I if X(I-1) = X(I), Y(I-1) = Y(I), and */
/*                     Z(I-1) = Z(I) for some I in the range */
/*                     2,...,N. */

/* Modules required by ARCL3D:  None */

/* Intrinsic function called by ARCL3D:  SQRT */

/* *********************************************************** */


    /* Parameter adjustments */
    --t;
    --z__;
    --y;
    --x;

    /* Function Body */
    nn = *n;
    if (nn < 2) {
	goto L2;
    }
    t[1] = 0.;
    i__1 = nn;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = x[i__] - x[i__ - 1];
/* Computing 2nd power */
	d__2 = y[i__] - y[i__ - 1];
/* Computing 2nd power */
	d__3 = z__[i__] - z__[i__ - 1];
	ds = d__1 * d__1 + d__2 * d__2 + d__3 * d__3;
	if (ds == 0.) {
	    goto L3;
	}
	t[i__] = t[i__ - 1] + sqrt(ds);
/* L1: */
    }
    *ier = 0;
    return 0;

/* N is outside its valid range. */

L2:
    *ier = 1;
    return 0;

/* Points I-1 and I coincide. */

L3:
    *ier = i__;
    return 0;
} /* arcl3d_ */

/* Subroutine */ int b2tri_(integer *n, doublereal *x, doublereal *y, 
	doublereal *w, doublereal *p, doublereal *d__, doublereal *sd, 
	doublereal *t11, doublereal *t12, doublereal *t21, doublereal *t22, 
	doublereal *ys, doublereal *yp, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal s11im1, s12im1, s22im1;
    static integer i__;
    static doublereal r1, r2, di;
    static integer nn;
    static doublereal dx, pp;
    static integer im1, nm1;
    static doublereal d11i, d12i, d22i, den, s11i, s12i, s22i, dim1;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/10/92 */

/*   This subroutine solves the order 2N symmetric positive- */
/* definite block tridiagonal linear system associated with */
/* minimizing the quadratic functional Q(YS,YP) described in */
/* Subroutine SMCRV. */

/* On input: */

/*       N = Number of data points.  N .GE. 2. */

/*       X,Y,W = Arrays of length N containing abscissae, */
/*               data values, and positive weights, respect- */
/*               ively.  The abscissae must be strictly in- */
/*               creasing. */

/*       P = Positive smoothing parameter defining Q. */

/*       D,SD = Arrays of length N-1 containing positive ma- */
/*              trix entries.  Letting DX and SIG denote the */
/*              width and tension factor associated with the */
/*              interval (X(I),X(I+1)), D(I) = SIG*(SIG* */
/*              COSHM(SIG) - SINHM(SIG))/(DX*E) and SD(I) = */
/*              SIG*SINHM(SIG)/(DX*E) where E = SIG*SINH(SIG) */
/*              - 2*COSHM(SIG). */

/* The above parameters are not altered by this routine. */

/*       T11,T12,T21,T22 = Arrays of length N-1 used as */
/*                         temporary work space. */

/* On output: */

/*       YS,YP = Arrays of length N containing solution com- */
/*               ponents:  function and derivative values, */
/*               respectively, at the abscissae. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if N or P is outside its valid range */
/*                     on input. */
/*             Note that no test is made for a nonpositive */
/*             value of X(I+1)-X(I), W(I), D(I), or SD(I). */

/* Modules required by B2TRI:  None */

/* *********************************************************** */


    /* Parameter adjustments */
    --yp;
    --ys;
    --t22;
    --t21;
    --t12;
    --t11;
    --sd;
    --d__;
    --w;
    --y;
    --x;

    /* Function Body */
    nn = *n;
    nm1 = nn - 1;
    pp = *p;
    *ier = 1;
    if (nn < 2 || pp <= 0.) {
	return 0;
    }

/* The forward elimination step consists of scaling a row by */
/*   the inverse of its diagonal block and eliminating the */
/*   subdiagonal block.  The superdiagonal is stored in T and */
/*   the right hand side in YS,YP.  For J = 11, 12, and 22, */
/*   SJI and SJIM1 denote the elements in position J of the */
/*   superdiagonal block in rows I and I-1, respectively. */
/*   Similarly, DJI denotes an element in the diagonal block */
/*   of row I. */

/* Initialize for I = 2. */

    dx = x[2] - x[1];
    dim1 = d__[1];
    s22im1 = sd[1];
    s12im1 = (dim1 + s22im1) / dx;
    s11im1 = s12im1 * -2. / dx;
    r1 = pp * w[1];
    d11i = r1 - s11im1;
    d12i = s12im1;
    d22i = dim1;
    den = d11i * d22i - d12i * d12i;
    t11[1] = (d22i * s11im1 + d12i * s12im1) / den;
    t12[1] = (d22i * s12im1 - d12i * s22im1) / den;
    t21[1] = -(d12i * s11im1 + d11i * s12im1) / den;
    t22[1] = (d11i * s22im1 - d12i * s12im1) / den;
    r1 = r1 * y[1] / den;
    ys[1] = d22i * r1;
    yp[1] = -d12i * r1;

/* I = 2,...,N-1: */

    i__1 = nm1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	im1 = i__ - 1;
	dx = x[i__ + 1] - x[i__];
	di = d__[i__];
	s22i = sd[i__];
	s12i = (di + s22i) / dx;
	s11i = s12i * -2. / dx;
	r1 = pp * w[i__];
	d11i = r1 - s11im1 - s11i - (s11im1 * t11[im1] - s12im1 * t21[im1]);
	d12i = s12i - s12im1 - (s11im1 * t12[im1] - s12im1 * t22[im1]);
	d22i = dim1 + di - (s12im1 * t12[im1] + s22im1 * t22[im1]);
	den = d11i * d22i - d12i * d12i;
	t11[i__] = (d22i * s11i + d12i * s12i) / den;
	t12[i__] = (d22i * s12i - d12i * s22i) / den;
	t21[i__] = -(d12i * s11i + d11i * s12i) / den;
	t22[i__] = (d11i * s22i - d12i * s12i) / den;
	r1 = r1 * y[i__] - s11im1 * ys[im1] + s12im1 * yp[im1];
	r2 = -s12im1 * ys[im1] - s22im1 * yp[im1];
	ys[i__] = (d22i * r1 - d12i * r2) / den;
	yp[i__] = (d11i * r2 - d12i * r1) / den;
	dim1 = di;
	s22im1 = s22i;
	s12im1 = s12i;
	s11im1 = s11i;
/* L1: */
    }

/* I = N: */

    r1 = pp * w[nn];
    d11i = r1 - s11im1 - (s11im1 * t11[nm1] - s12im1 * t21[nm1]);
    d12i = -s12im1 - (s11im1 * t12[nm1] - s12im1 * t22[nm1]);
    d22i = dim1 - (s12im1 * t12[nm1] + s22im1 * t22[nm1]);
    den = d11i * d22i - d12i * d12i;
    r1 = r1 * y[nn] - s11im1 * ys[nm1] + s12im1 * yp[nm1];
    r2 = -s12im1 * ys[nm1] - s22im1 * yp[nm1];
    ys[nn] = (d22i * r1 - d12i * r2) / den;
    yp[nn] = (d11i * r2 - d12i * r1) / den;

/* Back solve the system. */

    for (i__ = nm1; i__ >= 1; --i__) {
	ys[i__] -= t11[i__] * ys[i__ + 1] + t12[i__] * yp[i__ + 1];
	yp[i__] -= t21[i__] * ys[i__ + 1] + t22[i__] * yp[i__ + 1];
/* L2: */
    }
    *ier = 0;
    return 0;
} /* b2tri_ */

/* Subroutine */ int b2trip_(integer *n, doublereal *x, doublereal *y, 
	doublereal *w, doublereal *p, doublereal *d__, doublereal *sd, 
	doublereal *t11, doublereal *t12, doublereal *t21, doublereal *t22, 
	doublereal *u11, doublereal *u12, doublereal *u21, doublereal *u22, 
	doublereal *ys, doublereal *yp, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal s11im1, s12im1, s22im1, s11nm1, s12nm1, s22nm1, ypnm1, 
	    ysnm1;
    static integer i__;
    static doublereal r1, r2, di;
    static integer nn;
    static doublereal dx, pp;
    static integer im1, ip1, nm1, nm2, nm3;
    static doublereal d11i, d12i, d22i, den, s11i, s12i, s22i, su11, su12, 
	    su21, su22, dim1, dnm1;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/10/92 */

/*   This subroutine solves the order 2(N-1) symmetric posi- */
/* tive-definite linear system associated with minimizing the */
/* quadratic functional Q(YS,YP) (described in Subroutine */
/* SMCRV) with periodic end conditions.  The matrix is block */
/* tridiagonal except for nonzero blocks in the upper right */
/* and lower left corners. */

/* On input: */

/*       N = Number of data points.  N .GE. 3. */

/*       X = Array of length N containing a strictly in- */
/*           creasing sequence of abscissae. */

/*       Y,W = Arrays of length N-1 containing data values */
/*             and positive weights, respectively, associated */
/*             with the first N-1 abscissae. */

/*       P = Positive smoothing parameter defining Q. */

/*       D,SD = Arrays of length N-1 containing positive ma- */
/*              trix elements.  Letting DX and SIG denote the */
/*              width and tension factor associated with the */
/*              interval (X(I),X(I+1)), D(I) = SIG*(SIG* */
/*              COSHM(SIG) - SINHM(SIG))/(DX*E) and SD(I) = */
/*              SIG*SINHM(SIG)/(DX*E) where E = SIG*SINH(SIG) */
/*              - 2*COSHM(SIG). */

/* The above parameters are not altered by this routine. */

/*       T11,T12,T21,T22,U11,U12,U21,U22 = Arrays of length */
/*                                         N-2 used as temp- */
/*                                         orary work space. */

/* On output: */

/*       YS,YP = Arrays of length N containing solution com- */
/*               ponents:  function and derivative values, */
/*               respectively, at the abscissae.  YS(N) = */
/*               YS(1) and YP(N) = YP(1). */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if N or P is outside its valid range */
/*                     on input. */
/*             Note that no test is made for a nonpositive */
/*             value of X(I+1)-X(I), W(I), D(I), or SD(I). */

/* Modules required by B2TRIP:  None */

/* *********************************************************** */


    /* Parameter adjustments */
    --yp;
    --ys;
    --u22;
    --u21;
    --u12;
    --u11;
    --t22;
    --t21;
    --t12;
    --t11;
    --sd;
    --d__;
    --w;
    --y;
    --x;

    /* Function Body */
    nn = *n;
    nm1 = nn - 1;
    nm2 = nn - 2;
    nm3 = nn - 3;
    pp = *p;
    *ier = 1;
    if (nn < 3 || pp <= 0.) {
	return 0;
    }

/* The forward elimination step consists of scaling a row by */
/*   the inverse of its diagonal block and eliminating the */
/*   subdiagonal block for the first N-2 rows.  The super- */
/*   diagonal is stored in T, the negative of the last column */
/*   in U, and the right hand side in YS,YP.  For J = 11, 12, */
/*   and 22, SJI and SJIM1 denote the elements in position J */
/*   of the superdiagonal block in rows I and I-1, respect- */
/*   ively.  Similarly, DJI denotes an element in the diago- */
/*   nal block of row I. */

/* I = 1: */

    dx = x[nn] - x[nm1];
    dnm1 = d__[nm1];
    s22nm1 = sd[nm1];
    s12nm1 = -(dnm1 + s22nm1) / dx;
    s11nm1 = s12nm1 * 2. / dx;
    dx = x[2] - x[1];
    di = d__[1];
    s22i = sd[1];
    s12i = (di + s22i) / dx;
    s11i = s12i * -2. / dx;
    r1 = pp * w[1];
    d11i = r1 - s11nm1 - s11i;
    d12i = s12i + s12nm1;
    d22i = dnm1 + di;
    den = d11i * d22i - d12i * d12i;
    t11[1] = (d22i * s11i + d12i * s12i) / den;
    t12[1] = (d22i * s12i - d12i * s22i) / den;
    t21[1] = -(d12i * s11i + d11i * s12i) / den;
    t22[1] = (d11i * s22i - d12i * s12i) / den;
    u11[1] = -(d22i * s11nm1 + d12i * s12nm1) / den;
    u12[1] = (d12i * s22nm1 - d22i * s12nm1) / den;
    u21[1] = (d12i * s11nm1 + d11i * s12nm1) / den;
    u22[1] = (d12i * s12nm1 - d11i * s22nm1) / den;
    r1 = r1 * y[1] / den;
    ys[1] = d22i * r1;
    yp[1] = -d12i * r1;
    if (nn == 3) {
	goto L2;
    }

/* I = 2,...,N-2: */

    i__1 = nm2;
    for (i__ = 2; i__ <= i__1; ++i__) {
	im1 = i__ - 1;
	dim1 = di;
	s22im1 = s22i;
	s12im1 = s12i;
	s11im1 = s11i;
	dx = x[i__ + 1] - x[i__];
	di = d__[i__];
	s22i = sd[i__];
	s12i = (di + s22i) / dx;
	s11i = s12i * -2. / dx;
	r1 = pp * w[i__];
	d11i = r1 - s11im1 - s11i - (s11im1 * t11[im1] - s12im1 * t21[im1]);
	d12i = s12i - s12im1 - (s11im1 * t12[im1] - s12im1 * t22[im1]);
	d22i = dim1 + di - (s12im1 * t12[im1] + s22im1 * t22[im1]);
	den = d11i * d22i - d12i * d12i;
	t11[i__] = (d22i * s11i + d12i * s12i) / den;
	t12[i__] = (d22i * s12i - d12i * s22i) / den;
	t21[i__] = -(d12i * s11i + d11i * s12i) / den;
	t22[i__] = (d11i * s22i - d12i * s12i) / den;
	su11 = s11im1 * u11[im1] - s12im1 * u21[im1];
	su12 = s11im1 * u12[im1] - s12im1 * u22[im1];
	su21 = s12im1 * u11[im1] + s22im1 * u21[im1];
	su22 = s12im1 * u12[im1] + s22im1 * u22[im1];
	u11[i__] = (d12i * su21 - d22i * su11) / den;
	u12[i__] = (d12i * su22 - d22i * su12) / den;
	u21[i__] = (d12i * su11 - d11i * su21) / den;
	u22[i__] = (d12i * su12 - d11i * su22) / den;
	r1 = r1 * y[i__] - s11im1 * ys[im1] + s12im1 * yp[im1];
	r2 = -s12im1 * ys[im1] - s22im1 * yp[im1];
	ys[i__] = (d22i * r1 - d12i * r2) / den;
	yp[i__] = (d11i * r2 - d12i * r1) / den;
/* L1: */
    }

/* The backward elimination step zeros the first N-3 blocks */
/*   of the superdiagonal.  For I = N-2,N-3,...,1, T(I) and */
/*   (YS(I),YP(I)) are overwritten with the negative of the */
/*   last column and the new right hand side, respectively. */

L2:
    t11[nm2] = u11[nm2] - t11[nm2];
    t12[nm2] = u12[nm2] - t12[nm2];
    t21[nm2] = u21[nm2] - t21[nm2];
    t22[nm2] = u22[nm2] - t22[nm2];
    for (i__ = nm3; i__ >= 1; --i__) {
	ip1 = i__ + 1;
	ys[i__] = ys[i__] - t11[i__] * ys[ip1] - t12[i__] * yp[ip1];
	yp[i__] = yp[i__] - t21[i__] * ys[ip1] - t22[i__] * yp[ip1];
	t11[i__] = u11[i__] - t11[i__] * t11[ip1] - t12[i__] * t21[ip1];
	t12[i__] = u12[i__] - t11[i__] * t12[ip1] - t12[i__] * t22[ip1];
	t21[i__] = u21[i__] - t21[i__] * t11[ip1] - t22[i__] * t21[ip1];
	t22[i__] = u22[i__] - t21[i__] * t12[ip1] - t22[i__] * t22[ip1];
/* L3: */
    }

/* Solve the last equation for YS(N-1),YP(N-1).  SJI = SJNM2 */
/*   and DJI = DJNM1. */

    r1 = pp * w[nm1];
    d11i = r1 - s11i - s11nm1 + s11nm1 * t11[1] - s12nm1 * t21[1] + s11i * 
	    t11[nm2] - s12i * t21[nm2];
    d12i = -s12nm1 - s12i + s11nm1 * t12[1] - s12nm1 * t22[1] + s11i * t12[
	    nm2] - s12i * t22[nm2];
    d22i = di + dnm1 + s12nm1 * t12[1] + s22nm1 * t22[1] + s12i * t12[nm2] + 
	    s22i * t22[nm2];
    den = d11i * d22i - d12i * d12i;
    r1 = r1 * y[nm1] - s11nm1 * ys[1] + s12nm1 * yp[1] - s11i * ys[nm2] + 
	    s12i * yp[nm2];
    r2 = -s12nm1 * ys[1] - s22nm1 * yp[1] - s12i * ys[nm2] - s22i * yp[nm2];
    ysnm1 = (d22i * r1 - d12i * r2) / den;
    ypnm1 = (d11i * r2 - d12i * r1) / den;
    ys[nm1] = ysnm1;
    yp[nm1] = ypnm1;

/* Back substitute for the remainder of the solution */
/*   components. */

    i__1 = nm2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ys[i__] = ys[i__] + t11[i__] * ysnm1 + t12[i__] * ypnm1;
	yp[i__] = yp[i__] + t21[i__] * ysnm1 + t22[i__] * ypnm1;
/* L4: */
    }

/* YS(N) = YS(1) and YP(N) = YP(1). */

    ys[nn] = ys[1];
    yp[nn] = yp[1];
    *ier = 0;
    return 0;
} /* b2trip_ */

doublereal endslp_(doublereal *x1, doublereal *x2, doublereal *x3, doublereal 
	*y1, doublereal *y2, doublereal *y3, doublereal *sigma)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal sigs, t, dummy, s1, coshm1;
    extern /* Subroutine */ int snhcsh_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal coshms, dx1, dxs, sig1;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   11/17/96 */

/*   Given data values associated with a strictly increasing */
/* or decreasing sequence of three abscissae X1, X2, and X3, */
/* this function returns a derivative estimate at X1 based on */
/* the tension spline H(x) which interpolates the data points */
/* and has third derivative equal to zero at X1.  Letting S1 */
/* denote the slope defined by the first two points, the est- */
/* mate is obtained by constraining the derivative of H at X1 */
/* so that it has the same sign as S1 and its magnitude is */
/* at most 3*abs(S1).  If SIGMA = 0, H(x) is quadratic and */
/* the derivative estimate is identical to the value computed */
/* by Subroutine YPC1 at the first point (or the last point */
/* if the abscissae are decreasing). */

/* On input: */

/*       X1,X2,X3 = Abscissae satisfying either X1 < X2 < X3 */
/*                  or X1 > X2 > X3. */

/*       Y1,Y2,Y3 = Data values associated with the abscis- */
/*                  sae.  H(X1) = Y1, H(X2) = Y2, and H(X3) */
/*                  = Y3. */

/*       SIGMA = Tension factor associated with H in inter- */
/*               val (X1,X2) or (X2,X1). */

/* Input parameters are not altered by this function. */

/* On output: */

/*       ENDSLP = (Constrained) derivative of H at X1, or */
/*                zero if the abscissae are not strictly */
/*                monotonic. */

/* Module required by ENDSLP:  SNHCSH */

/* Intrinsic functions called by ENDSLP:  ABS, EXP, MAX, MIN */

/* *********************************************************** */


    dx1 = *x2 - *x1;
    dxs = *x3 - *x1;
    if (dx1 * (dxs - dx1) <= 0.) {
	goto L2;
    }
    sig1 = abs(*sigma);
    if (sig1 < 1e-9) {

/* SIGMA = 0:  H is the quadratic interpolant. */

/* Computing 2nd power */
	d__1 = dx1 / dxs;
	t = d__1 * d__1;
	goto L1;
    }
    sigs = sig1 * dxs / dx1;
    if (sigs <= .5) {

/* 0 < SIG1 < SIGS .LE. .5:  compute approximations to */
/*   COSHM1 = COSH(SIG1)-1 and COSHMS = COSH(SIGS)-1. */

	snhcsh_(&sig1, &dummy, &coshm1, &dummy);
	snhcsh_(&sigs, &dummy, &coshms, &dummy);
	t = coshm1 / coshms;
    } else {

/* SIGS > .5:  compute T = COSHM1/COSHMS. */

/* Computing 2nd power */
	d__1 = (1. - exp(-sig1)) / (1. - exp(-sigs));
	t = exp(sig1 - sigs) * (d__1 * d__1);
    }

/* The derivative of H at X1 is */
/*   T = ((Y3-Y1)*COSHM1-(Y2-Y1)*COSHMS)/ */
/*       (DXS*COSHM1-DX1*COSHMS). */

/* ENDSLP = T unless T*S1 < 0 or abs(T) > 3*abs(S1). */

L1:
    t = ((*y3 - *y1) * t - *y2 + *y1) / (dxs * t - dx1);
    s1 = (*y2 - *y1) / dx1;
    if (s1 >= 0.) {
/* Computing MIN */
	d__1 = max(0.,t), d__2 = s1 * 3.;
	ret_val = min(d__1,d__2);
    } else {
/* Computing MAX */
	d__1 = min(0.,t), d__2 = s1 * 3.;
	ret_val = max(d__1,d__2);
    }
    return ret_val;

/* Error in the abscissae. */

L2:
    ret_val = 0.;
    return ret_val;
} /* endslp_ */

doublereal hppval_(doublereal *t, integer *n, doublereal *x, doublereal *y, 
	doublereal *yp, doublereal *sigma, integer *ier)
{
    /* Initialized data */

    static doublereal sbig = 85.;

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal cosh2, sinh2, e;
    static integer i__;
    static doublereal s, b1, b2, d1, d2, e1, e2, dummy, cm, dx, sm, tm;
    extern /* Subroutine */ int snhcsh_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal cm2, sb1, sb2;
    static integer ip1;
    extern integer intrvl_(doublereal *, integer *, doublereal *);
    static doublereal sm2, cmm, sig, ems;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   11/17/96 */

/*   This function evaluates the second derivative HPP of a */
/* Hermite interpolatory tension spline H at a point T. */

/* On input: */

/*       T = Point at which HPP is to be evaluated.  Extrap- */
/*           olation is performed if T < X(1) or T > X(N). */

/*       N = Number of data points.  N .GE. 2. */

/*       X = Array of length N containing the abscissae. */
/*           These must be in strictly increasing order: */
/*           X(I) < X(I+1) for I = 1,...,N-1. */

/*       Y = Array of length N containing data values. */
/*           H(X(I)) = Y(I) for I = 1,...,N. */

/*       YP = Array of length N containing first deriva- */
/*            tives.  HP(X(I)) = YP(I) for I = 1,...,N, where */
/*            HP denotes the derivative of H. */

/*       SIGMA = Array of length N-1 containing tension fac- */
/*               tors whose absolute values determine the */
/*               balance between cubic and linear in each */
/*               interval.  SIGMA(I) is associated with int- */
/*               erval (I,I+1) for I = 1,...,N-1. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered and */
/*                     X(1) .LE. T .LE. X(N). */
/*             IER = 1 if no errors were encountered and */
/*                     extrapolation was necessary. */
/*             IER = -1 if N < 2. */
/*             IER = -2 if the abscissae are not in strictly */
/*                      increasing order.  (This error will */
/*                      not necessarily be detected.) */

/*       HPPVAL = Second derivative value HPP(T), or zero if */
/*                IER < 0. */

/* Modules required by HPPVAL:  INTRVL, SNHCSH */

/* Intrinsic functions called by HPPVAL:  ABS, EXP */

/* *********************************************************** */


    /* Parameter adjustments */
    --sigma;
    --yp;
    --y;
    --x;

    /* Function Body */
    if (*n < 2) {
	goto L1;
    }

/* Find the index of the left end of an interval containing */
/*   T.  If T < X(1) or T > X(N), extrapolation is performed */
/*   using the leftmost or rightmost interval. */

    if (*t < x[1]) {
	i__ = 1;
	*ier = 1;
    } else if (*t > x[*n]) {
	i__ = *n - 1;
	*ier = 1;
    } else {
	i__ = intrvl_(t, n, &x[1]);
	*ier = 0;
    }
    ip1 = i__ + 1;

/* Compute interval width DX, local coordinates B1 and B2, */
/*   and second differences D1 and D2. */

    dx = x[ip1] - x[i__];
    if (dx <= 0.) {
	goto L2;
    }
    b1 = (x[ip1] - *t) / dx;
    b2 = 1. - b1;
    s = (y[ip1] - y[i__]) / dx;
    d1 = s - yp[i__];
    d2 = yp[ip1] - s;
    sig = (d__1 = sigma[i__], abs(d__1));
    if (sig < 1e-9) {

/* SIG = 0:  H is the Hermite cubic interpolant. */

	ret_val = (d1 + d2 + (b2 - b1) * 3. * (d2 - d1)) / dx;
    } else if (sig <= .5) {

/* 0 .LT. SIG .LE. .5:  use approximations designed to avoid */
/*   cancellation error in the hyperbolic functions. */

	sb2 = sig * b2;
	snhcsh_(&sig, &sm, &cm, &cmm);
	snhcsh_(&sb2, &sm2, &cm2, &dummy);
	sinh2 = sm2 + sb2;
	cosh2 = cm2 + 1.;
	e = sig * sm - cmm - cmm;
	ret_val = sig * ((cm * sinh2 - sm * cosh2) * (d1 + d2) + sig * (cm * 
		cosh2 - (sm + sig) * sinh2) * d1) / (dx * e);
    } else {

/* SIG > .5:  use negative exponentials in order to avoid */
/*   overflow.  Note that EMS = EXP(-SIG).  In the case of */
/*   extrapolation (negative B1 or B2), H is approximated by */
/*   a linear function if -SIG*B1 or -SIG*B2 is large. */

	sb1 = sig * b1;
	sb2 = sig - sb1;
	if (-sb1 > sbig || -sb2 > sbig) {
	    ret_val = 0.;
	} else {
	    e1 = exp(-sb1);
	    e2 = exp(-sb2);
	    ems = e1 * e2;
	    tm = 1. - ems;
	    e = tm * (sig * (ems + 1.) - tm - tm);
	    ret_val = sig * (sig * ((e1 * ems + e2) * d1 + (e1 + e2 * ems) * 
		    d2) - tm * (e1 + e2) * (d1 + d2)) / (dx * e);
	}
    }
    return ret_val;

/* N is outside its valid range. */

L1:
    ret_val = 0.;
    *ier = -1;
    return ret_val;

/* X(I) .GE. X(I+1). */

L2:
    ret_val = 0.;
    *ier = -2;
    return ret_val;
} /* hppval_ */

doublereal hpval_(doublereal *t, integer *n, doublereal *x, doublereal *y, 
	doublereal *yp, doublereal *sigma, integer *ier)
{
    /* Initialized data */

    static doublereal sbig = 85.;

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal sinh2, e;
    static integer i__;
    static doublereal s, b1, b2, d1, d2, e1, e2, dummy, s1, cm, dx, sm, tm;
    extern /* Subroutine */ int snhcsh_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal cm2, sb1, sb2;
    static integer ip1;
    extern integer intrvl_(doublereal *, integer *, doublereal *);
    static doublereal sm2, cmm, sig, ems;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   11/17/96 */

/*   This function evaluates the first derivative HP of a */
/* Hermite interpolatory tension spline H at a point T. */

/* On input: */

/*       T = Point at which HP is to be evaluated.  Extrapo- */
/*           lation is performed if T < X(1) or T > X(N). */

/*       N = Number of data points.  N .GE. 2. */

/*       X = Array of length N containing the abscissae. */
/*           These must be in strictly increasing order: */
/*           X(I) < X(I+1) for I = 1,...,N-1. */

/*       Y = Array of length N containing data values. */
/*           H(X(I)) = Y(I) for I = 1,...,N. */

/*       YP = Array of length N containing first deriva- */
/*            tives.  HP(X(I)) = YP(I) for I = 1,...,N. */

/*       SIGMA = Array of length N-1 containing tension fac- */
/*               tors whose absolute values determine the */
/*               balance between cubic and linear in each */
/*               interval.  SIGMA(I) is associated with int- */
/*               erval (I,I+1) for I = 1,...,N-1. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       IER = Error indicator: */
/*             IER = 0  if no errors were encountered and */
/*                      X(1) .LE. T .LE. X(N). */
/*             IER = 1  if no errors were encountered and */
/*                      extrapolation was necessary. */
/*             IER = -1 if N < 2. */
/*             IER = -2 if the abscissae are not in strictly */
/*                      increasing order.  (This error will */
/*                      not necessarily be detected.) */

/*       HPVAL = Derivative value HP(T), or zero if IER < 0. */

/* Modules required by HPVAL:  INTRVL, SNHCSH */

/* Intrinsic functions called by HPVAL:  ABS, EXP */

/* *********************************************************** */


    /* Parameter adjustments */
    --sigma;
    --yp;
    --y;
    --x;

    /* Function Body */
    if (*n < 2) {
	goto L1;
    }

/* Find the index of the left end of an interval containing */
/*   T.  If T < X(1) or T > X(N), extrapolation is performed */
/*   using the leftmost or rightmost interval. */

    if (*t < x[1]) {
	i__ = 1;
	*ier = 1;
    } else if (*t > x[*n]) {
	i__ = *n - 1;
	*ier = 1;
    } else {
	i__ = intrvl_(t, n, &x[1]);
	*ier = 0;
    }
    ip1 = i__ + 1;

/* Compute interval width DX, local coordinates B1 and B2, */
/*   and second differences D1 and D2. */

    dx = x[ip1] - x[i__];
    if (dx <= 0.) {
	goto L2;
    }
    b1 = (x[ip1] - *t) / dx;
    b2 = 1. - b1;
    s1 = yp[i__];
    s = (y[ip1] - y[i__]) / dx;
    d1 = s - s1;
    d2 = yp[ip1] - s;
    sig = (d__1 = sigma[i__], abs(d__1));
    if (sig < 1e-9) {

/* SIG = 0:  H is the Hermite cubic interpolant. */

	ret_val = s1 + b2 * (d1 + d2 - b1 * 3. * (d2 - d1));
    } else if (sig <= .5) {

/* 0 .LT. SIG .LE. .5:  use approximations designed to avoid */
/*   cancellation error in the hyperbolic functions. */

	sb2 = sig * b2;
	snhcsh_(&sig, &sm, &cm, &cmm);
	snhcsh_(&sb2, &sm2, &cm2, &dummy);
	sinh2 = sm2 + sb2;
	e = sig * sm - cmm - cmm;
	ret_val = s1 + ((cm * cm2 - sm * sinh2) * (d1 + d2) + sig * (cm * 
		sinh2 - (sm + sig) * cm2) * d1) / e;
    } else {

/* SIG > .5:  use negative exponentials in order to avoid */
/*   overflow.  Note that EMS = EXP(-SIG).  In the case of */
/*   extrapolation (negative B1 or B2), H is approximated by */
/*   a linear function if -SIG*B1 or -SIG*B2 is large. */

	sb1 = sig * b1;
	sb2 = sig - sb1;
	if (-sb1 > sbig || -sb2 > sbig) {
	    ret_val = s;
	} else {
	    e1 = exp(-sb1);
	    e2 = exp(-sb2);
	    ems = e1 * e2;
	    tm = 1. - ems;
	    e = tm * (sig * (ems + 1.) - tm - tm);
	    ret_val = s + (tm * ((e2 - e1) * (d1 + d2) + tm * (d1 - d2)) + 
		    sig * ((e1 * ems - e2) * d1 + (e1 - e2 * ems) * d2)) / e;
	}
    }
    return ret_val;

/* N is outside its valid range. */

L1:
    ret_val = 0.;
    *ier = -1;
    return ret_val;

/* X(I) .GE. X(I+1). */

L2:
    ret_val = 0.;
    *ier = -2;
    return ret_val;
} /* hpval_ */

doublereal hval_(doublereal *t, integer *n, doublereal *x, doublereal *y, 
	doublereal *yp, doublereal *sigma, integer *ier)
{
    /* Initialized data */

    static doublereal sbig = 85.;

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal e;
    static integer i__;
    static doublereal s, u, b1, b2, d1, d2, e1, e2, dummy, s1, y1, cm, dx, sm,
	     tm, tp, ts;
    extern /* Subroutine */ int snhcsh_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal cm2, sb1, sb2;
    static integer ip1;
    extern integer intrvl_(doublereal *, integer *, doublereal *);
    static doublereal sm2, cmm, sig, ems;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   11/17/96 */

/*   This function evaluates a Hermite interpolatory tension */
/* spline H at a point T.  Note that a large value of SIGMA */
/* may cause underflow.  The result is assumed to be zero. */

/*   Given arrays X, Y, YP, and SIGMA of length NN, if T is */
/* known to lie in the interval (X(I),X(J)) for some I < J, */
/* a gain in efficiency can be achieved by calling this */
/* function with N = J+1-I (rather than NN) and the I-th */
/* components of the arrays (rather than the first) as par- */
/* ameters. */

/* On input: */

/*       T = Point at which H is to be evaluated.  Extrapo- */
/*           lation is performed if T < X(1) or T > X(N). */

/*       N = Number of data points.  N .GE. 2. */

/*       X = Array of length N containing the abscissae. */
/*           These must be in strictly increasing order: */
/*           X(I) < X(I+1) for I = 1,...,N-1. */

/*       Y = Array of length N containing data values. */
/*           H(X(I)) = Y(I) for I = 1,...,N. */

/*       YP = Array of length N containing first deriva- */
/*            tives.  HP(X(I)) = YP(I) for I = 1,...,N, where */
/*            HP denotes the derivative of H. */

/*       SIGMA = Array of length N-1 containing tension fac- */
/*               tors whose absolute values determine the */
/*               balance between cubic and linear in each */
/*               interval.  SIGMA(I) is associated with int- */
/*               erval (I,I+1) for I = 1,...,N-1. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       IER = Error indicator: */
/*             IER = 0  if no errors were encountered and */
/*                      X(1) .LE. T .LE. X(N). */
/*             IER = 1  if no errors were encountered and */
/*                      extrapolation was necessary. */
/*             IER = -1 if N < 2. */
/*             IER = -2 if the abscissae are not in strictly */
/*                      increasing order.  (This error will */
/*                      not necessarily be detected.) */

/*       HVAL = Function value H(T), or zero if IER < 0. */

/* Modules required by HVAL:  INTRVL, SNHCSH */

/* Intrinsic functions called by HVAL:  ABS, EXP */

/* *********************************************************** */


    /* Parameter adjustments */
    --sigma;
    --yp;
    --y;
    --x;

    /* Function Body */
    if (*n < 2) {
	goto L1;
    }

/* Find the index of the left end of an interval containing */
/*   T.  If T < X(1) or T > X(N), extrapolation is performed */
/*   using the leftmost or rightmost interval. */

    if (*t < x[1]) {
	i__ = 1;
	*ier = 1;
    } else if (*t > x[*n]) {
	i__ = *n - 1;
	*ier = 1;
    } else {
	i__ = intrvl_(t, n, &x[1]);
	*ier = 0;
    }
    ip1 = i__ + 1;

/* Compute interval width DX, local coordinates B1 and B2, */
/*   and second differences D1 and D2. */

    dx = x[ip1] - x[i__];
    if (dx <= 0.) {
	goto L2;
    }
    u = *t - x[i__];
    b2 = u / dx;
    b1 = 1. - b2;
    y1 = y[i__];
    s1 = yp[i__];
    s = (y[ip1] - y1) / dx;
    d1 = s - s1;
    d2 = yp[ip1] - s;
    sig = (d__1 = sigma[i__], abs(d__1));
    if (sig < 1e-9) {

/* SIG = 0:  H is the Hermite cubic interpolant. */

	ret_val = y1 + u * (s1 + b2 * (d1 + b1 * (d1 - d2)));
    } else if (sig <= .5) {

/* 0 .LT. SIG .LE. .5:  use approximations designed to avoid */
/*   cancellation error in the hyperbolic functions. */

	sb2 = sig * b2;
	snhcsh_(&sig, &sm, &cm, &cmm);
	snhcsh_(&sb2, &sm2, &cm2, &dummy);
	e = sig * sm - cmm - cmm;
	ret_val = y1 + s1 * u + dx * ((cm * sm2 - sm * cm2) * (d1 + d2) + sig 
		* (cm * cm2 - (sm + sig) * sm2) * d1) / (sig * e);
    } else {

/* SIG > .5:  use negative exponentials in order to avoid */
/*   overflow.  Note that EMS = EXP(-SIG).  In the case of */
/*   extrapolation (negative B1 or B2), H is approximated by */
/*   a linear function if -SIG*B1 or -SIG*B2 is large. */

	sb1 = sig * b1;
	sb2 = sig - sb1;
	if (-sb1 > sbig || -sb2 > sbig) {
	    ret_val = y1 + s * u;
	} else {
	    e1 = exp(-sb1);
	    e2 = exp(-sb2);
	    ems = e1 * e2;
	    tm = 1. - ems;
	    ts = tm * tm;
	    tp = ems + 1.;
	    e = tm * (sig * tp - tm - tm);
	    ret_val = y1 + s * u + dx * (tm * (tp - e1 - e2) * (d1 + d2) + 
		    sig * ((e2 + ems * (e1 - 2.) - b1 * ts) * d1 + (e1 + ems *
		     (e2 - 2.) - b2 * ts) * d2)) / (sig * e);
	}
    }
    return ret_val;

/* N is outside its valid range. */

L1:
    ret_val = 0.;
    *ier = -1;
    return ret_val;

/* X(I) .GE. X(I+1). */

L2:
    ret_val = 0.;
    *ier = -2;
    return ret_val;
} /* hval_ */

integer intrvl_(doublereal *t, integer *n, doublereal *x)
{
    /* Initialized data */

    static integer il = 1;

    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer k, ih;
    static doublereal tt;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/10/92 */

/*   This function returns the index of the left end of an */
/* interval (defined by an increasing sequence X) which */
/* contains the value T.  The method consists of first test- */
/* ing the interval returned by a previous call, if any, and */
/* then using a binary search if necessary. */

/* On input: */

/*       T = Point to be located. */

/*       N = Length of X.  N .GE. 2. */

/*       X = Array of length N assumed (without a test) to */
/*           contain a strictly increasing sequence of */
/*           values. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       INTRVL = Index I defined as follows: */

/*                  I = 1    if  T .LT. X(2) or N .LE. 2, */
/*                  I = N-1  if  T .GE. X(N-1), and */
/*                  X(I) .LE. T .LT. X(I+1) otherwise. */

/* Modules required by INTRVL:  None */

/* *********************************************************** */


    /* Parameter adjustments */
    --x;

    /* Function Body */
    tt = *t;
    if (il >= 1 && il < *n) {
	if (x[il] <= tt && tt < x[il + 1]) {
	    goto L2;
	}
    }

/* Initialize low and high indexes. */

    il = 1;
    ih = *n;

/* Binary search: */

L1:
    if (ih <= il + 1) {
	goto L2;
    }
    k = (il + ih) / 2;
    if (tt < x[k]) {
	ih = k;
    } else {
	il = k;
    }
    goto L1;

/* X(IL) .LE. T .LT. X(IL+1)  or  (T .LT. X(1) and IL=1) */
/*                            or  (T .GE. X(N) and IL=N-1) */

L2:
    ret_val = il;
    return ret_val;
} /* intrvl_ */

doublereal sig0_(doublereal *x1, doublereal *x2, doublereal *y1, doublereal *
	y2, doublereal *y1p, doublereal *y2p, integer *ifl, doublereal *hbnd, 
	doublereal *tol, integer *ier)
{
    /* Initialized data */

    static doublereal sbig = 85.;

    /* System generated locals */
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *), exp(
	    doublereal), log(doublereal);

    /* Local variables */
    static doublereal d1pd2, fneg, dsig, dmax__, fmax, sneg, ftol, rsig, rtol,
	     stol, a, b, c__, d__, e, f, r__, s, t, coshm, sinhm, a0, b0, c1, 
	    c2, d0, d2, f0, ssinh;
    extern doublereal store_(doublereal *);
    static doublereal s1, s2, t0, t1, t2, aa, rf, dx, tm, coshmm;
    extern /* Subroutine */ int snhcsh_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal y1l, y2l, bnd, scm, sig, ems;
    static integer nit;
    static doublereal ssm;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   11/18/96 */

/*   Given a pair of abscissae with associated ordinates and */
/* slopes, this function determines the smallest (nonnega- */
/* tive) tension factor SIGMA such that the Hermite interpo- */
/* latory tension spline H(x), defined by SIGMA and the data, */
/* is bounded (either above or below) by HBND for all x in */
/* (X1,X2). */

/* On input: */

/*       X1,X2 = Abscissae.  X1 < X2. */

/*       Y1,Y2 = Values of H at X1 and X2. */

/*       Y1P,Y2P = Derivative values of H at X1 and X2. */

/*       IFL = Option indicator: */
/*             IFL = -1 if HBND is a lower bound on H. */
/*             IFL = 1 if HBND is an upper bound on H. */

/*       HBND = Bound on H.  If IFL = -1, HBND .LE. min(Y1, */
/*              Y2).  If IFL = 1, HBND .GE. max(Y1,Y2). */

/*       TOL = Tolerance whose magnitude determines how close */
/*             SIGMA is to its optimal value when nonzero */
/*             finite tension is necessary and sufficient to */
/*             satisfy the constraint.  For a lower bound, */
/*             SIGMA is chosen so that HBND .LE. HMIN .LE. */
/*             HBND + abs(TOL), where HMIN is the minimum */
/*             value of H in the interval, and for an upper */
/*             bound, the maximum of H satisfies HBND - */
/*             abs(TOL) .LE. HMAX .LE. HBND.  Thus, the con- */
/*             straint is satisfied but possibly with more */
/*             tension than necessary. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered and the */
/*                     constraint can be satisfied with fin- */
/*                     ite tension. */
/*             IER = 1 if no errors were encountered but in- */
/*                     finite tension is required to satisfy */
/*                     the constraint (e.g., IFL = -1, HBND */
/*                     = Y1, and Y1P < 0.). */
/*             IER = -1 if X2 .LE. X1 or abs(IFL) .NE. 1. */
/*             IER = -2 if HBND is outside its valid range */
/*                      on input. */

/*       SIG0 = Minimum tension factor defined above unless */
/*              IER < 0, in which case SIG0 = -1.  If IER = */
/*              1, SIG0 = 85, resulting in an approximation */
/*              to the linear interpolant of the endpoint */
/*              values.  Note, however, that SIG0 may be */
/*              larger than 85 if IER = 0. */

/* Modules required by SIG0:  SNHCSH, STORE */

/* Intrinsic functions called by SIG0:  ABS, DBLE, EXP, LOG, */
/*                                        MAX, MIN, SIGN, */
/*                                        SQRT */

/* *********************************************************** */



/* Store local parameters and test for errors. */

    rf = (doublereal) (*ifl);
    dx = *x2 - *x1;
    if (abs(rf) != 1. || dx <= 0.) {
	goto L8;
    }
    y1l = *y1;
    y2l = *y2;
    bnd = *hbnd;

/* Test for a valid constraint. */

    if (rf < 0. && min(y1l,y2l) < bnd || rf > 0. && bnd < max(y1l,y2l)) {
	goto L9;
    }

/* Test for infinite tension required. */

    s1 = *y1p;
    s2 = *y2p;
    if (y1l == bnd && rf * s1 > 0. || y2l == bnd && rf * s2 < 0.) {
	goto L7;
    }

/* Test for SIG = 0 sufficient. */

    sig = 0.;
    if (rf * s1 <= 0. && rf * s2 >= 0.) {
	goto L6;
    }

/*   Compute coefficients A0 and B0 of the Hermite cubic in- */
/*     terpolant H0(x) = Y2 - DX*(S2*R + B0*R**2 + A0*R**3/3) */
/*     where R = (X2-x)/DX. */

    s = (y2l - y1l) / dx;
    t0 = s * 3. - s1 - s2;
    a0 = (s - t0) * 3.;
    b0 = t0 - s2;
    d0 = t0 * t0 - s1 * s2;

/*   H0 has local extrema in (X1,X2) iff S1*S2 < 0 or */
/*     (T0*(S1+S2) < 0 and D0 .GE. 0). */

    if (s1 * s2 >= 0. && (t0 * (s1 + s2) >= 0. || d0 < 0.)) {
	goto L6;
    }
    if (a0 == 0.) {

/*   H0 is quadratic and has an extremum at R = -S2/(2*B0). */
/*     H0(R) = Y2 + DX*S2**2/(4*B0).  Note that A0 = 0 im- */
/*     plies 2*B0 = S1-S2, and S1*S2 < 0 implies B0 .NE. 0. */
/*     Also, the extremum is a min iff HBND is a lower bound. */

	f0 = (bnd - y2l - dx * s2 * s2 / (b0 * 4.)) * rf;
    } else {

/*   A0 .NE. 0 and H0 has extrema at R = (-B0 +/- SQRT(D0))/ */
/*     A0 = S2/(-B0 -/+ SQRT(D0)), where the negative root */
/*     corresponds to a min.  The expression for R is chosen */
/*     to avoid cancellation error.  H0(R) = Y2 + DX*(S2*B0 + */
/*     2*D0*R)/(3*A0). */

	d__1 = sqrt(d0);
	t = -b0 - d_sign(&d__1, &b0);
	r__ = t / a0;
	if (rf * b0 > 0.) {
	    r__ = s2 / t;
	}
	f0 = (bnd - y2l - dx * (s2 * b0 + d0 * 2. * r__) / (a0 * 3.)) * rf;
    }

/*   F0 .GE. 0 iff SIG = 0 is sufficient to satisfy the */
/*     constraint. */

    if (f0 >= 0.) {
	goto L6;
    }

/* Find a zero of F(SIG) = (BND-H(R))*RF where the derivative */
/*   of H, HP, vanishes at R.  F is a nondecreasing function, */
/*   F(0) < 0, and F = FMAX for SIG sufficiently large. */

/* Initialize parameters for the secant method.  The method */
/*   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG), */
/*   where SG0 and SNEG are defined implicitly by DSIG = SIG */
/*   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and */
/*   SNEG is initialized to a sufficiently large value that */
/*   FNEG > 0.  This value is used only if the initial value */
/*   of F is negative. */

/* Computing MAX */
/* Computing MIN */
    d__5 = (d__1 = y1l - bnd, abs(d__1)), d__6 = (d__2 = y2l - bnd, abs(d__2))
	    ;
    d__3 = .001, d__4 = min(d__5,d__6);
    fmax = max(d__3,d__4);
/* Computing MAX */
    d__3 = (d__1 = y1l - bnd, abs(d__1)), d__4 = (d__2 = y2l - bnd, abs(d__2))
	    ;
    t = max(d__3,d__4);
/* Computing MAX */
    d__1 = abs(s1), d__2 = abs(s2);
    sig = dx * max(d__1,d__2) / t;
    dmax__ = sig * (1. - t / fmax);
    sneg = sig - dmax__;
    dsig = sig;
    fneg = fmax;
    d2 = s2 - s;
    d1pd2 = s2 - s1;
    nit = 0;

/* Compute an absolute tolerance FTOL = abs(TOL), and a */
/*   relative tolerance RTOL = 100*MACHEPS. */

    ftol = abs(*tol);
    rtol = 1.;
L1:
    rtol /= 2.;
    d__1 = rtol + 1.;
    if (store_(&d__1) > 1.) {
	goto L1;
    }
    rtol *= 200.;

/* Top of loop:  compute F. */

L2:
    ems = exp(-sig);
    if (sig <= .5) {

/*   SIG .LE. .5:  use approximations designed to avoid can- */
/*                 cellation error (associated with small */
/*                 SIG) in the modified hyperbolic functions. */

	snhcsh_(&sig, &sinhm, &coshm, &coshmm);
	c1 = sig * coshm * d2 - sinhm * d1pd2;
	c2 = sig * (sinhm + sig) * d2 - coshm * d1pd2;
	a = c2 - c1;
	aa = a / ems;
	e = sig * sinhm - coshmm - coshmm;
    } else {

/*   SIG > .5:  scale SINHM and COSHM by 2*EXP(-SIG) in order */
/*              to avoid overflow. */

	tm = 1. - ems;
	ssinh = tm * (ems + 1.);
	ssm = ssinh - sig * 2. * ems;
	scm = tm * tm;
	c1 = sig * scm * d2 - ssm * d1pd2;
	c2 = sig * ssinh * d2 - scm * d1pd2;
	aa = (sig * tm * d2 + (tm - sig) * d1pd2) * 2.;
	a = ems * aa;
	e = sig * ssinh - scm - scm;
    }

/*   HP(R) = S2 - (C1*SINH(SIG*R) - C2*COSHM(SIG*R))/E = 0 */
/*     for ESR = (-B +/- SQRT(D))/A = C/(-B -/+ SQRT(D)), */
/*     where ESR = EXP(SIG*R), A = C2-C1, D = B**2 - A*C, and */
/*     B and C are defined below. */

    b = e * s2 - c2;
    c__ = c2 + c1;
    d__ = b * b - a * c__;
    f = 0.;
    if (aa * c__ == 0. && b == 0.) {
	goto L3;
    }
    f = fmax;
    if (d__ < 0.) {
	goto L3;
    }
    t1 = sqrt(d__);
    t = -b - d_sign(&t1, &b);
    rsig = 0.;
    if (rf * b < 0. && aa != 0.f) {
	if (t / aa > 0.) {
	    rsig = sig + log(t / aa);
	}
    }
    if ((rf * b > 0. || aa == 0.) && c__ / t > 0.) {
	rsig = log(c__ / t);
    }
    if ((rsig <= 0. || rsig >= sig) && b != 0.) {
	goto L3;
    }

/*   H(R) = Y2 - DX*(B*SIG*R + C1 + RF*SQRT(D))/(SIG*E). */

    f = (bnd - y2l + dx * (b * rsig + c1 + rf * t1) / (sig * e)) * rf;

/*   Update the number of iterations NIT. */

L3:
    ++nit;
    if (f0 * f < 0.) {

/*   F0*F < 0.  Update (SNEG,FNEG) to (SG0,F0) so that F */
/*     and FNEG always have opposite signs.  If SIG is */
/*     closer to SNEG than SG0, then swap (SNEG,FNEG) with */
/*     (SG0,F0). */

	t1 = dmax__;
	t2 = fneg;
	dmax__ = dsig;
	fneg = f0;
	if (abs(dsig) > abs(t1)) {
	    dsig = t1;
	    f0 = t2;
	}
    }

/*   Test for convergence. */

    stol = rtol * sig;
    if (abs(dmax__) <= stol || f >= 0. && f <= ftol || abs(f) <= rtol) {
	goto L6;
    }

/*   Test for F0 = F = FMAX or F < 0 on the first iteration. */

    if (f0 != f && (nit > 1 || f > 0.)) {
	goto L5;
    }

/*   F*F0 > 0 and either the new estimate would be outside of */
/*     the bracketing interval of length abs(DMAX) or F < 0 */
/*     on the first iteration.  Reset (SG0,F0) to (SNEG,FNEG). */

L4:
    dsig = dmax__;
    f0 = fneg;

/*   Compute the change in SIG by linear interpolation be- */
/*     tween (SG0,F0) and (SIG,F). */

L5:
    dsig = -f * dsig / (f - f0);
    if (abs(dsig) > abs(dmax__) || dsig * dmax__ > 0.f) {
	goto L4;
    }

/*   Restrict the step-size such that abs(DSIG) .GE. STOL/2. */
/*     Note that DSIG and DMAX have opposite signs. */

    if (abs(dsig) < stol / 2.) {
	d__1 = stol / 2.;
	dsig = -d_sign(&d__1, &dmax__);
    }

/*   Bottom of loop:  update SIG, DMAX, and F0. */

    sig += dsig;
    dmax__ += dsig;
    f0 = f;
    goto L2;

/* No errors encountered and SIGMA finite. */

L6:
    *ier = 0;
    ret_val = sig;
    return ret_val;

/* Infinite tension required. */

L7:
    *ier = 1;
    ret_val = sbig;
    return ret_val;

/* Error in an input parameter. */

L8:
    *ier = -1;
    ret_val = -1.;
    return ret_val;

/* Invalid constraint. */

L9:
    *ier = -2;
    ret_val = -1.;
    return ret_val;
} /* sig0_ */

doublereal sig1_(doublereal *x1, doublereal *x2, doublereal *y1, doublereal *
	y2, doublereal *y1p, doublereal *y2p, integer *ifl, doublereal *hpbnd,
	 doublereal *tol, integer *ier)
{
    /* Initialized data */

    static doublereal sbig = 85.;

    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double exp(doublereal), sqrt(doublereal), d_sign(doublereal *, doublereal 
	    *);

    /* Local variables */
    static doublereal d1pd2, fneg, dsig, dmax__, fmax, sinh__, ftol, rtol, 
	    stol, a, e, f, s, coshm, sinhm, a0, b0, c0, c1, c2, d0, d1, d2, 
	    f0;
    extern doublereal store_(doublereal *);
    static doublereal s1, s2, t0, t1, t2, rf, dx, tm, coshmm;
    extern /* Subroutine */ int snhcsh_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal bnd, sig, ems;
    static integer nit;
    static doublereal ems2;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   11/17/96 */

/*   Given a pair of abscissae with associated ordinates and */
/* slopes, this function determines the smallest (nonnega- */
/* tive) tension factor SIGMA such that the derivative HP(x) */
/* of the Hermite interpolatory tension spline H(x), defined */
/* by SIGMA and the data, is bounded (either above or below) */
/* by HPBND for all x in (X1,X2). */

/* On input: */

/*       X1,X2 = Abscissae.  X1 < X2. */

/*       Y1,Y2 = Values of H at X1 and X2. */

/*       Y1P,Y2P = Values of HP at X1 and X2. */

/*       IFL = Option indicator: */
/*             IFL = -1 if HPBND is a lower bound on HP. */
/*             IFL = 1 if HPBND is an upper bound on HP. */

/*       HPBND = Bound on HP.  If IFL = -1, HPBND .LE. */
/*               min(Y1P,Y2P,S) for S = (Y2-Y1)/(X2-X1).  If */
/*               IFL = 1, HPBND .GE. max(Y1P,Y2P,S). */

/*       TOL = Tolerance whose magnitude determines how close */
/*             SIGMA is to its optimal value when nonzero */
/*             finite tension is necessary and sufficient to */
/*             satisfy the constraint.  For a lower bound, */
/*             SIGMA is chosen so that HPBND .LE. HPMIN .LE. */
/*             HPBND + abs(TOL), where HPMIN is the minimum */
/*             value of HP in the interval, and for an upper */
/*             bound, the maximum of HP satisfies HPBND - */
/*             abs(TOL) .LE. HPMAX .LE. HPBND.  Thus, the */
/*             constraint is satisfied but possibly with more */
/*             tension than necessary. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered and the */
/*                     constraint can be satisfied with fin- */
/*                     ite tension. */
/*             IER = 1 if no errors were encountered but in- */
/*                     finite tension is required to satisfy */
/*                     the constraint (e.g., IFL = -1, HPBND */
/*                     = S, and Y1P > S). */
/*             IER = -1 if X2 .LE. X1 or abs(IFL) .NE. 1. */
/*             IER = -2 if HPBND is outside its valid range */
/*                      on input. */

/*       SIG1 = Minimum tension factor defined above unless */
/*              IER < 0, in which case SIG1 = -1.  If IER = */
/*              1, SIG1 = 85, resulting in an approximation */
/*              to the linear interpolant of the endpoint */
/*              values.  Note, however, that SIG1 may be */
/*              larger than 85 if IER = 0. */

/* Modules required by SIG1:  SNHCSH, STORE */

/* Intrinsic functions called by SIG1:  ABS, DBLE, EXP, MAX, */
/*                                        MIN, SIGN, SQRT */

/* *********************************************************** */



/* Store local parameters and test for errors. */

    rf = (doublereal) (*ifl);
    dx = *x2 - *x1;
    if (abs(rf) != 1. || dx <= 0.) {
	goto L7;
    }
    s1 = *y1p;
    s2 = *y2p;
    s = (*y2 - *y1) / dx;
    bnd = *hpbnd;

/* Test for a valid constraint. */

/* Computing MIN */
    d__1 = min(s1,s2);
/* Computing MAX */
    d__2 = max(s1,s2);
    if (rf < 0. && min(d__1,s) < bnd || rf > 0. && bnd < max(d__2,s)) {
	goto L8;
    }

/* Test for infinite tension required. */

    if (s == bnd && (s1 != s || s2 != s)) {
	goto L6;
    }

/* Test for SIG = 0 sufficient.  The Hermite cubic interpo- */
/*   land H0 has derivative HP0(x) = S2 + 2*B0*R + A0*R**2, */
/*   where R = (X2-x)/DX. */

    sig = 0.;
    t0 = s * 3. - s1 - s2;
    b0 = t0 - s2;
    c0 = t0 - s1;
    a0 = -b0 - c0;

/*   HP0(R) has an extremum (at R = -B0/A0) in (0,1) iff */
/*     B0*C0 > 0 and the third derivative of H0 has the */
/*     sign of A0. */

    if (b0 * c0 <= 0. || a0 * rf > 0.) {
	goto L5;
    }

/*   A0*RF < 0 and HP0(R) = -D0/A0 at R = -B0/A0. */

    d0 = t0 * t0 - s1 * s2;
    f0 = (bnd + d0 / a0) * rf;
    if (f0 >= 0.) {
	goto L5;
    }

/* Find a zero of F(SIG) = (BND-HP(R))*RF, where HP has an */
/*   extremum at R.  F has a unique zero, F(0) = F0 < 0, and */
/*   F = (BND-S)*RF > 0 for SIG sufficiently large. */

/* Initialize parameters for the secant method.  The method */
/*   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG), */
/*   where SG0 and SNEG are defined implicitly by DSIG = SIG */
/*   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and */
/*   SIG is initialized to the zero of (BND - (SIG*S-S1-S2)/ */
/*   (SIG-2.))*RF -- a value for which F(SIG) .GE. 0 and */
/*   F(SIG) = 0 for SIG sufficiently large that 2*SIG is in- */
/*   significant relative to EXP(SIG). */

    fmax = (bnd - s) * rf;
    sig = 2. - a0 / ((bnd - s) * 3.);
    d__1 = sig * exp(-sig) + .5;
    if (store_(&d__1) == .5) {
	goto L5;
    }
    dsig = sig;
    dmax__ = sig * -2.;
    fneg = fmax;
    d1 = s - s1;
    d2 = s2 - s;
    d1pd2 = d1 + d2;
    nit = 0;

/* Compute an absolute tolerance FTOL = abs(TOL), and a */
/*   relative tolerance RTOL = 100*MACHEPS. */

    ftol = abs(*tol);
    rtol = 1.;
L1:
    rtol /= 2.;
    d__1 = rtol + 1.;
    if (store_(&d__1) > 1.) {
	goto L1;
    }
    rtol *= 200.;

/* Top of loop:  compute F. */

L2:
    if (sig <= .5) {

/*   Use approximations designed to avoid cancellation error */
/*     (associated with small SIG) in the modified hyperbolic */
/*     functions. */

	snhcsh_(&sig, &sinhm, &coshm, &coshmm);
	c1 = sig * coshm * d2 - sinhm * d1pd2;
	c2 = sig * (sinhm + sig) * d2 - coshm * d1pd2;
	a = c2 - c1;
	e = sig * sinhm - coshmm - coshmm;
    } else {

/*   Scale SINHM and COSHM by 2*EXP(-SIG) in order to avoid */
/*     overflow. */

	ems = exp(-sig);
	ems2 = ems + ems;
	tm = 1. - ems;
	sinh__ = tm * (ems + 1.);
	sinhm = sinh__ - sig * ems2;
	coshm = tm * tm;
	c1 = sig * coshm * d2 - sinhm * d1pd2;
	c2 = sig * sinh__ * d2 - coshm * d1pd2;
	a = ems2 * (sig * tm * d2 + (tm - sig) * d1pd2);
	e = sig * sinh__ - coshm - coshm;
    }

/*   The second derivative of H(R) has a zero at EXP(SIG*R) = */
/*     SQRT((C2+C1)/A) and R is in (0,1) and well-defined */
/*     iff HPP(X1)*HPP(X2) < 0. */

    f = fmax;
    t1 = a * (c2 + c1);
    if (t1 >= 0.) {
	if (c1 * (sig * coshm * d1 - sinhm * d1pd2) < 0.) {

/*   HP(R) = (B+SIGN(A)*SQRT(A*C))/E at the critical value */
/*     of R, where A = C2-C1, B = E*S2-C2, and C = C2+C1. */
/*     NOTE THAT RF*A < 0. */

	    f = (bnd - (e * s2 - c2 - rf * sqrt(t1)) / e) * rf;
	}
    }

/*   Update the number of iterations NIT. */

    ++nit;
    if (f0 * f < 0.) {

/*   F0*F < 0.  Update (SNEG,FNEG) to (SG0,F0) so that F */
/*     and FNEG always have opposite signs.  If SIG is closer */
/*     to SNEG than SG0 and abs(F) < abs(FNEG), then swap */
/*     (SNEG,FNEG) with (SG0,F0). */

	t1 = dmax__;
	t2 = fneg;
	dmax__ = dsig;
	fneg = f0;
	if (abs(dsig) > abs(t1) && abs(f) < abs(t2)) {

	    dsig = t1;
	    f0 = t2;
	}
    }

/*   Test for convergence. */

    stol = rtol * sig;
    if (abs(dmax__) <= stol || f >= 0. && f <= ftol || abs(f) <= rtol) {
	goto L5;
    }
    if (f0 * f < 0. || abs(f) < abs(f0)) {
	goto L4;
    }

/*   F*F0 > 0 and the new estimate would be outside of the */
/*     bracketing interval of length abs(DMAX).  Reset */
/*     (SG0,F0) to (SNEG,FNEG). */

L3:
    dsig = dmax__;
    f0 = fneg;

/*   Compute the change in SIG by linear interpolation be- */
/*     tween (SG0,F0) and (SIG,F). */

L4:
    dsig = -f * dsig / (f - f0);
    if (abs(dsig) > abs(dmax__) || dsig * dmax__ > 0.f) {
	goto L3;
    }

/*   Restrict the step-size such that abs(DSIG) .GE. STOL/2. */
/*     Note that DSIG and DMAX have opposite signs. */

    if (abs(dsig) < stol / 2.) {
	d__1 = stol / 2.;
	dsig = -d_sign(&d__1, &dmax__);
    }

/*   Bottom of loop:  update SIG, DMAX, and F0. */

    sig += dsig;
    dmax__ += dsig;
    f0 = f;
    goto L2;

/* No errors encountered and SIGMA finite. */

L5:
    *ier = 0;
    ret_val = sig;
    return ret_val;

/* Infinite tension required. */

L6:
    *ier = 1;
    ret_val = sbig;
    return ret_val;

/* Error in an input parameter. */

L7:
    *ier = -1;
    ret_val = -1.;
    return ret_val;

/* Invalid constraint. */

L8:
    *ier = -2;
    ret_val = -1.;
    return ret_val;
} /* sig1_ */

doublereal sig2_(doublereal *x1, doublereal *x2, doublereal *y1, doublereal *
	y2, doublereal *y1p, doublereal *y2p, integer *ifl, doublereal *tol, 
	integer *ier)
{
    /* Initialized data */

    static doublereal sbig = 85.;

    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal);

    /* Local variables */
    static doublereal dsig, ftol, rtol, f, s, t, coshm, sinhm, d1, d2, dummy;
    extern doublereal store_(doublereal *);
    static doublereal t1, fp, dx;
    extern /* Subroutine */ int snhcsh_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal tp1, sig, ems;
    static integer nit;
    static doublereal ssm;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/08/92 */

/*   Given a pair of abscissae with associated ordinates and */
/* slopes, this function determines the smallest (nonnega- */
/* tive) tension factor SIGMA such that the Hermite interpo- */
/* latory tension spline H(x) preserves convexity (or con- */
/* cavity) of the data;  i.e., */

/*   Y1P .LE. S .LE. Y2P implies HPP(x) .GE. 0  or */
/*   Y1P .GE. S .GE. Y2P implies HPP(x) .LE. 0 */

/* for all x in the open interval (X1,X2), where S = (Y2-Y1)/ */
/* (X2-X1) and HPP denotes the second derivative of H.  Note, */
/* however, that infinite tension is required if Y1P = S or */
/* Y2P = S (unless Y1P = Y2P = S). */

/* On input: */

/*       X1,X2 = Abscissae.  X1 < X2. */

/*       Y1,Y2 = Values of H at X1 and X2. */

/*       Y1P,Y2P = Derivative values of H at X1 and X2. */

/*       IFL = Option indicator (sign of HPP): */
/*             IFL = -1 if HPP is to be bounded above by 0. */
/*             IFL = 1 if HPP is to be bounded below by 0 */
/*                     (preserve convexity of the data). */

/*       TOL = Tolerance whose magnitude determines how close */
/*             SIGMA is to its optimal value when nonzero */
/*             finite tension is necessary and sufficient to */
/*             satisfy convexity or concavity.  In the case */
/*             of convexity, SIGMA is chosen so that 0 .LE. */
/*             HPPMIN .LE. abs(TOL), where HPPMIN is the min- */
/*             imum value of HPP in the interval.  In the */
/*             case of concavity, the maximum value of HPP */
/*             satisfies -abs(TOL) .LE. HPPMAX .LE. 0.  Thus, */
/*             the constraint is satisfied but possibly with */
/*             more tension than necessary. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered and fin- */
/*                     ite tension is sufficient to satisfy */
/*                     the constraint. */
/*             IER = 1 if no errors were encountered but in- */
/*                     finite tension is required to satisfy */
/*                     the constraint. */
/*             IER = -1 if X2 .LE. X1 or abs(IFL) .NE. 1. */
/*             IER = -2 if the constraint cannot be satis- */
/*                      fied:  the sign of S-Y1P or Y2P-S */
/*                      does not agree with IFL. */

/*       SIG2 = Tension factor defined above unless IER < 0, */
/*              in which case SIG2 = -1.  If IER = 1, SIG2 */
/*              is set to 85, resulting in an approximation */
/*              to the linear interpolant of the endpoint */
/*              values.  Note, however, that SIG2 may be */
/*              larger than 85 if IER = 0. */

/* Modules required by SIG2:  SNHCSH, STORE */

/* Intrinsic functions called by SIG2:  ABS, EXP, MAX, MIN, */
/*                                        SQRT */

/* *********************************************************** */



/* Test for an errors in the input parameters. */

    dx = *x2 - *x1;
    if ((doublereal) abs(*ifl) != 1. || dx <= 0.) {
	goto L5;
    }

/* Compute the slope and second differences, and test for */
/*   an invalid constraint. */

    s = (*y2 - *y1) / dx;
    d1 = s - *y1p;
    d2 = *y2p - s;
    if ((doublereal) (*ifl) > 0. && min(d1,d2) < 0. || (doublereal) (*ifl) < 
	    0. && max(d1,d2) > 0.) {
	goto L6;
    }

/* Test for infinite tension required. */

    if (d1 * d2 == 0. && d1 != d2) {
	goto L4;
    }

/* Test for SIG = 0 sufficient. */

    sig = 0.;
    if (d1 * d2 == 0.) {
	goto L3;
    }
/* Computing MAX */
    d__1 = d1 / d2, d__2 = d2 / d1;
    t = max(d__1,d__2);
    if (t <= 2.) {
	goto L3;
    }

/* Find a zero of F(SIG) = SIG*COSHM(SIG)/SINHM(SIG) - (T+1). */
/*   Since the derivative of F vanishes at the origin, a */
/*   quadratic approximation is used to obtain an initial */
/*   estimate for the Newton method. */

    tp1 = t + 1.;
    sig = sqrt(t * 10. - 20.);
    nit = 0;

/*   Compute an absolute tolerance FTOL = abs(TOL) and a */
/*     relative tolerance RTOL = 100*MACHEPS. */

    ftol = abs(*tol);
    rtol = 1.;
L1:
    rtol /= 2.;
    d__1 = rtol + 1.;
    if (store_(&d__1) > 1.) {
	goto L1;
    }
    rtol *= 200.;

/* Evaluate F and its derivative FP. */

L2:
    if (sig <= .5) {

/*   Use approximations designed to avoid cancellation error */
/*     in the hyperbolic functions. */

	snhcsh_(&sig, &sinhm, &coshm, &dummy);
	t1 = coshm / sinhm;
	fp = t1 + sig * (sig / sinhm - t1 * t1 + 1.);
    } else {

/*   Scale SINHM and COSHM by 2*exp(-SIG) in order to avoid */
/*     overflow. */

	ems = exp(-sig);
	ssm = 1. - ems * (ems + sig + sig);
	t1 = (1. - ems) * (1. - ems) / ssm;
	fp = t1 + sig * (sig * 2. * ems / ssm - t1 * t1 + 1.);
    }

    f = sig * t1 - tp1;
    ++nit;

/*   Test for convergence. */

    if (fp <= 0.) {
	goto L3;
    }
    dsig = -f / fp;
    if (abs(dsig) <= rtol * sig || f >= 0. && f <= ftol || abs(f) <= rtol) {
	goto L3;
    }

/*   Update SIG. */

    sig += dsig;
    goto L2;

/* No errors encountered, and SIGMA is finite. */

L3:
    *ier = 0;
    ret_val = sig;
    return ret_val;

/* Infinite tension required. */

L4:
    *ier = 1;
    ret_val = sbig;
    return ret_val;

/* X2 .LE. X1 or abs(IFL) .NE. 1. */

L5:
    *ier = -1;
    ret_val = -1.;
    return ret_val;

/* The constraint cannot be satisfied. */

L6:
    *ier = -2;
    ret_val = -1.;
    return ret_val;
} /* sig2_ */

/* Subroutine */ int sigbi_(integer *n, doublereal *x, doublereal *y, 
	doublereal *yp, doublereal *tol, doublereal *b, doublereal *bmax, 
	doublereal *sigma, integer *icflg, doublereal *dsmax, integer *ier)
{
    /* Initialized data */

    static doublereal sbig = 85.;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer icfk;
    static doublereal dsig;
    static integer icnt, ierr, i__, k;
    static doublereal s, sigin;
    static integer nm1;
    static doublereal bnd;
    static integer ifl;
    static doublereal sig, dsm, bmx;
    extern doublereal sig0_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *), sig1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *), sig2_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/10/92 */

/*   Given a set of abscissae X with associated data values Y */
/* and derivatives YP, this subroutine determines the small- */
/* est (nonnegative) tension factors SIGMA such that the Her- */
/* mite interpolatory tension spline H(x) satisfies a set of */
/* user-specified constraints. */

/*   SIGBI may be used in conjunction with Subroutine YPC2 */
/* (or YPC2P) in order to produce a C-2 interpolant which */
/* satisfies the constraints.  This is achieved by calling */
/* YPC2 with SIGMA initialized to the zero vector, and then */
/* alternating calls to SIGBI with calls to YPC2 until the */
/* change in SIGMA is small (refer to the parameter descrip- */
/* tions for SIGMA, DSMAX and IER), or the maximum relative */
/* change in YP is bounded by a tolerance (a reasonable value */
/* is .01).  A similar procedure may be used to produce a C-2 */
/* shape-preserving smoothing curve (Subroutine SMCRV). */

/*   Refer to Subroutine SIGS for a means of selecting mini- */
/* mum tension factors to preserve shape properties of the */
/* data. */

/* On input: */

/*       N = Number of data points.  N .GE. 2. */

/*       X = Array of length N containing a strictly in- */
/*           creasing sequence of abscissae:  X(I) < X(I+1) */
/*           for I = 1,...,N-1. */

/*       Y = Array of length N containing data values (or */
/*           function values computed by SMCRV) associated */
/*           with the abscissae.  H(X(I)) = Y(I) for I = */
/*           1,...,N. */

/*       YP = Array of length N containing first derivatives */
/*            of H at the abscissae.  Refer to Subroutines */
/*            YPC1, YPC1P, YPC2, YPC2P, and SMCRV. */

/*       TOL = Tolerance whose magnitude determines how close */
/*             each tension factor is to its optimal value */
/*             when nonzero finite tension is necessary and */
/*             sufficient to satisfy a constraint.  Refer to */
/*             functions SIG0, SIG1, and SIG2.  TOL should be */
/*             set to 0 for optimal tension. */

/*       B = Array dimensioned 5 by N-1 containing bounds or */
/*           flags which define the constraints.  For I = 1 */
/*           to N-1, column I defines the constraints associ- */
/*           ated with interval I (X(I),X(I+1)) as follows: */

/*             B(1,I) is an upper bound on H */
/*             B(2,I) is a lower bound on H */
/*             B(3,I) is an upper bound on HP */
/*             B(4,I) is a lower bound on HP */
/*             B(5,I) specifies the required sign of HPP */

/*           where HP and HPP denote the first and second */
/*           derivatives of H, respectively.  A null con- */
/*           straint is specified by abs(B(K,I)) .GE. BMAX */
/*           for K < 5, or B(5,I) = 0:   B(1,I) .GE. BMAX, */
/*           B(2,I) .LE. -BMAX, B(3,I) .GE. BMAX, B(4,I) .LE. */
/*           -BMAX, or B(5,I) = 0.  Any positive value of */
/*           B(5,I) specifies that H should be convex, a */
/*           negative values specifies that H should be con- */
/*           cave, and 0 specifies that no restriction be */
/*           placed on HPP.  Refer to Functions SIG0, SIG1, */
/*           and SIG2 for definitions of valid constraints. */

/*       BMAX = User-defined value of infinity which, when */
/*              used as an upper bound in B (or when its */
/*              negative is used as a lower bound), specifies */
/*              that no constraint is to be enforced. */

/* The above parameters are not altered by this routine. */

/*       SIGMA = Array of length N-1 containing minimum val- */
/*               ues of the tension factors.  SIGMA(I) is as- */
/*               sociated with interval (I,I+1) and SIGMA(I) */
/*               .GE. 0 for I = 1,...,N-1.  SIGMA should be */
/*               set to the zero vector if minimal tension */
/*               is desired, and should be unchanged from a */
/*               previous call in order to ensure convergence */
/*               of the C-2 iterative procedure. */

/*       ICFLG = Array of length .GE. N-1. */

/* On output: */

/*       SIGMA = Array containing tension factors for which */
/*               H(x) satisfies the constraints defined by B, */
/*               with the restriction that SIGMA(I) .LE. 85 */
/*               for all I (unless the input value is larger). */
/*               The factors are as small as possible (within */
/*               the tolerance), but not less than their */
/*               input values.  If infinite tension is re- */
/*               quired in interval (X(I),X(I+1)), then */
/*               SIGMA(I) = 85 (and H is an approximation to */
/*               the linear interpolant on the interval), */
/*               and if no constraint is specified in the */
/*               interval, then SIGMA(I) = 0 (unless the */
/*               input value is positive), and thus H is */
/*               cubic.  Invalid constraints are treated as */
/*               null constraints. */

/*       ICFLG = Array of invalid constraint flags associated */
/*               with intervals.  For I = 1 to N-1, ICFLG(I) */
/*               is a 5-bit value b5b4b3b2b1, where bK = 1 if */
/*               and only if constraint K cannot be satis- */
/*               fied.  Thus, all constraints in interval I */
/*               are satisfied if and only if ICFLG(I) = 0 */
/*               (and IER .GE. 0). */

/*       DSMAX = Maximum increase in a component of SIGMA */
/*               from its input value.  The increase is a */
/*               relative change if the input value is */
/*               positive, and an absolute change otherwise. */

/*       IER = Error indicator and information flag: */
/*             IER = I if no errors (other than invalid con- */
/*                     straints) were encountered and I */
/*                     components of SIGMA were altered from */
/*                     their input values for 0 .LE. I .LE. */
/*                     N-1. */
/*             IER = -1 if N < 2.  SIGMA and ICFLG are not */
/*                      altered in this case. */
/*             IER = -I if X(I) .LE. X(I-1) for some I in the */
/*                      range 2,...,N.  SIGMA(J) and ICFLG(J) */
/*                      are unaltered for J .GE. I-1 in this */
/*                      case. */

/* Modules required by SIGBI:  SIG0, SIG1, SIG2, SNHCSH, */
/*                               STORE */

/* Intrinsic functions called by SIGBI:  ABS, MAX, MIN */

/* *********************************************************** */


    /* Parameter adjustments */
    --icflg;
    --sigma;
    b -= 6;
    --yp;
    --y;
    --x;

    /* Function Body */
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L4;
    }
    bmx = *bmax;

/* Initialize change counter ICNT and maximum change DSM for */
/*   loop on intervals. */

    icnt = 0;
    dsm = 0.;
    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (x[i__] >= x[i__ + 1]) {
	    goto L5;
	}
	icflg[i__] = 0;

/* Loop on constraints for interval I.  SIG is set to the */
/*   largest tension factor required to satisfy all five */
/*   constraints.  ICFK = 2**(K-1) is the increment for */
/*   ICFLG(I) when constraint K is invalid. */

	sig = 0.;
	icfk = 1;
	for (k = 1; k <= 5; ++k) {
	    bnd = b[k + i__ * 5];
	    if (k < 5 && abs(bnd) >= bmx) {
		goto L1;
	    }
	    if (k <= 2) {
		ifl = 3 - (k << 1);
		s = sig0_(&x[i__], &x[i__ + 1], &y[i__], &y[i__ + 1], &yp[i__]
			, &yp[i__ + 1], &ifl, &bnd, tol, &ierr);
	    } else if (k <= 4) {
		ifl = 7 - (k << 1);
		s = sig1_(&x[i__], &x[i__ + 1], &y[i__], &y[i__ + 1], &yp[i__]
			, &yp[i__ + 1], &ifl, &bnd, tol, &ierr);
	    } else {
		if (bnd == 0.) {
		    goto L1;
		}
		ifl = -1;
		if (bnd > 0.) {
		    ifl = 1;
		}
		s = sig2_(&x[i__], &x[i__ + 1], &y[i__], &y[i__ + 1], &yp[i__]
			, &yp[i__ + 1], &ifl, tol, &ierr);
	    }
	    if (ierr == -2) {

/*   An invalid constraint was encountered.  Increment */
/*     ICFLG(I). */

		icflg[i__] += icfk;
	    } else {

/*   Update SIG. */

		sig = max(sig,s);
	    }

/*   Bottom of loop on constraints K:  update ICFK. */

L1:
	    icfk <<= 1;
/* L2: */
	}

/* Bottom of loop on intervals:  update SIGMA(I), ICNT, and */
/*   DSM if necessary. */

	sig = min(sig,sbig);
	sigin = sigma[i__];
	if (sig > sigin) {
	    sigma[i__] = sig;
	    ++icnt;
	    dsig = sig - sigin;
	    if (sigin > 0.) {
		dsig /= sigin;
	    }
	    dsm = max(dsm,dsig);
	}
/* L3: */
    }

/* No errors (other than invalid constraints) encountered. */

    *dsmax = dsm;
    *ier = icnt;
    return 0;

/* N < 2. */

L4:
    *dsmax = 0.;
    *ier = -1;
    return 0;

/* X(I) .GE. X(I+1). */

L5:
    *dsmax = dsm;
    *ier = -(i__ + 1);
    return 0;
} /* sigbi_ */

/* Subroutine */ int sigbp_(integer *n, doublereal *x, doublereal *y, 
	doublereal *xp, doublereal *yp, doublereal *tol, doublereal *bl, 
	doublereal *bu, doublereal *bmax, doublereal *sigma, doublereal *
	dsmax, integer *ier)
{
    /* Initialized data */

    static doublereal sbig = 85.;

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *), exp(
	    doublereal), log(doublereal);

    /* Local variables */
    static doublereal fneg, dsig, dmax__, fmax, sneg;
    static integer icnt;
    static doublereal sinh__, ftol, rtol, stol, a, b, c__, d__, e, f;
    static integer i__;
    static doublereal s, t, coshm, sigin, sinhm, a1, a2, b0, d0, f0;
    extern doublereal store_(doublereal *);
    static doublereal t1, t2, v1, v2, aa, eb, dm, dp, dx, dy, rm, tm, rp, 
	    coshmm;
    extern /* Subroutine */ int snhcsh_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer ip1, nm1;
    static doublereal bhi, blo, sig, dsm, ems, v2m1, bmx;
    static integer nit;
    static doublereal rsm, rsp;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   11/18/96 */

/*   Given an ordered sequence of points C(I) = (X(I),Y(I)) */
/* with associated derivative vectors CP(I) = (XP(I),YP(I)), */
/* this subroutine determines the smallest (nonnegative) ten- */
/* sion factors SIGMA such that a parametric planar curve */
/* C(t) satisfies a set of user-specified constraints.  The */
/* components x(t) and y(t) of C(t) are the Hermite interpo- */
/* latory tension splines defined by the data and tension */
/* factors:  C(t(I)) = C(I) and C'(t(I)) = CP(I) for para- */
/* meter values t(1), t(2), ..., t(N).  In each subinterval */
/* [t1,t2], the signed perpendicular distance from the */
/* corresponding line segment C1-C2 to the curve C(t) is */
/* given by the vector cross product */

/*     d(t) = (C2-C1)/DC X (C(t)-C1) */

/* where DC = abs(C2-C1) is the length of the line segment. */
/* The associated tension factor SIGMA is chosen to satisfy */
/* an upper bound on the maximum of d(t) and a lower bound on */
/* the minimum of d(t) over t in [t1,t2].  Thus, the upper */
/* bound is associated with distance to the left of the line */
/* segment as viewed from C1 toward C2.  Note that the curve */
/* is assumed to be parameterized by arc length (Subroutine */
/* ARCL2D) so that t2-t1 = DC.  If this is not the case, the */
/* required bounds should be scaled by DC/(t2-t1) to obtain */
/* the input parameters BL and BU. */

/*   SIGBP may be used in conjunction with Subroutine YPC2 */
/* (or YPC2P) in order to produce a C-2 interpolant which */
/* satisfies the constraints.  This is achieved by calling */
/* YPC2 with SIGMA initialized to the zero vector, and then */
/* alternating calls to SIGBP with calls to YPC2 until the */
/* change in SIGMA is small (refer to the parameter descrip- */
/* tions for SIGMA, DSMAX and IER), or the maximum relative */
/* change in YP is bounded by a tolerance (a reasonable value */
/* is .01).  A similar procedure may be used to produce a C-2 */
/* shape-preserving smoothing curve (Subroutine SMCRV). */

/* On input: */

/*       N = Number of data points.  N .GE. 2. */

/*       X,Y = Arrays of length N containing the Cartesian */
/*             coordinates of the points C(I), I = 1 to N. */

/*       XP,YP = Arrays of length N containing the components */
/*               of the derivative (velocity) vectors CP(I). */
/*               Refer to Subroutines YPC1, YPC1P, YPC2, */
/*               YPC2P, and SMCRV. */

/*       TOL = Tolerance whose magnitude determines how close */
/*             each tension factor SIGMA is to its optimal */
/*             value when nonzero finite tension is necessary */
/*             and sufficient to satisfy a constraint. */
/*             SIGMA(I) is chosen so that BL(I) .LE. dmin */
/*             .LE. BL(I) + abs(TOL) and BU(I) - abs(TOL) */
/*             .LE. dmax .LE. BU(I), where dmin and dmax are */
/*             the minimum and maximum values of d(t) in the */
/*             interval [t(I),t(I+1)].  Thus, a large toler- */
/*             ance might increase execution efficiency but */
/*             may result in more tension than is necessary. */
/*             TOL may be set to 0 for optimal tension. */

/*       BL,BU = Arrays of length N-1 containing lower and */
/*               upper bounds, respectively, which define */
/*               the constraints as described above.  BL(I) */
/*               < 0 and BU(I) > 0 for I = 1 to N-1.  A null */
/*               straint is specified by BL(I) .LE. -BMAX or */
/*               BU(I) .GE. BMAX. */

/*       BMAX = User-defined value of infinity which, when */
/*              used as an upper bound in BU (or when its */
/*              negative is used as a lower bound in BL), */
/*              specifies that no constraint is to be en- */
/*              forced. */

/* The above parameters are not altered by this routine. */

/*       SIGMA = Array of length N-1 containing minimum val- */
/*               ues of the tension factors.  SIGMA(I) is as- */
/*               sociated with interval (I,I+1) and SIGMA(I) */
/*               .GE. 0 for I = 1,...,N-1.  SIGMA should be */
/*               set to the zero vector if minimal tension */
/*               is desired, and should be unchanged from a */
/*               previous call in order to ensure convergence */
/*               of the C-2 iterative procedure. */

/* On output: */

/*       SIGMA = Array containing tension factors for which */
/*               d(t) satisfies the constraints defined by */
/*               BL and BU, with the restriction that */
/*               SIGMA(I) .LE. 85 for all I (unless the input */
/*               value is larger).  The factors are as small */
/*               as possible (within the tolerance), but not */
/*               less than their input values.  If no con- */
/*               straint is specified in interval I, then */
/*               SIGMA(I) = 0 (unless the input value is */
/*               positive), and thus x(t) and y(t) are cubic */
/*               polynomials. */

/*       DSMAX = Maximum increase in a component of SIGMA */
/*               from its input value.  The increase is a */
/*               relative change if the input value is */
/*               positive, and an absolute change otherwise. */

/*       IER = Error indicator and information flag: */
/*             IER = I if no errors were encountered and I */
/*                     components of SIGMA were altered from */
/*                     their input values for 0 .LE. I .LE. */
/*                     N-1. */
/*             IER = -1 if N < 2.  SIGMA is not altered in */
/*                      this case. */
/*             IER = -I if BL(I-1) .GE. 0 or BU(I-1) .LE. 0 */
/*                      for some I in the range 2 to N. */
/*                      SIGMA(J) is unaltered for J .GE. I-1 */
/*                      in this case. */

/* Modules required by SIGBP:  SNHCSH, STORE */

/* Intrinsic functions called by SIGBP:  ABS, EXP, LOG, MAX, */
/*                                         MIN, SIGN, SQRT */

/* *********************************************************** */


    /* Parameter adjustments */
    --sigma;
    --bu;
    --bl;
    --yp;
    --xp;
    --y;
    --x;

    /* Function Body */
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L8;
    }
    bmx = *bmax;

/* Compute an absolute tolerance FTOL = abs(TOL), and a */
/*   relative tolerance RTOL = 100*MACHEPS. */

    ftol = abs(*tol);
    rtol = 1.;
L1:
    rtol /= 2.;
    d__1 = rtol + 1.;
    if (store_(&d__1) > 1.) {
	goto L1;
    }
    rtol *= 200.;

/* Initialize change counter ICNT and maximum change DSM for */
/*   loop on intervals. */

    icnt = 0;
    dsm = 0.;
    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ip1 = i__ + 1;
	blo = bl[i__];
	bhi = bu[i__];
	sigin = sigma[i__];
	if (blo >= 0. || bhi <= 0.) {
	    goto L9;
	}
	if (sigin >= sbig) {
	    goto L7;
	}

/* Initialize SIG to 0 and test for a null constraint. */

	sig = 0.;
	if (blo <= -bmx && bhi >= bmx) {
	    goto L6;
	}

/* Test for SIG = 0 sufficient. */

/*   The signed orthogonal distance is d(b) = b*(1-b)* */
/*     (b*V1 - (1-b)*V2), where b = (t2-t)/(t2-t1), */
/*     V1 = (C2-C1) X CP(1), and V2 = (C2-C1) X CP(2). */

	dx = x[ip1] - x[i__];
	dy = y[ip1] - y[i__];
	v1 = dx * yp[i__] - dy * xp[i__];
	v2 = dx * yp[ip1] - dy * xp[ip1];

/*   Set DP and DM to the maximum and minimum values of d(b) */
/*     for b in [0,1].  Note that DP .GE. 0 and DM .LE. 0. */

	s = v1 + v2;
	if (s == 0.) {

/*   The derivative d'(b) is zero at the midpoint b = .5. */

	    if (v1 >= 0.) {
		dp = v1 / 4.;
		dm = 0.;
	    } else {
		dp = 0.;
		dm = v1 / 4.;
	    }
	} else {

/*   Set RP/RM to the roots of the quadratic equation d'(b) = */
/*     (B0 +/- SQRT(D0))/(3*S) = V2/(B0 -/+ SQRT(D0)) = 0, */
/*     where B0 = V1 + 2*V2 and D0 = V1**2 + V1*V2 + V2**2. */
/*     The expression is chosen to avoid cancellation error. */

	    b0 = s + v2;
	    d0 = s * s - v1 * v2;
	    d__1 = sqrt(d0);
	    t = b0 + d_sign(&d__1, &b0);
	    if (b0 >= 0.) {
		rp = t / (s * 3.);
		rm = v2 / t;
	    } else {
		rp = v2 / t;
		rm = t / (s * 3.);
	    }
	    if (v1 <= 0. && v2 >= 0.) {

/*   The maximum is DP = 0 at the endpoints. */

		dp = 0.;
	    } else {
		dp = rp * (1. - rp) * (rp * s - v2);
	    }
	    if (v1 >= 0. && v2 <= 0.) {

/*   The minimum is DM = 0 at the endpoints. */

		dm = 0.;
	    } else {
		dm = rm * (1. - rm) * (rm * s - v2);
	    }
	}

/*   SIG = 0 is sufficient to satisfy the constraints iff */
/*     DP .LE. BHI and DM .GE. BLO iff F0 .GE. 0. */

/* Computing MIN */
	d__1 = bhi - dp, d__2 = dm - blo;
	f0 = min(d__1,d__2);
	if (f0 >= 0.) {
	    goto L6;
	}

/* Find a zero of F(SIG) = min(BHI-DP,DM-BLO), where DP and */
/*   DM are the maximum and minimum values of d(b).  F is an */
/*   increasing function, F(0) = F0 < 0, and F = FMAX = */
/*   min(BHI,-BLO) for SIG sufficiently large.  Note that F */
/*   has a discontinuity in its first derivative if the */
/*   curves BHI-DP and DM-BLO (as functions of SIG) inter- */
/*   sect, and the rate of convergence of the zero finder is */
/*   reduced to linear if such an intersection occurs near */
/*   the zero of F. */

/* Initialize parameters for the secant method.  The method */
/*   uses three points:  (SG0,F0), (SIG,F), and (SNEG,FNEG), */
/*   where SG0 and SNEG are defined implicitly by DSIG = SIG */
/*   - SG0 and DMAX = SIG - SNEG.  SG0 is initially zero and */
/*   SNEG is initialized to a sufficiently large value that */
/*   FNEG > 0.  This value is used only if the initial value */
/*   of F is negative. */

/* Computing MIN */
	d__1 = bhi, d__2 = -blo;
	t = min(d__1,d__2);
	fmax = max(.001,t);
/* Computing MAX */
	d__1 = abs(v1), d__2 = abs(v2);
	sig = max(d__1,d__2) / t;
	dmax__ = sig * (1. - t / fmax);
	sneg = sig - dmax__;
	dsig = sig;
	fneg = fmax;
	v2m1 = v2 - v1;
	nit = 0;

/* Top of loop:  compute F. */

L2:
	ems = exp(-sig);
	if (sig <= .5) {

/*   SIG .LE. .5:  use approximations designed to avoid can- */
/*                 cellation error (associated with small */
/*                 SIG) in the modified hyperbolic functions. */

	    snhcsh_(&sig, &sinhm, &coshm, &coshmm);
	    sinh__ = sinhm + sig;
	    a1 = sig * coshm * v2 - sinhm * v2m1;
	    a2 = sig * sinh__ * v2 - coshm * v2m1;
	    a = a2 - a1;
	    aa = a / ems;
	    e = sig * sinhm - coshmm - coshmm;
	} else {

/*   SIG > .5:  scale SINHM and COSHM by 2*EXP(-SIG) in order */
/*              to avoid overflow. */

	    tm = 1. - ems;
	    sinh__ = tm * (ems + 1.);
	    sinhm = sinh__ - sig * 2. * ems;
	    coshm = tm * tm;
	    a1 = sig * coshm * v2 - sinhm * v2m1;
	    a2 = sig * sinh__ * v2 - coshm * v2m1;
	    aa = (sig * tm * v2 + (tm - sig) * v2m1) * 2.;
	    a = ems * aa;
	    e = sig * sinh__ - coshm - coshm;
	}
	if (s == 0.) {

/*   The derivative d'(b) is zero at the midpoint b = .5. */

	    eb = sig * coshm - sinhm - sinhm;
	    if (v1 >= 0.) {
		dp = e * v1 / (sig * (sqrt(eb * eb - e * e) + eb));
		dm = 0.;
	    } else {
		dp = 0.;
		dm = e * v1 / (sig * (sqrt(eb * eb - e * e) + eb));
	    }
	} else {

/*   d'(b)*DC = V2 - (A1*sinh(SIG*b) - A2*coshm(SIG*b))/E = 0 */
/*     for ESB = (-B +/- sqrt(D))/A = C/(-B -/+ sqrt(D)), */
/*     where ESB = exp(SIG*b), A = A2-A1, D = B**2 - A*C, and */
/*     B and C are defined below. */

	    b = -coshm * s;
	    c__ = a2 + a1;
	    d__ = b * b - a * c__;
	    f = fmax;
	    if (d__ < 0.) {
		goto L3;
	    }
	    t1 = sqrt(d__);
	    t = -b - d_sign(&t1, &b);

	    rsp = 0.;
	    if (b < 0. && aa != 0.) {
		if (t / aa > 0.) {
		    rsp = sig + log(t / aa);
		}
	    }
	    if ((b > 0. || aa == 0.) && c__ / t > 0.) {
		rsp = log(c__ / t);
	    }
	    if ((rsp <= 0. || rsp >= sig) && b != 0.) {

/*   The maximum is DP = 0 at the endpoints. */

		dp = 0.;
	    } else {
		dp = -(b * rsp + a1 + t1) / (sig * e);
	    }

	    rsm = 0.;
	    if (b > 0. && aa != 0.) {
		if (t / aa > 0.) {
		    rsm = sig + log(t / aa);
		}
	    }
	    if ((b < 0. || aa == 0.) && c__ / t > 0.) {
		rsm = log(c__ / t);
	    }
	    if ((rsm <= 0. || rsm >= sig) && b != 0.) {

/*   The minimum is DM = 0 at the endpoints. */

		dm = 0.;
	    } else {
		dm = -(b * rsm + a1 - t1) / (sig * e);
	    }
	}

/* Computing MIN */
	d__1 = bhi - dp, d__2 = dm - blo;
	f = min(d__1,d__2);

/*   Update the number of iterations NIT. */

L3:
	++nit;
	if (f0 * f < 0.) {

/*   F0*F < 0.  Update (SNEG,FNEG) to (SG0,F0) so that F */
/*     and FNEG always have opposite signs.  If SIG is */
/*     closer to SNEG than SG0, then swap (SNEG,FNEG) with */
/*     (SG0,F0). */

	    t1 = dmax__;
	    t2 = fneg;
	    dmax__ = dsig;
	    fneg = f0;
	    if (abs(dsig) > abs(t1)) {
		dsig = t1;
		f0 = t2;
	    }
	}

/*   Test for convergence. */

	stol = rtol * sig;
	if (abs(dmax__) <= stol || f >= 0. && f <= ftol || abs(f) <= rtol) {
	    goto L6;
	}

/*   Test for F0 = F = FMAX or F < 0 on the first iteration. */

	if (f0 != f && (nit > 1 || f > 0.)) {
	    goto L5;
	}

/*   F*F0 > 0 and either the new estimate would be outside of */
/*     the bracketing interval of length abs(DMAX) or F < 0 */
/*     on the first iteration.  Reset (SG0,F0) to (SNEG,FNEG). */

L4:
	dsig = dmax__;
	f0 = fneg;

/*   Compute the change in SIG by linear interpolation be- */
/*     tween (SG0,F0) and (SIG,F). */

L5:
	dsig = -f * dsig / (f - f0);
	if (abs(dsig) > abs(dmax__) || dsig * dmax__ > 0.f) {
	    goto L4;
	}

/*   Restrict the step-size such that abs(DSIG) .GE. STOL/2. */
/*     Note that DSIG and DMAX have opposite signs. */

	if (abs(dsig) < stol / 2.) {
	    d__1 = stol / 2.;
	    dsig = -d_sign(&d__1, &dmax__);
	}

/*   Bottom of loop:  update SIG, DMAX, and F0. */

	sig += dsig;
	dmax__ += dsig;
	f0 = f;
	goto L2;

/* Bottom of loop on intervals:  update SIGMA(I), ICNT, and */
/*   DSM if necessary. */

L6:
	sig = min(sig,sbig);
	if (sig > sigin) {
	    sigma[i__] = sig;
	    ++icnt;
	    dsig = sig - sigin;
	    if (sigin > 0.) {
		dsig /= sigin;
	    }
	    dsm = max(dsm,dsig);
	}
L7:
	;
    }

/* No errors encountered. */

    *dsmax = dsm;
    *ier = icnt;
    return 0;

/* N < 2. */

L8:
    *dsmax = 0.;
    *ier = -1;
    return 0;

/* BL(I) .GE. 0 or BU(I) .LE. 0. */

L9:
    *dsmax = dsm;
    *ier = -ip1;
    return 0;
} /* sigbp_ */

/* Subroutine */ int sigs_(integer *n, doublereal *x, doublereal *y, 
	doublereal *yp, doublereal *tol, doublereal *sigma, doublereal *dsmax,
	 integer *ier)
{
    /* Initialized data */

    static doublereal sbig = 85.;

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal), d_sign(doublereal *, doublereal 
	    *);

    /* Local variables */
    static doublereal d1pd2, fneg, dsig, dmax__, fmax;
    static integer icnt;
    static doublereal ftol, rtol, stol, a, e, f;
    static integer i__;
    static doublereal s, t, coshm, sigin, sinhm, c1, c2, d0, d1, d2, f0, 
	    ssinh;
    extern doublereal store_(doublereal *);
    static doublereal s1, s2, t0, t1, t2, fp, dx, tm, coshmm;
    extern /* Subroutine */ int snhcsh_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer ip1, nm1;
    static doublereal tp1, d1d2, scm, dsm, ems, sig, sgn;
    static integer nit;
    static doublereal ssm, ems2;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   11/17/96 */

/*   Given a set of abscissae X with associated data values Y */
/* and derivatives YP, this subroutine determines the small- */
/* est (nonnegative) tension factors SIGMA such that the Her- */
/* mite interpolatory tension spline H(x) preserves local */
/* shape properties of the data.  In an interval (X1,X2) with */
/* data values Y1,Y2 and derivatives YP1,YP2, the properties */
/* of the data are */

/*       Monotonicity:  S, YP1, and YP2 are nonnegative or */
/*                        nonpositive, */
/*  and */
/*       Convexity:     YP1 .LE. S .LE. YP2  or  YP1 .GE. S */
/*                        .GE. YP2, */

/* where S = (Y2-Y1)/(X2-X1).  The corresponding properties */
/* of H are constant sign of the first and second deriva- */
/* tives, respectively.  Note that, unless YP1 = S = YP2, in- */
/* finite tension is required (and H is linear on the inter- */
/* val) if S = 0 in the case of monotonicity, or if YP1 = S */
/* or YP2 = S in the case of convexity. */

/*   SIGS may be used in conjunction with Subroutine YPC2 */
/* (or YPC2P) in order to produce a C-2 interpolant which */
/* preserves the shape properties of the data.  This is */
/* achieved by calling YPC2 with SIGMA initialized to the */
/* zero vector, and then alternating calls to SIGS with */
/* calls to YPC2 until the change in SIGMA is small (refer to */
/* the parameter descriptions for SIGMA, DSMAX and IER), or */
/* the maximum relative change in YP is bounded by a toler- */
/* ance (a reasonable value is .01).  A similar procedure may */
/* be used to produce a C-2 shape-preserving smoothing curve */
/* (Subroutine SMCRV). */

/*   Refer to Subroutine SIGBI for a means of selecting mini- */
/* mum tension factors to satisfy more general constraints. */

/* On input: */

/*       N = Number of data points.  N .GE. 2. */

/*       X = Array of length N containing a strictly in- */
/*           creasing sequence of abscissae:  X(I) < X(I+1) */
/*           for I = 1,...,N-1. */

/*       Y = Array of length N containing data values (or */
/*           function values computed by SMCRV) associated */
/*           with the abscissae.  H(X(I)) = Y(I) for I = */
/*           1,...,N. */

/*       YP = Array of length N containing first derivatives */
/*            of H at the abscissae.  Refer to Subroutines */
/*            YPC1, YPC1P, YPC2, YPC2P, and SMCRV. */

/*       TOL = Tolerance whose magnitude determines how close */
/*             each tension factor is to its optimal value */
/*             when nonzero finite tension is necessary and */
/*             sufficient to satisfy the constraint: */
/*             abs(TOL) is an upper bound on the magnitude */
/*             of the smallest (nonnegative) or largest (non- */
/*             positive) value of the first or second deriva- */
/*             tive of H in the interval.  Thus, the con- */
/*             straint is satisfied, but possibly with more */
/*             tension than necessary.  TOL should be set to */
/*             0 for optimal tension. */

/* The above parameters are not altered by this routine. */

/*       SIGMA = Array of length N-1 containing minimum val- */
/*               ues of the tension factors.  SIGMA(I) is as- */
/*               sociated with interval (I,I+1) and SIGMA(I) */
/*               .GE. 0 for I = 1,...,N-1.  SIGMA should be */
/*               set to the zero vector if minimal tension */
/*               is desired, and should be unchanged from a */
/*               previous call in order to ensure convergence */
/*               of the C-2 iterative procedure. */

/* On output: */

/*       SIGMA = Array containing tension factors for which */
/*               H(x) preserves the properties of the data, */
/*               with the restriction that SIGMA(I) .LE. 85 */
/*               for all I (unless the input value is larger). */
/*               The factors are as small as possible (within */
/*               the tolerance), but not less than their */
/*               input values.  If infinite tension is re- */
/*               quired in interval (X(I),X(I+1)), then */
/*               SIGMA(I) = 85 (and H is an approximation to */
/*               the linear interpolant on the interval), */
/*               and if neither property is satisfied by the */
/*               data, then SIGMA(I) = 0 (unless the input */
/*               value is positive), and thus H is cubic in */
/*               the interval. */

/*       DSMAX = Maximum increase in a component of SIGMA */
/*               from its input value.  The increase is a */
/*               relative change if the input value is */
/*               nonzero, and an absolute change otherwise. */

/*       IER = Error indicator and information flag: */
/*             IER = I if no errors were encountered and I */
/*                     components of SIGMA were altered from */
/*                     their input values for 0 .LE. I .LE. */
/*                     N-1. */
/*             IER = -1 if N < 2.  SIGMA is not altered in */
/*                      this case. */
/*             IER = -I if X(I) .LE. X(I-1) for some I in the */
/*                      range 2,...,N.  SIGMA(J-1) is unal- */
/*                      tered for J = I,...,N in this case. */

/* Modules required by SIGS:  SNHCSH, STORE */

/* Intrinsic functions called by SIGS:  ABS, EXP, MAX, MIN, */
/*                                        SIGN, SQRT */

/* *********************************************************** */


    /* Parameter adjustments */
    --sigma;
    --yp;
    --y;
    --x;

    /* Function Body */
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L9;
    }

/* Compute an absolute tolerance FTOL = abs(TOL) and a */
/*   relative tolerance RTOL = 100*MACHEPS. */

    ftol = abs(*tol);
    rtol = 1.;
L1:
    rtol /= 2.;
    d__1 = rtol + 1.;
    if (store_(&d__1) > 1.) {
	goto L1;
    }
    rtol *= 200.;

/* Initialize change counter ICNT and maximum change DSM for */
/*   loop on intervals. */

    icnt = 0;
    dsm = 0.;
    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ip1 = i__ + 1;
	dx = x[ip1] - x[i__];
	if (dx <= 0.) {
	    goto L10;
	}
	sigin = sigma[i__];
	if (sigin >= sbig) {
	    goto L8;
	}

/* Compute first and second differences. */

	s1 = yp[i__];
	s2 = yp[ip1];
	s = (y[ip1] - y[i__]) / dx;
	d1 = s - s1;
	d2 = s2 - s;
	d1d2 = d1 * d2;

/* Test for infinite tension required to satisfy either */
/*   property. */

	sig = sbig;
	if (d1d2 == 0. && s1 != s2 || s == 0. && s1 * s2 > 0.) {
	    goto L7;
	}

/* Test for SIGMA = 0 sufficient.  The data satisfies convex- */
/*   ity iff D1D2 .GE. 0, and D1D2 = 0 implies S1 = S = S2. */

	sig = 0.;
	if (d1d2 < 0.) {
	    goto L3;
	}
	if (d1d2 == 0.) {
	    goto L7;
	}
/* Computing MAX */
	d__1 = d1 / d2, d__2 = d2 / d1;
	t = max(d__1,d__2);
	if (t <= 2.) {
	    goto L7;
	}
	tp1 = t + 1.;

/* Convexity:  Find a zero of F(SIG) = SIG*COSHM(SIG)/ */
/*   SINHM(SIG) - TP1. */

/*   F(0) = 2-T < 0, F(TP1) .GE. 0, the derivative of F */
/*     vanishes at SIG = 0, and the second derivative of F is */
/*     .2 at SIG = 0.  A quadratic approximation is used to */
/*     obtain a starting point for the Newton method. */

	sig = sqrt(t * 10. - 20.);
	nit = 0;

/*   Top of loop: */

L2:
	if (sig <= .5) {
	    snhcsh_(&sig, &sinhm, &coshm, &coshmm);
	    t1 = coshm / sinhm;
	    fp = t1 + sig * (sig / sinhm - t1 * t1 + 1.);
	} else {

/*   Scale SINHM and COSHM by 2*EXP(-SIG) in order to avoid */
/*     overflow with large SIG. */

	    ems = exp(-sig);
	    ssm = 1. - ems * (ems + sig + sig);
	    t1 = (1. - ems) * (1. - ems) / ssm;
	    fp = t1 + sig * (sig * 2. * ems / ssm - t1 * t1 + 1.);
	}

	f = sig * t1 - tp1;
	++nit;

/*   Test for convergence. */

	if (fp <= 0.) {
	    goto L7;
	}
	dsig = -f / fp;
	if (abs(dsig) <= rtol * sig || f >= 0. && f <= ftol || abs(f) <= rtol)
		 {
	    goto L7;
	}

/*   Update SIG. */

	sig += dsig;
	goto L2;

/* Convexity cannot be satisfied.  Monotonicity can be satis- */
/*   fied iff S1*S .GE. 0 and S2*S .GE. 0 since S .NE. 0. */

L3:
	if (s1 * s < 0. || s2 * s < 0.) {
	    goto L7;
	}
	t0 = s * 3. - s1 - s2;
	d0 = t0 * t0 - s1 * s2;

/* SIGMA = 0 is sufficient for monotonicity iff S*T0 .GE. 0 */
/*   or D0 .LE. 0. */

	if (d0 <= 0. || s * t0 >= 0.) {
	    goto L7;
	}

/* Monotonicity:  find a zero of F(SIG) = SIGN(S)*HP(R), */
/*   where HPP(R) = 0 and HP, HPP denote derivatives of H. */
/*   F has a unique zero, F(0) < 0, and F approaches abs(S) */
/*   as SIG increases. */

/*   Initialize parameters for the secant method.  The method */
/*     uses three points:  (SG0,F0), (SIG,F), and */
/*     (SNEG,FNEG), where SG0 and SNEG are defined implicitly */
/*     by DSIG = SIG - SG0 and DMAX = SIG - SNEG. */

	sgn = d_sign(&c_b82, &s);
	sig = sbig;
	fmax = sgn * (sig * s - s1 - s2) / (sig - 2.);
	if (fmax <= 0.) {
	    goto L7;
	}
	stol = rtol * sig;
	f = fmax;
	f0 = sgn * d0 / ((d1 - d2) * 3.);
	fneg = f0;
	dsig = sig;
	dmax__ = sig;
	d1pd2 = d1 + d2;
	nit = 0;

/*   Top of loop:  compute the change in SIG by linear */
/*     interpolation. */

L4:
	dsig = -f * dsig / (f - f0);
	if (abs(dsig) > abs(dmax__) || dsig * dmax__ > 0.f) {
	    goto L6;
	}

/*   Restrict the step-size such that abs(DSIG) .GE. STOL/2. */
/*     Note that DSIG and DMAX have opposite signs. */

	if (abs(dsig) < stol / 2.) {
	    d__1 = stol / 2.;
	    dsig = -d_sign(&d__1, &dmax__);
	}

/*   Update SIG, F0, and F. */

	sig += dsig;
	f0 = f;
	if (sig <= .5) {

/*   Use approximations to the hyperbolic functions designed */
/*     to avoid cancellation error with small SIG. */

	    snhcsh_(&sig, &sinhm, &coshm, &coshmm);
	    c1 = sig * coshm * d2 - sinhm * d1pd2;
	    c2 = sig * (sinhm + sig) * d2 - coshm * d1pd2;
	    a = c2 - c1;
	    e = sig * sinhm - coshmm - coshmm;
	} else {

/*   Scale SINHM and COSHM by 2*EXP(-SIG) in order to avoid */
/*     overflow with large SIG. */

	    ems = exp(-sig);
	    ems2 = ems + ems;
	    tm = 1. - ems;
	    ssinh = tm * (ems + 1.);
	    ssm = ssinh - sig * ems2;
	    scm = tm * tm;
	    c1 = sig * scm * d2 - ssm * d1pd2;
	    c2 = sig * ssinh * d2 - scm * d1pd2;

/*   R is in (0,1) and well-defined iff HPP(X1)*HPP(X2) < 0. */

	    f = fmax;
	    if (c1 * (sig * scm * d1 - ssm * d1pd2) >= 0.) {
		goto L5;
	    }
	    a = ems2 * (sig * tm * d2 + (tm - sig) * d1pd2);
	    if (a * (c2 + c1) < 0.) {
		goto L5;
	    }
	    e = sig * ssinh - scm - scm;
	}

	f = (sgn * (e * s2 - c2) + sqrt(a * (c2 + c1))) / e;

/*   Update number of iterations NIT. */

L5:
	++nit;

/*   Test for convergence. */

	stol = rtol * sig;
	if (abs(dmax__) <= stol || f >= 0. && f <= ftol || abs(f) <= rtol) {
	    goto L7;
	}
	dmax__ += dsig;
	if (f0 * f > 0. && abs(f) >= abs(f0)) {
	    goto L6;
	}
	if (f0 * f <= 0.) {

/*   F and F0 have opposite signs.  Update (SNEG,FNEG) to */
/*     (SG0,F0) so that F and FNEG always have opposite */
/*     signs.  If SIG is closer to SNEG than SG0 and abs(F) < */
/*     abs(FNEG), then swap (SNEG,FNEG) with (SG0,F0). */

	    t1 = dmax__;
	    t2 = fneg;
	    dmax__ = dsig;
	    fneg = f0;
	    if (abs(dsig) > abs(t1) && abs(f) < abs(t2)) {

		dsig = t1;
		f0 = t2;
	    }
	}
	goto L4;

/*   Bottom of loop:  F0*F > 0 and the new estimate would */
/*     be outside of the bracketing interval of length */
/*     abs(DMAX).  Reset (SG0,F0) to (SNEG,FNEG). */

L6:
	dsig = dmax__;
	f0 = fneg;
	goto L4;

/*  Update SIGMA(I), ICNT, and DSM if necessary. */

L7:
	sig = min(sig,sbig);
	if (sig > sigin) {
	    sigma[i__] = sig;
	    ++icnt;
	    dsig = sig - sigin;
	    if (sigin > 0.) {
		dsig /= sigin;
	    }
	    dsm = max(dsm,dsig);
	}
L8:
	;
    }

/* No errors encountered. */

    *dsmax = dsm;
    *ier = icnt;
    return 0;

/* N < 2. */

L9:
    *dsmax = 0.;
    *ier = -1;
    return 0;

/* X(I+1) .LE. X(I). */

L10:
    *dsmax = dsm;
    *ier = -ip1;
    return 0;
} /* sigs_ */

/* Subroutine */ int smcrv_(integer *n, doublereal *x, doublereal *y, 
	doublereal *sigma, logical *period, doublereal *w, doublereal *sm, 
	doublereal *smtol, doublereal *wk, doublereal *ys, doublereal *yp, 
	integer *ier)
{
    /* System generated locals */
    integer wk_dim1, wk_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal gneg, dmax__;
    static integer ierr, iter;
    static doublereal wixi;
    extern /* Subroutine */ int b2tri_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static doublereal q2min, q2max, d__, g;
    static integer i__;
    static doublereal p, s, g0, h0, p0, q2, r1, r2;
    extern /* Subroutine */ int b2trip_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static doublereal c11, c12, c22, dp, sd, hp;
    static integer nn;
    static doublereal dx, wi, xi, yi;
    extern /* Subroutine */ int ypcoef_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer nm1;
    static doublereal sig;
    static logical per;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   10/05/98 */

/*   Given a sequence of abscissae X with associated data */
/* values Y and tension factors SIGMA, this routine deter- */
/* mines a set of function values YS and first derivatives YP */
/* associated with a Hermite interpolatory tension spline */
/* H(x) which smoothes the data.  H(x) has two continuous */
/* derivatives for all x and satisfies either natural or per- */
/* iodic end conditions.  The values and derivatives are */
/* chosen to minimize a quadratic functional Q1(YS,YP) sub- */
/* ject to the constraint Q2(YS) .LE. SM for Q2(YS) = */
/* (Y-YS)**T*W*(Y-YS), where **T denotes transpose and W is a */
/* diagonal matrix of positive weights. */

/*   Functions HVAL, HPVAL, HPPVAL, and TSINTL may be called */
/* to compute values, derivatives, and integrals of H.  The */
/* function values YS must be used as data values in those */
/* subprograms. */

/*   The smoothing procedure is an extension of the method */
/* for cubic spline smoothing due to C. Reinsch:  Numer. */
/* Math., 10 (1967) and 16 (1971).  Q1 is defined as the sum */
/* of integrals over the intervals (X(I),X(I+1)) of HPP**2 + */
/* (SIGMA(I)/DX)**2*(HP-S)**2, where DX = X(I+1)-X(I), HP and */
/* HPP denote first and second derivatives of H, and S = */
/* (YS(I+1)-YS(I))/DX.  Introducing a smoothing parameter P, */
/* and assuming the constraint is active, the problem is */
/* equivalent to minimizing Q(P,YS,YP) = Q1(YS,YP) + */
/* P*(Q2(YS)-SM).  The secant method is used to find a zero */
/* of G(P) = 1/SQRT(Q2) - 1/SQRT(SM), where YS and YP satisfy */
/* the order 2N symmetric positive-definite linear system */
/* obtained by setting the gradient of Q (treated as a func- */
/* tion of YS and YP) to zero. */

/*   Note that the interpolation problem corresponding to */
/* YS = Y, SM = 0, and P infinite is solved by Subroutine */
/* YPC2 or YPC2P. */

/* On input: */

/*       N = Number of data points.  N .GE. 2 if PERIOD = */
/*           FALSE, and N .GE. 3 if PERIOD = TRUE. */

/*       X = Array of length N containing a strictly in- */
/*           creasing sequence of abscissae:  X(I) < X(I+1) */
/*           for I = 1,...,N-1. */

/*       Y = Array of length N containing data values assoc- */
/*           iated with the abscissae.  If PERIOD = TRUE, it */
/*           is assumed that Y(N) = Y(1). */

/*       SIGMA = Array of length N-1 containing tension */
/*               factors.  SIGMA(I) is associated with inter- */
/*               val (X(I),X(I+1)) for I = 1,...,N-1.  If */
/*               SIGMA(I) = 0, H is cubic, and as SIGMA in- */
/*               creases, H approaches linear in the inter- */
/*               val. */

/*       PERIOD = Periodic end condition flag: */
/*                PERIOD = .F. if H is to satisfy natural end */
/*                             conditions:  zero second der- */
/*                             ivatives at X(1) and X(N). */
/*                PERIOD = .T. if H is to satisfy periodic */
/*                             end conditions:  the values */
/*                             and first two derivatives at */
/*                             X(1) agree with those at X(N), */
/*                             and a period thus has length */
/*                             X(N)-X(1). */

/*       W = Array of length N containing positive weights */
/*           associated with the data values.  The recommend- */
/*           ed value of W(I) is 1/DY**2, where DY is the */
/*           standard deviation associated with Y(I).  If */
/*           nothing is known about the errors in Y, a con- */
/*           stant (estimated value) should be used for DY. */
/*           If PERIOD = TRUE, it is assumed that W(N) = */
/*           W(1). */

/*       SM = Positive parameter specifying an upper bound on */
/*            Q2(YS).  H(x) is linear (and Q2 is minimized) */
/*            if SM is sufficiently large that the constraint */
/*            is not active.  It is recommended that SM sat- */
/*            isfy N-SQRT(2N) .LE. SM .LE. N+SQRT(2N) and */
/*            SM = N is reasonable if W(I) = 1/DY**2. */

/*       SMTOL = Parameter in the range (0,1) specifying the */
/*               relative error allowed in satisfying the */
/*               constraint:  the constraint is assumed to */
/*               be satisfied if SM*(1-SMTOL) .LE. Q2 .LE. */
/*               SM*(1+SMTOL).  A reasonable value for SMTOL */
/*               is SQRT(2/N). */

/* The above parameters are not altered by this routine. */

/*       WK = Work space of length at least 6N if PERIOD = */
/*            FALSE, and 10N if PERIOD = TRUE. */

/* On output: */

/*       YS = Array of length N containing values of H at the */
/*            abscissae unless IER < 0.  YS(N) = YS(1) if */
/*            PERIOD = TRUE. */

/*       YP = Array of length N containing first derivative */
/*            values of H at the abscissae unless IER < 0. */
/*            YP(N) = YP(1) if PERIOD = TRUE. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered and the */
/*                     constraint is active:  Q2(YS) is ap- */
/*                     proximately equal to SM. */
/*             IER = 1 if no errors were encountered but the */
/*                     constraint is not active:  YS and YP */
/*                     are the values and derivatives of the */
/*                     linear function (constant function if */
/*                     PERIOD = TRUE) which minimizes Q2, and */
/*                     Q1 = 0. */
/*             IER = -1 if N, W, SM, or SMTOL is outside its */
/*                      valid range.  YS and YP are unaltered */
/*                      in this case. */
/*             IER = -I if X(I) .LE. X(I-1) for some I in the */
/*                      range 2,...,N.  YS and YP are unal- */
/*                      tered in this case. */

/* Modules required by SMCRV:  B2TRI or B2TRIP, SNHCSH, */
/*                               YPCOEF */

/* Intrinsic functions called by SMCRV:  ABS, SQRT */

/* *********************************************************** */


    /* Parameter adjustments */
    --yp;
    --ys;
    wk_dim1 = *n;
    wk_offset = wk_dim1 + 1;
    wk -= wk_offset;
    --w;
    --sigma;
    --y;
    --x;

    /* Function Body */
    nn = *n;
    per = *period;

/* Test for errors, and compute the components of the system */
/*   (normal equations) for the weighted least squares linear */
/*   fit. */

    *ier = -1;
    if (nn < 2 || nn < 3 && per || *sm <= 0. || *smtol <= 0. || *smtol >= 1.) 
	    {
	return 0;
    }
    c11 = 0.;
    c12 = 0.;
    c22 = 0.;
    r1 = 0.;
    r2 = 0.;
    xi = x[1] - 1.;
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wi = w[i__];
	if (wi <= 0.) {
	    return 0;
	}
	if (x[i__] <= xi) {
	    *ier = -i__;
	    return 0;
	}
	xi = x[i__];
	yi = y[i__];
	c22 += wi;
	r2 += wi * yi;
	if (! per) {
	    wixi = wi * xi;
	    c11 += wixi * xi;
	    c12 += wixi;
	    r1 += wixi * yi;
	}
/* L1: */
    }

/* Solve the system for (HP,H0), where HP is the derivative */
/*   (constant) and H0 = H(0). */

    if (per) {
	h0 = r2 / c22;
	hp = 0.;
    } else {
	h0 = (c11 * r2 - c12 * r1) / (c11 * c22 - c12 * c12);
	hp = (r1 - c12 * h0) / c11;
    }

/* Store function values and derivatives, and accumulate */
/*   Q2 = (Y-YS)**T*W*(Y-YS). */

    q2 = 0.;
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ys[i__] = hp * x[i__] + h0;
	yp[i__] = hp;
/* Computing 2nd power */
	d__1 = y[i__] - ys[i__];
	q2 += w[i__] * (d__1 * d__1);
/* L2: */
    }

/* Compute bounds on Q2 defined by SMTOL, and test for the */
/*   constraint satisfied by the linear fit. */

    q2min = *sm * (1. - *smtol);
    q2max = *sm * (*smtol + 1.);
    if (q2 <= q2max) {

/*   The constraint is satisfied by a linear function. */

	*ier = 1;
	return 0;
    }

/* Compute the matrix components for the linear system. */

    *ier = 0;
    nm1 = nn - 1;
    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dx = x[i__ + 1] - x[i__];
	sig = (d__1 = sigma[i__], abs(d__1));
	ypcoef_(&sig, &dx, &d__, &sd);
	wk[i__ + wk_dim1] = d__;
	wk[i__ + (wk_dim1 << 1)] = sd;
/* L3: */
    }

/* Compute G0 = G(0), and print a heading. */

    s = 1. / sqrt(*sm);
    g0 = 1. / sqrt(q2) - s;

/* G(P) is strictly increasing and concave, and G(0) < 0. */

/* Initialize parameters for the secant method.  The method */
/*   uses three points:  (P0,G0), (P,G), and (PNEG,GNEG), */
/*   where P0 and PNEG are defined implicitly by DP = P - P0 */
/*   and DMAX = P - PNEG. */

    p = *sm * 10.;
    dp = p;
    dmax__ = 0.;
    iter = 0;

/* Top of loop:  compute G and print a message.  For each */
/*               secant iteration, the following values are */
/*               printed:  P, G(P), and DP, where DP is the */
/*               change in P computed by linear interpolation */
/*               between the current point (P,G) and a previ- */
/*               ous point. */


L4:
    if (! per) {
	b2tri_(&nn, &x[1], &y[1], &w[1], &p, &wk[wk_offset], &wk[(wk_dim1 << 
		1) + 1], &wk[wk_dim1 * 3 + 1], &wk[(wk_dim1 << 2) + 1], &wk[
		wk_dim1 * 5 + 1], &wk[wk_dim1 * 6 + 1], &ys[1], &yp[1], &ierr)
		;
    } else {
	b2trip_(&nn, &x[1], &y[1], &w[1], &p, &wk[wk_offset], &wk[(wk_dim1 << 
		1) + 1], &wk[wk_dim1 * 3 + 1], &wk[(wk_dim1 << 2) + 1], &wk[
		wk_dim1 * 5 + 1], &wk[wk_dim1 * 6 + 1], &wk[wk_dim1 * 7 + 1], 
		&wk[(wk_dim1 << 3) + 1], &wk[wk_dim1 * 9 + 1], &wk[wk_dim1 * 
		10 + 1], &ys[1], &yp[1], &ierr);
    }
    q2 = 0.;
    i__1 = nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = y[i__] - ys[i__];
	q2 += w[i__] * (d__1 * d__1);
/* L5: */
    }
    g = 1. / sqrt(q2) - s;
    ++iter;
    p0 = p - dp;

/*   Test for convergence. */

    if (g == g0 || q2min <= q2 && q2 <= q2max) {
	return 0;
    }
    if (dmax__ != 0. || g > 0.) {
	goto L6;
    }

/*   Increase P until G(P) > 0. */

    p *= 10.;
    dp = p;
    goto L4;

/*   A bracketing interval [P0,P] has been found. */

L6:
    if (g0 * g <= 0.) {

/*   G0*G < 0.  Update (PNEG,GNEG) to (P0,G0) so that G */
/*     and GNEG always have opposite signs. */

	dmax__ = dp;
	gneg = g0;
    }

/*   Compute the change in P by linear interpolation between */
/*     (P0,G0) and (P,G). */

L7:
    dp = -g * dp / (g - g0);
    if (abs(dp) > abs(dmax__)) {

/*   G0*G > 0, and the new estimate would be outside of the */
/*     bracketing interval of length abs(DMAX).  Reset */
/*     (P0,G0) to (PNEG,GNEG). */

	dp = dmax__;
	g0 = gneg;
	goto L7;
    }

/*   Bottom of loop:  update P, DMAX, and G0. */

    p += dp;
    dmax__ += dp;
    g0 = g;
    goto L4;
} /* smcrv_ */

/* Subroutine */ int snhcsh_(doublereal *x, doublereal *sinhm, doublereal *
	coshm, doublereal *coshmm)
{
    /* Initialized data */

    static doublereal p1 = -351754.9648081513948;
    static doublereal p2 = -11561.4435765005216044;
    static doublereal p3 = -163.725857525983828727;
    static doublereal p4 = -.789474443963537015605;
    static doublereal q1 = -2110529.78884890840399;
    static doublereal q2 = 36157.8279834431989373;
    static doublereal q3 = -277.711081420602794433;
    static doublereal q4 = 1.;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal expx, f, p, q, ax, xc, xs, xsd2, xsd4;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   11/20/96 */

/*   This subroutine computes approximations to the modified */
/* hyperbolic functions defined below with relative error */
/* bounded by 3.4E-20 for a floating point number system with */
/* sufficient precision. */

/*   Note that the 21-digit constants in the data statements */
/* below may not be acceptable to all compilers. */

/* On input: */

/*       X = Point at which the functions are to be */
/*           evaluated. */

/* X is not altered by this routine. */

/* On output: */

/*       SINHM = sinh(X) - X. */

/*       COSHM = cosh(X) - 1. */

/*       COSHMM = cosh(X) - 1 - X*X/2. */

/* Modules required by SNHCSH:  None */

/* Intrinsic functions called by SNHCSH:  ABS, EXP */

/* *********************************************************** */


    ax = abs(*x);
    xs = ax * ax;
    if (ax <= .5) {

/* Approximations for small X: */

	xc = *x * xs;
	p = ((p4 * xs + p3) * xs + p2) * xs + p1;
	q = ((q4 * xs + q3) * xs + q2) * xs + q1;
	*sinhm = xc * (p / q);
	xsd4 = xs * .25;
	xsd2 = xsd4 + xsd4;
	p = ((p4 * xsd4 + p3) * xsd4 + p2) * xsd4 + p1;
	q = ((q4 * xsd4 + q3) * xsd4 + q2) * xsd4 + q1;
	f = xsd4 * (p / q);
	*coshmm = xsd2 * f * (f + 2.);
	*coshm = *coshmm + xsd2;
    } else {

/* Approximations for large X: */

	expx = exp(ax);
	*sinhm = -(1. / expx + ax + ax - expx) / 2.;
	if (*x < 0.) {
	    *sinhm = -(*sinhm);
	}
	*coshm = (1. / expx - 2. + expx) / 2.;
	*coshmm = *coshm - xs / 2.;
    }
    return 0;
} /* snhcsh_ */

doublereal store_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/10/92 */

/*   This function forces its argument X to be stored in a */
/* memory location, thus providing a means of determining */
/* floating point number characteristics (such as the machine */
/* precision) when it is necessary to avoid computation in */
/* high precision registers. */

/* On input: */

/*       X = Value to be stored. */

/* X is not altered by this function. */

/* On output: */

/*       STORE = Value of X after it has been stored and */
/*               possibly truncated or rounded to the single */
/*               precision word length. */

/* Modules required by STORE:  None */

/* *********************************************************** */

    stcom_1.y = *x;
    ret_val = stcom_1.y;
    return ret_val;
} /* store_ */

doublereal tsintl_(doublereal *a, doublereal *b, integer *n, doublereal *x, 
	doublereal *y, doublereal *yp, doublereal *sigma, integer *ier)
{
    /* Initialized data */

    static doublereal sbig = 85.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer imin, imax;
    static doublereal e;
    static integer i__;
    static doublereal s, t, u, b1, b2, d1, d2, e1, e2, s1, s2, y1, y2, cm;
    static integer il;
    static doublereal dx;
    static integer iu;
    static doublereal sm, tm, tp, xl, xu;
    extern /* Subroutine */ int snhcsh_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal cm1, cm2, sb1, sb2;
    static integer ip1;
    extern integer intrvl_(doublereal *, integer *, doublereal *);
    static doublereal sm1, sm2, cmm, sig, ems, sum, cmm1, cmm2;
    static integer ilp1, iup1;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   11/17/96 */

/*   This function computes the integral from A to B of a */
/* Hermite interpolatory tension spline H. */

/* On input: */

/*       A,B = Lower and upper limits of integration, re- */
/*             spectively.  Note that -TSINTL(B,A,...) = */
/*             TSINTL(A,B,...). */

/*       N = Number of data points.  N .GE. 2. */

/*       X = Array of length N containing the abscissae. */
/*           These must be in strictly increasing order: */
/*           X(I) < X(I+1) for I = 1,...,N-1. */

/*       Y = Array of length N containing data values. */
/*           H(X(I)) = Y(I) for I = 1,...,N. */

/*       YP = Array of length N containing first deriva- */
/*            tives.  HP(X(I)) = YP(I) for I = 1,...,N, where */
/*            HP denotes the derivative of H. */

/*       SIGMA = Array of length N-1 containing tension fac- */
/*               tors whose absolute values determine the */
/*               balance between cubic and linear in each */
/*               interval.  SIGMA(I) is associated with int- */
/*               erval (I,I+1) for I = 1,...,N-1. */

/* Input parameters are not altered by this function. */

/* On output: */

/*       IER = Error indicator: */
/*             IER = 0  if no errors were encountered and */
/*                      X(1) .LE. T .LE. X(N) for T = A and */
/*                      T = B, or A = B. */
/*             IER = 1  if no errors were encountered but */
/*                      extrapolation was necessary:  A or B */
/*                      not in the interval (X(1),X(N)). */
/*             IER = -1 IF N < 2. */
/*             IER = -2 if the abscissae are not in strictly */
/*                      increasing order.  Only those in or */
/*                      adjacent to the interval of integra- */
/*                      tion are tested. */

/*       TSINTL = Integral of H from A to B, or zero if */
/*                IER < 0. */

/* Modules required by TSINTL:  INTRVL, SNHCSH */

/* Intrinsic functions called by TSINTL:  ABS, EXP, MAX, MIN */

/* *********************************************************** */


    /* Parameter adjustments */
    --sigma;
    --yp;
    --y;
    --x;

    /* Function Body */
    if (*n < 2) {
	goto L7;
    }

/* Accumulate the integral from XL to XU in SUM. */

    xl = min(*a,*b);
    xu = max(*a,*b);
    sum = 0.;
    *ier = 0;
    if (xl == xu) {
	goto L6;
    }

/* Find left-end indexes of intervals containing XL and XU. */
/*   If XL < X(1) or XU > X(N), extrapolation is performed */
/*   using the leftmost or rightmost interval. */

    il = intrvl_(&xl, n, &x[1]);
    iu = intrvl_(&xu, n, &x[1]);
    if (xl < x[1] || xu > x[*n]) {
	*ier = 1;
    }
    ilp1 = il + 1;
    imin = il;
    if (xl == x[il]) {
	goto L2;
    }

/* Compute the integral from XL to X(IL+1). */

    dx = x[ilp1] - x[il];
    if (dx <= 0.) {
	goto L8;
    }
    u = x[ilp1] - xl;
    if (u == 0.) {
	goto L1;
    }
    b1 = u / dx;
    y2 = y[ilp1];
    s = (y2 - y[il]) / dx;
    s2 = yp[ilp1];
    d1 = s - yp[il];
    d2 = s2 - s;
    sig = (d__1 = sigma[il], abs(d__1));
    if (sig < 1e-9) {

/*   SIG = 0. */

	sum += u * (y2 - u * (s2 * 6. - b1 * (d2 * 4. + (b1 * 3. - 4.) * (d1 
		- d2))) / 12.);
    } else if (sig <= .5) {

/*   0 .LT. SIG .LE. .5. */

	sb1 = sig * b1;
	snhcsh_(&sig, &sm, &cm, &cmm);
	snhcsh_(&sb1, &sm1, &cm1, &cmm1);
	e = sig * sm - cmm - cmm;
/* Computing 2nd power */
	d__1 = sig / dx;
	sum = sum + u * (y2 - s2 * u / 2.) + ((cm * cmm1 - sm * sm1) * (d1 + 
		d2) + sig * (cm * sm1 - (sm + sig) * cmm1) * d2) / (d__1 * 
		d__1 * e);
    } else {

/*   SIG > .5. */

	sb1 = sig * b1;
	sb2 = sig - sb1;
	if (-sb1 > sbig || -sb2 > sbig) {
	    sum += u * (y2 - s * u / 2.);
	} else {
	    e1 = exp(-sb1);
	    e2 = exp(-sb2);
	    ems = e1 * e2;
	    tm = 1. - ems;
	    tp = ems + 1.;
	    t = sb1 * sb1 / 2. + 1.;
	    e = tm * (sig * tp - tm - tm);
/* Computing 2nd power */
	    d__1 = sig / dx;
	    sum = sum + u * (y2 - s2 * u / 2.) + (sig * tm * (tp * t - e1 - 
		    e2 - tm * sb1) * d2 - (tm * (tm * t - e1 + e2 - tp * sb1) 
		    + sig * (e1 * ems - e2 + sb1 * 2. * ems)) * (d1 + d2)) / (
		    d__1 * d__1 * e);
	}
    }

/* Test for XL and XU in the same interval. */

L1:
    imin = ilp1;
    if (il == iu) {
	goto L5;
    }

/* Add in the integral from X(IMIN) to X(J) for J = */
/*   Max(IL+1,IU). */

L2:
/* Computing MAX */
    i__1 = il, i__2 = iu - 1;
    imax = max(i__1,i__2);
    i__1 = imax;
    for (i__ = imin; i__ <= i__1; ++i__) {
	ip1 = i__ + 1;
	dx = x[ip1] - x[i__];
	if (dx <= 0.) {
	    goto L8;
	}
	sig = (d__1 = sigma[i__], abs(d__1));
	if (sig < 1e-9) {

/*   SIG = 0. */

	    sum += dx * ((y[i__] + y[ip1]) / 2. - dx * (yp[ip1] - yp[i__]) / 
		    12.);
	} else if (sig <= .5) {

/*   0 .LT. SIG .LE. .5. */

	    snhcsh_(&sig, &sm, &cm, &cmm);
	    e = sig * sm - cmm - cmm;
	    sum += dx * (y[i__] + y[ip1] - dx * e * (yp[ip1] - yp[i__]) / (
		    sig * sig * cm)) / 2.;
	} else {

/*   SIG > .5. */

	    ems = exp(-sig);
	    sum += dx * (y[i__] + y[ip1] - dx * (sig * (ems + 1.) / (1. - ems)
		     - 2.) * (yp[ip1] - yp[i__]) / (sig * sig)) / 2.;
	}
/* L3: */
    }

/* Add in the integral from X(IU) to XU if IU > IL. */

    if (il == iu) {
	goto L4;
    }
    iup1 = iu + 1;
    dx = x[iup1] - x[iu];
    if (dx <= 0.) {
	goto L8;
    }
    u = xu - x[iu];
    if (u == 0.) {
	goto L6;
    }
    b2 = u / dx;
    y1 = y[iu];
    s = (y[iup1] - y1) / dx;
    s1 = yp[iu];
    d1 = s - s1;
    d2 = yp[iup1] - s;
    sig = (d__1 = sigma[iu], abs(d__1));
    if (sig < 1e-9) {

/*   SIG = 0. */

	sum += u * (y1 + u * (s1 * 6. + b2 * (d1 * 4. + (4. - b2 * 3.) * (d1 
		- d2))) / 12.);
    } else if (sig <= .5) {

/*   0 .LT. SIG .LE. .5. */

	sb2 = sig * b2;
	snhcsh_(&sig, &sm, &cm, &cmm);
	snhcsh_(&sb2, &sm2, &cm2, &cmm2);
	e = sig * sm - cmm - cmm;
/* Computing 2nd power */
	d__1 = sig / dx;
	sum = sum + u * (y1 + s1 * u / 2.) + ((cm * cmm2 - sm * sm2) * (d1 + 
		d2) + sig * (cm * sm2 - (sm + sig) * cmm2) * d1) / (d__1 * 
		d__1 * e);
    } else {

/*   SIG > .5. */

	sb2 = sig * b2;
	sb1 = sig - sb2;
	if (-sb1 > sbig || -sb2 > sbig) {
	    sum += u * (y1 + s * u / 2.);
	} else {
	    e1 = exp(-sb1);
	    e2 = exp(-sb2);
	    ems = e1 * e2;
	    tm = 1. - ems;
	    tp = ems + 1.;
	    t = sb2 * sb2 / 2. + 1.;
	    e = tm * (sig * tp - tm - tm);
/* Computing 2nd power */
	    d__1 = sig / dx;
	    sum = sum + u * (y1 + s1 * u / 2.) + (sig * tm * (tp * t - e1 - 
		    e2 - tm * sb2) * d1 - (tm * (tm * t - e2 + e1 - tp * sb2) 
		    + sig * (e2 * ems - e1 + sb2 * 2. * ems)) * (d1 + d2)) / (
		    d__1 * d__1 * e);
	}
    }
    goto L6;

/* IL = IU and SUM contains the integral from XL to X(IL+1). */
/*   Subtract off the integral from XU to X(IL+1).  DX and */
/*   SIG were computed above. */

L4:
    y2 = y[ilp1];
    s = (y2 - y[il]) / dx;
    s2 = yp[ilp1];
    d1 = s - yp[il];
    d2 = s2 - s;

L5:
    u = x[ilp1] - xu;
    if (u == 0.) {
	goto L6;
    }
    b1 = u / dx;
    if (sig < 1e-9) {

/*   SIG = 0. */

	sum -= u * (y2 - u * (s2 * 6. - b1 * (d2 * 4. + (b1 * 3. - 4.) * (d1 
		- d2))) / 12.);
    } else if (sig <= .5) {

/*   0 .LT. SIG .LE. .5. */

	sb1 = sig * b1;
	snhcsh_(&sig, &sm, &cm, &cmm);
	snhcsh_(&sb1, &sm1, &cm1, &cmm1);
	e = sig * sm - cmm - cmm;
/* Computing 2nd power */
	d__1 = sig / dx;
	sum = sum - u * (y2 - s2 * u / 2.) - ((cm * cmm1 - sm * sm1) * (d1 + 
		d2) + sig * (cm * sm1 - (sm + sig) * cmm1) * d2) / (d__1 * 
		d__1 * e);
    } else {

/*   SIG > .5. */

	sb1 = sig * b1;
	sb2 = sig - sb1;
	if (-sb1 > sbig || -sb2 > sbig) {
	    sum -= u * (y2 - s * u / 2.);
	} else {
	    e1 = exp(-sb1);
	    e2 = exp(-sb2);
	    ems = e1 * e2;
	    tm = 1. - ems;
	    tp = ems + 1.;
	    t = sb1 * sb1 / 2. + 1.;
	    e = tm * (sig * tp - tm - tm);
/* Computing 2nd power */
	    d__1 = sig / dx;
	    sum = sum - u * (y2 - s2 * u / 2.) - (sig * tm * (tp * t - e1 - 
		    e2 - tm * sb1) * d2 - (tm * (tm * t - e1 + e2 - tp * sb1) 
		    + sig * (e1 * ems - e2 + sb1 * 2. * ems)) * (d1 + d2)) / (
		    d__1 * d__1 * e);
	}
    }

/* No errors were encountered.  Adjust the sign of SUM. */

L6:
    if (xl == *b) {
	sum = -sum;
    }
    ret_val = sum;
    return ret_val;

/* N < 2. */

L7:
    *ier = -1;
    ret_val = 0.;
    return ret_val;

/* Abscissae not strictly increasing. */

L8:
    *ier = -2;
    ret_val = 0.;
    return ret_val;
} /* tsintl_ */

/* Subroutine */ int tspbi_(integer *n, doublereal *x, doublereal *y, integer 
	*ncd, integer *iendc, logical *per, doublereal *b, doublereal *bmax, 
	integer *lwk, doublereal *wk, doublereal *yp, doublereal *sigma, 
	integer *icflg, integer *ier)
{
    /* Initialized data */

    static doublereal stol = 0.;
    static integer maxit = 49;
    static doublereal dyptol = .01;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer icnt, ierr, iter;
    static logical loop2;
    extern /* Subroutine */ int ypc1p_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *), ypc2p_(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;
    static doublereal e;
    static integer i__;
    extern /* Subroutine */ int sigbi_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal dsmax;
    static integer nn, nm1;
    static doublereal yp1, dyp, ypn;
    extern /* Subroutine */ int ypc1_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *), ypc2_(integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/08/92 */

/*   This subroutine computes a set of parameter values which */
/* define a Hermite interpolatory tension spline H(x).  The */
/* parameters consist of knot derivative values YP computed */
/* by Subroutine YPC1, YPC1P, YPC2, or YPC2P, and tension */
/* factors SIGMA chosen to satisfy user-specified constraints */
/* (by Subroutine SIGBI).  Refer to Subroutine TSPSI for an */
/* alternative method of computing tension factors. */

/*   Refer to Subroutine TSPSS for a means of computing */
/* parameters which define a smoothing curve rather than an */
/* interpolatory curve. */

/*   The tension spline may be evaluated by Subroutine TSVAL1 */
/* or Functions HVAL (values), HPVAL (first derivatives), */
/* HPPVAL (second derivatives), and TSINTL (integrals). */

/* On input: */

/*       N = Number of data points.  N .GE. 2 and N .GE. 3 if */
/*           PER = TRUE. */

/*       X = Array of length N containing a strictly in- */
/*           creasing sequence of abscissae:  X(I) < X(I+1) */
/*           for I = 1,...,N-1. */

/*       Y = Array of length N containing data values asso- */
/*           ciated with the abscissae.  H(X(I)) = Y(I) for */
/*           I = 1,...,N.  If NCD = 1 and PER = TRUE, Y(N) */
/*           is set to Y(1). */

/*       NCD = Number of continuous derivatives at the knots. */
/*             NCD = 1 or NCD = 2.  If NCD = 1, the YP values */
/*             are computed by local monotonicity-constrained */
/*             quadratic fits.  Otherwise, a linear system is */
/*             solved for the derivative values which result */
/*             in second derivative continuity.  This re- */
/*             quires iterating on calls to YPC2 or YPC2P and */
/*             calls to SIGBI, and generally results in more */
/*             nonzero tension factors (hence more expensive */
/*             evaluation). */

/*       IENDC = End condition indicator for NCD = 2 and PER */
/*               = FALSE (or dummy parameter otherwise): */
/*               IENDC = 0 if YP(1) and YP(N) are to be com- */
/*                         puted by monotonicity-constrained */
/*                         parabolic fits to the first three */
/*                         and last three points, respective- */
/*                         ly.  This is identical to the */
/*                         values computed by YPC1. */
/*               IENDC = 1 if the first derivatives of H at */
/*                         X(1) and X(N) are user-specified */
/*                         in YP(1) and YP(N), respectively. */
/*               IENDC = 2 if the second derivatives of H at */
/*                         X(1) and X(N) are user-specified */
/*                         in YP(1) and YP(N), respectively. */
/*               IENDC = 3 if the end conditions are to be */
/*                         computed by Subroutine ENDSLP and */
/*                         vary with SIGMA(1) and SIGMA(N-1). */

/*       PER = Logical variable with value TRUE if and only */
/*             H(x) is to be a periodic function with period */
/*             X(N)-X(1).  It is assumed without a test that */
/*             Y(N) = Y(1) in this case.  On output, YP(N) = */
/*             YP(1).  If H(x) is one of the components of a */
/*             parametric curve, this option may be used to */
/*             obtained a closed curve. */

/*       B = Array dimensioned 5 by N-1 containing bounds or */
/*           flags which define the constraints.  For I = 1 */
/*           to N-1, column I defines the constraints associ- */
/*           ated with interval (X(I),X(I+1)) as follows: */

/*             B(1,I) is an upper bound on H */
/*             B(2,I) is a lower bound on H */
/*             B(3,I) is an upper bound on HP */
/*             B(4,I) is a lower bound on HP */
/*             B(5,I) specifies the required sign of HPP */

/*           where HP and HPP denote the first and second */
/*           derivatives of H, respectively.  A null con- */
/*           straint is specified by abs(B(K,I)) .GE. BMAX */
/*           for K < 5, or B(5,I) = 0:   B(1,I) .GE. BMAX, */
/*           B(2,I) .LE. -BMAX, B(3,I) .GE. BMAX, B(4,I) .LE. */
/*           -BMAX, or B(5,I) = 0.  Any positive value of */
/*           B(5,I) specifies that H should be convex, a */
/*           negative values specifies that H should be con- */
/*           cave, and 0 specifies that no restriction be */
/*           placed on HPP.  Refer to Functions SIG0, SIG1, */
/*           and SIG2 for definitions of valid constraints. */

/*       BMAX = User-defined value of infinity which, when */
/*              used as an upper bound in B (or when when */
/*              its negative is used as a lower bound), */
/*              specifies that no constraint is to be en- */
/*              forced. */

/*       LWK = Length of work space WK: */
/*             LWK GE 2N-2 if NCD = 2 and PER = FALSE */
/*             LWK GE 3N-3 if NCD = 2 and PER = TRUE */

/*   The above parameters, except possibly Y(N), are not */
/* altered by this routine. */

/*       WK = Array of length at least LWK to be used as */
/*            temporary work space. */

/*       YP = Array of length .GE. N containing end condition */
/*            values in positions 1 and N if NCD = 2 and */
/*            IENDC = 1 or IENDC = 2. */

/*       SIGMA = Array of length .GE. N-1. */

/*       ICFLG = Array of length .GE. N-1. */

/* On output: */

/*       WK = Array containing convergence parameters in the */
/*            first two locations if IER > 0 (NCD = 2 and */
/*            no error other than invalid constraints was */
/*            encountered): */
/*            WK(1) = Maximum relative change in a component */
/*                    of YP on the last iteration. */
/*            WK(2) = Maximum relative change in a component */
/*                    of SIGMA on the last iteration. */

/*       YP = Array containing derivatives of H at the */
/*            abscissae.  YP is not altered if -3 < IER < 0, */
/*            and YP is only partially defined if IER = -4. */

/*       SIGMA = Array containing tension factors for which */
/*               H(x) satisfies the constraints defined by B. */
/*               SIGMA(I) is associated with interval (X(I), */
/*               X(I+1)) for I = 1,...,N-1.  If infinite ten- */
/*               sion is required in interval I, then */
/*               SIGMA(I) = 85 (and H is an approximation to */
/*               the linear interpolant on the interval), */
/*               and if no constraint is specified in the */
/*               interval, then SIGMA(I) = 0, and thus H is */
/*               cubic.  Invalid constraints are treated as */
/*               null constraints.  SIGMA is not altered if */
/*               -3 < IER < 0 (unless IENDC is invalid), and */
/*               SIGMA is the zero vector if IER = -4 or */
/*               IENDC (if used) is invalid. */

/*       ICFLG = Array of invalid constraint flags associated */
/*               with intervals.  For I = 1 to N-1, ICFLG(I) */
/*               is a 5-bit value b5b4b3b2b1, where bK = 1 if */
/*               and only if constraint K cannot be satis- */
/*               fied.  Thus, all constraints in interval I */
/*               are satisfied if and only if ICFLG(I) = 0 */
/*               (and IER .GE. 0).  ICFLG is not altered if */
/*               IER < 0. */

/*       IER = Error indicator or iteration count: */
/*             IER = IC .GE. 0 if no errors were encountered */
/*                      (other than invalid constraints) and */
/*                      IC calls to SIGBI and IC+1 calls to */
/*                      YPC1, YPC1P, YPC2 or YPC2P were */
/*                      employed.  (IC = 0 if NCD = 1). */
/*             IER = -1 if N, NCD, or IENDC is outside its */
/*                      valid range. */
/*             IER = -2 if LWK is too small. */
/*             IER = -4 if the abscissae X are not strictly */
/*                      increasing. */

/* Modules required by TSPBI:  ENDSLP, SIG0, SIG1, SIG2, */
/*                               SIGBI, SNHCSH, STORE, */
/*                               YPCOEF, YPC1, YPC1P, YPC2, */
/*                               YPC2P */

/* Intrinsic functions called by TSPBI:  ABS, MAX */

/* *********************************************************** */


    /* Parameter adjustments */
    --icflg;
    --sigma;
    --yp;
    b -= 6;
    --y;
    --x;
    --wk;

    /* Function Body */

/* Convergence parameters: */

/*   STOL = Absolute tolerance for SIGBI. */
/*   MAXIT = Maximum number of YPC2/SIGBI iterations for */
/*             each loop if NCD = 2. */
/*   DYPTOL = Bound on the maximum relative change in a */
/*              component of YP defining convergence of */
/*              the YPC2/SIGBI iteration when NCD = 2. */

    nn = *n;
    nm1 = nn - 1;

/* Test for invalid input parameters N, NCD, or LWK. */

    if (nn < 2 || *per && nn < 3 || *ncd < 1 || *ncd > 2) {
	goto L11;
    }
    if (*ncd == 2 && (*lwk < nm1 << 1 || *per && *lwk < nm1 * 3)) {
	goto L12;
    }

/* Initialize iteration count ITER, and initialize SIGMA to */
/*   zeros. */

    iter = 0;
    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sigma[i__] = 0.;
/* L1: */
    }
    if (*ncd == 1) {

/* NCD = 1. */

	if (! (*per)) {
	    ypc1_(&nn, &x[1], &y[1], &yp[1], &ierr);
	} else {
	    ypc1p_(&nn, &x[1], &y[1], &yp[1], &ierr);
	}
	if (ierr != 0) {
	    goto L14;
	}

/*   Compute tension factors. */

	sigbi_(&nn, &x[1], &y[1], &yp[1], &stol, &b[6], bmax, &sigma[1], &
		icflg[1], &dsmax, &ierr);
	goto L10;
    }

/* NCD = 2. */

    if (! (*per)) {

/*   Nonperiodic case:  call YPC2 and test for IENDC or X */
/*     invalid. */

	yp1 = yp[1];
	ypn = yp[nn];
	ypc2_(&nn, &x[1], &y[1], &sigma[1], iendc, iendc, &yp1, &ypn, &wk[1], 
		&yp[1], &ierr);
	if (ierr == 1) {
	    goto L11;
	}
	if (ierr > 1) {
	    goto L14;
	}
    } else {

/*   Periodic fit:  call YPC2P. */

	ypc2p_(&nn, &x[1], &y[1], &sigma[1], &wk[1], &yp[1], &ierr);
	if (ierr > 1) {
	    goto L14;
	}
    }
    loop2 = FALSE_;

/*   Iterate on calls to SIGBI and YPC2 (or YPC2P).  The */
/*     first N-1 WK locations are used to store the deriva- */
/*     tive estimates YP from the previous iteration. */

/*   LOOP2 is TRUE iff tension factors are not allowed to */
/*         decrease between iterations (loop 1 failed to */
/*         converge with MAXIT iterations). */
/*   DYP is the maximum relative change in a component of YP. */
/*   ICNT is the number of tension factors which were altered */
/*        by SIGBI. */
/*   DSMAX is the maximum relative change in a component of */
/*         SIGMA. */

L2:
    i__1 = maxit;
    for (iter = 1; iter <= i__1; ++iter) {
	dyp = 0.;
	i__2 = nm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    wk[i__] = yp[i__];
/* L3: */
	}
	sigbi_(&nn, &x[1], &y[1], &yp[1], &stol, &b[6], bmax, &sigma[1], &
		icflg[1], &dsmax, &icnt);
	if (! (*per)) {
	    ypc2_(&nn, &x[1], &y[1], &sigma[1], iendc, iendc, &yp1, &ypn, &wk[
		    nn], &yp[1], &ierr);
	} else {
	    ypc2p_(&nn, &x[1], &y[1], &sigma[1], &wk[nn], &yp[1], &ierr);
	}
	i__2 = nm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    e = (d__1 = yp[i__] - wk[i__], abs(d__1));
	    if (wk[i__] != 0.) {
		e /= (d__1 = wk[i__], abs(d__1));
	    }
	    dyp = max(dyp,e);
/* L4: */
	}
	if (icnt == 0 || dyp <= dyptol) {
	    goto L7;
	}
	if (! loop2) {

/*   Loop 1:  reinitialize SIGMA to zeros. */

	    i__2 = nm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		sigma[i__] = 0.;
/* L5: */
	    }
	}
/* L6: */
    }

/*   The loop failed to converge within MAXIT iterations. */

    iter = maxit;
    if (! loop2) {
	loop2 = TRUE_;
	goto L2;
    }

/* Store convergence parameters. */

L7:
    wk[1] = dyp;
    wk[2] = dsmax;
    if (loop2) {
	iter += maxit;
    }

/* No error encountered. */

L10:
    *ier = iter;
    return 0;

/* Invalid input parameter N, NCD, or IENDC. */

L11:
    *ier = -1;
    return 0;

/* LWK too small. */

L12:
    *ier = -2;
    return 0;

/* Abscissae are not strictly increasing. */

L14:
    *ier = -4;
    return 0;
} /* tspbi_ */

/* Subroutine */ int tspbp_(integer *n, doublereal *x, doublereal *y, integer 
	*ncd, integer *iendc, logical *per, doublereal *bl, doublereal *bu, 
	doublereal *bmax, integer *lwk, doublereal *wk, doublereal *t, 
	doublereal *xp, doublereal *yp, doublereal *sigma, integer *ier)
{
    /* Initialized data */

    static doublereal stol = 0.;
    static integer maxit = 49;
    static doublereal dyptol = .01;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer icnt, ierr, iter;
    static logical loop2;
    extern /* Subroutine */ int ypc1p_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *), ypc2p_(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;
    static integer i__;
    extern /* Subroutine */ int sigbp_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;
    static doublereal dsmax;
    extern /* Subroutine */ int arcl2d_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *);
    static integer nn;
    static doublereal ex, ey;
    static integer nm1;
    static doublereal xp1, yp1;
    static integer n2m1;
    static doublereal dyp, xpn, ypn;
    extern /* Subroutine */ int ypc1_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *), ypc2_(integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/08/92 */

/*   This subroutine computes a set of values which define a */
/* parametric planar curve C(t) = (H1(t),H2(t)) whose compo- */
/* nents are Hermite interpolatory tension splines.  The */
/* output values consist of parameters (knots) T computed by */
/* ARCL2D, knot derivative values XP and YP computed by Sub- */
/* routine YPC1, YPC1P, YPC2, or YPC2P, and tension factors */
/* SIGMA chosen (by Subroutine SIGBP) to satisfy user- */
/* specified bounds on the distance between C(t) and the */
/* polygonal curve associated with the data points (refer to */
/* BL and BU below). */

/*   Refer to Subroutine TSPSP for an alternative method of */
/* computing tension factors. */

/*   The tension splines may be evaluated by Subroutine */
/* TSVAL2 or Functions HVAL (values), HPVAL (first deriva- */
/* tives), HPPVAL (second derivatives), and TSINTL */
/* (integrals). */

/* On input: */

/*       N = Number of knots and data points.  N .GE. 2 and */
/*           N .GE. 3 if PER = TRUE. */

/*       X,Y = Arrays of length N containing the Cartesian */
/*             coordinates of an ordered sequence of data */
/*             points C(I), I = 1 to N, such that C(I) .NE. */
/*             C(I+1).  C(t) is constrained to pass through */
/*             these points.  In the case of a closed curve */
/*             (PER = TRUE), the first and last points should */
/*             coincide.  (X(N) and Y(N) are set to X(1) and */
/*             Y(1) if NCD = 1, but not altered if NCD = 2, */
/*             in this case.) */

/*       NCD = Number of continuous derivatives at the knots. */
/*             NCD = 1 or NCD = 2.  If NCD = 1, XP and YP are */
/*             computed by local monotonicity-constrained */
/*             quadratic fits.  Otherwise, a linear system is */
/*             solved for the derivative values which result */
/*             in second derivative continuity.  This re- */
/*             quires iterating on calls to YPC2 or YPC2P and */
/*             calls to SIGBP, and generally results in more */
/*             nonzero tension factors (hence more expensive */
/*             evaluation). */

/*       IENDC = End condition indicator for NCD = 2 and PER */
/*               = FALSE (or dummy parameter otherwise): */
/*               IENDC = 0 if XP(1), XP(N), YP(1), and YP(N) */
/*                         are to be computed by monotonicity- */
/*                         constrained parabolic fits (YPC1). */
/*               IENDC = 1 if the first derivatives of H1 at */
/*                         the left and right endpoints are */
/*                         user-specified in XP(1) and XP(N), */
/*                         respectively, and the first deriv- */
/*                         atives of H2 at the ends are */
/*                         specified in YP(1) and YP(N). */
/*               IENDC = 2 if the second derivatives of H1 */
/*                         and H2 at the endpoints are user- */
/*                         specified in XP(1), XP(N), YP(1), */
/*                         and YP(N). */
/*               IENDC = 3 if the end conditions are to be */
/*                         computed by Subroutine ENDSLP and */
/*                         vary with SIGMA(1) and SIGMA(N-1). */

/*       PER = Logical variable with value TRUE if and only */
/*             a closed curve is to be constructed -- H1(t) */
/*             and H2(t) are to be periodic functions with */
/*             period T(N)-T(1), where T(1) and T(N) are the */
/*             parameter values associated with the first and */
/*             last data points.  It is assumed that X(N) = */
/*             X(1) and Y(N) = Y(1) in this case, and, on */
/*             output, XP(N) = XP(1) and YP(N) = YP(1). */

/*       BL,BU = Arrays of length N-1 containing (for each */
/*               knot subinterval) lower and upper bounds, */
/*               respectively, on the signed perpendicular */
/*               distance d(t) = (C2-C1)/DC X (C(t)-C1), */
/*               where C1 and C2 are the ordered data points */
/*               associated with the interval, and DC is the */
/*               interval length (and length of the line seg- */
/*               ment C1-C2).  Note that d(t) > 0 iff C(t) */
/*               lies strictly to the left of the line seg- */
/*               ment as viewed from C1 toward C2.  For I = */
/*               1 to N-1, SIGMA(I) is chosen to be as small */
/*               as possible within the constraint that */
/*               BL(I) .LE. d(t) .LE. BU(I) for all t in the */
/*               interval.  BL(I) < 0 and BU(I) > 0 for I = 1 */
/*               to N-1.  A null constraint is specified by */
/*               BL(I) .LE. -BMAX or BU(I) .GE. BMAX. */

/*       BMAX = User-defined value of infinity which, when */
/*              used as an upper bound in BU (or when its */
/*              negative is used as a lower bound in BL), */
/*              specifies that no constraint is to be en- */
/*              forced. */

/*       LWK = Length of work space WK: */
/*             LWK GE 3N-3 if NCD = 2 and PER = FALSE */
/*             LWK GE 4N-4 if NCD = 2 and PER = TRUE */

/*   The above parameters, except possibly X(N) and Y(N), are */
/* not altered by this routine. */

/*       WK = Array of length .GE. LWK to be used as tempor- */
/*            ary work space. */

/*       T = Array of length .GE. N. */

/*       XP = Array of length .GE. N containing end condition */
/*            values in positions 1 and N if NCD = 2 and */
/*            IENDC = 1 or IENDC = 2. */

/*       YP = Array of length .GE. N containing end condition */
/*            values in positions 1 and N if NCD = 2 and */
/*            IENDC = 1 or IENDC = 2. */

/*       SIGMA = Array of length .GE. N-1. */

/* On output: */

/*       WK = Array containing convergence parameters in the */
/*            first two locations if IER > 0 (NCD = 2 and */
/*            no error was encountered): */
/*            WK(1) = Maximum relative change in a component */
/*                    of XP or YP on the last iteration. */
/*            WK(2) = Maximum relative change in a component */
/*                    of SIGMA on the last iteration. */

/*       T = Array containing parameter values computed by */
/*           Subroutine ARCL2D unless IER = -1 or IER = -2. */
/*           T is only partially defined if IER = -4. */

/*       XP = Array containing derivatives of H1 at the */
/*            knots.  XP is not altered if -5 < IER < 0, */
/*            and XP is only partially defined if IER = -6. */

/*       YP = Array containing derivatives of H2 at the */
/*            knots.  YP is not altered if -5 < IER < 0, */
/*            and YP is only partially defined if IER = -6. */

/*       SIGMA = Array containing tension factors for which */
/*               C(t) satisfies the constraints defined by */
/*               BL and BU.  SIGMA(I) is associated with */
/*               interval (T(I),T(I+1)) for I = 1,...,N-1. */
/*               SIGMA(I) is limited to 85 (in which case */
/*               C(t) is indistinguishable from the line */
/*               segment associated with the interval), and */
/*               if no constraint is specified in the */
/*               interval, then SIGMA(I) = 0, and thus H1 and */
/*               H2 are cubic functions of t.  SIGMA is not */
/*               altered if -5 < IER < 0 (unless IENDC in */
/*               invalid), and SIGMA is the zero vector if */
/*               IER = -6 or IENDC (if used) is invalid. */

/*       IER = Error indicator or iteration count: */
/*             IER = IC .GE. 0 if no errors were encountered */
/*                      and IC calls to SIGBP and IC+1 calls */
/*                      to YPC1, YPC1P, YPC2 or YPC2P were */
/*                      employed.  (IC = 0 if NCD = 1). */
/*             IER = -1 if N, NCD, or IENDC is outside its */
/*                      valid range. */
/*             IER = -2 if LWK is too small. */
/*             IER = -4 if a pair of adjacent data points */
/*                      coincide:  X(I) = X(I+1) and Y(I) = */
/*                      Y(I+1) for some I in the range 1 to */
/*                      N-1. */
/*             IER = -5 if BL(I) .GE. 0 or BU(I) .LE. 0 for */
/*                      some I in the range 1 to N-1. */
/*                      SIGMA(J) = 0 for J .GE. I in this */
/*                      case. */
/*             IER = -6 if invalid knots T were returned by */
/*                      ARCL2D.  This should not occur. */

/* Modules required by TSPBP:  ARCL2D, ENDSLP, SIGBP, SNHCSH, */
/*                               STORE, YPCOEF, YPC1, YPC1P, */
/*                               YPC2, YPC2P */

/* Intrinsic functions called by TSPBP:  ABS, MAX */

/* *********************************************************** */


    /* Parameter adjustments */
    --sigma;
    --yp;
    --xp;
    --t;
    --bu;
    --bl;
    --y;
    --x;
    --wk;

    /* Function Body */

/* Convergence parameters: */

/*   STOL = Absolute tolerance for SIGBP. */
/*   MAXIT = Maximum number of YPC2/SIGBP iterations for each */
/*             loop in NCD = 2. */
/*   DYPTOL = Bound on the maximum relative change in a */
/*              component of XP or YP defining convergence */
/*              of the YPC2/SIGBP iteration when NCD = 2. */

    nn = *n;
    nm1 = nn - 1;

/* Test for invalid input parameters N, NCD, or LWK. */

    n2m1 = (nn << 1) - 1;
    if (nn < 2 || *per && nn < 3 || *ncd < 1 || *ncd > 2) {
	goto L11;
    }
    if (*ncd == 2 && (*lwk < nm1 * 3 || *per && *lwk < nm1 << 2)) {
	goto L12;
    }

/* Compute the sequence of parameter values T. */

    arcl2d_(&nn, &x[1], &y[1], &t[1], &ierr);
    if (ierr > 0) {
	goto L14;
    }

/* Initialize iteration count ITER, and initialize SIGMA to */
/*   zeros. */

    iter = 0;
    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sigma[i__] = 0.;
/* L1: */
    }
    if (*ncd == 1) {

/* NCD = 1. */

	if (! (*per)) {
	    ypc1_(&nn, &t[1], &x[1], &xp[1], &ierr);
	    ypc1_(&nn, &t[1], &y[1], &yp[1], &ierr);
	} else {
	    ypc1p_(&nn, &t[1], &x[1], &xp[1], &ierr);
	    ypc1p_(&nn, &t[1], &y[1], &yp[1], &ierr);
	}
	if (ierr != 0) {
	    goto L16;
	}

/*   Compute tension factors. */

	sigbp_(&nn, &x[1], &y[1], &xp[1], &yp[1], &stol, &bl[1], &bu[1], bmax,
		 &sigma[1], &dsmax, &ierr);
	if (ierr < 0) {
	    goto L15;
	}
	goto L10;
    }

/* NCD = 2. */

    if (! (*per)) {

/*   Nonperiodic case:  call YPC2 and test for IENDC invalid. */

	xp1 = xp[1];
	xpn = xp[nn];
	ypc2_(&nn, &t[1], &x[1], &sigma[1], iendc, iendc, &xp1, &xpn, &wk[1], 
		&xp[1], &ierr);
	yp1 = yp[1];
	ypn = yp[nn];
	ypc2_(&nn, &t[1], &y[1], &sigma[1], iendc, iendc, &yp1, &ypn, &wk[1], 
		&yp[1], &ierr);
	if (ierr == 1) {
	    goto L11;
	}
	if (ierr > 1) {
	    goto L16;
	}
    } else {

/*   Periodic fit:  call YPC2P. */

	ypc2p_(&nn, &t[1], &x[1], &sigma[1], &wk[1], &xp[1], &ierr);
	ypc2p_(&nn, &t[1], &y[1], &sigma[1], &wk[1], &yp[1], &ierr);
	if (ierr != 0) {
	    goto L16;
	}
    }
    loop2 = FALSE_;

/*   Iterate on calls to SIGBP and YPC2 (or YPC2P).  The */
/*     first 2N-2 WK locations are used to store the deriva- */
/*     tive estimates XP and YP from the previous iteration. */

/*   LOOP2 is TRUE iff tension factors are not allowed to */
/*         decrease between iterations (loop 1 failed to */
/*         converge with MAXIT iterations). */
/*   DYP is the maximum relative change in a component of XP */
/*       or YP. */
/*   ICNT is the number of tension factors which were altered */
/*        by SIGBP. */
/*   DSMAX is the maximum relative change in a component of */
/*         SIGMA. */

L2:
    i__1 = maxit;
    for (iter = 1; iter <= i__1; ++iter) {
	dyp = 0.;
	i__2 = nm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    wk[i__] = xp[i__];
	    wk[nm1 + i__] = yp[i__];
/* L3: */
	}
	sigbp_(&nn, &x[1], &y[1], &xp[1], &yp[1], &stol, &bl[1], &bu[1], bmax,
		 &sigma[1], &dsmax, &icnt);
	if (icnt < 0) {
	    goto L15;
	}
	if (! (*per)) {
	    ypc2_(&nn, &t[1], &x[1], &sigma[1], iendc, iendc, &xp1, &xpn, &wk[
		    n2m1], &xp[1], &ierr);
	    ypc2_(&nn, &t[1], &y[1], &sigma[1], iendc, iendc, &yp1, &ypn, &wk[
		    n2m1], &yp[1], &ierr);
	} else {
	    ypc2p_(&nn, &t[1], &x[1], &sigma[1], &wk[n2m1], &xp[1], &ierr);
	    ypc2p_(&nn, &t[1], &y[1], &sigma[1], &wk[n2m1], &yp[1], &ierr);
	}
	i__2 = nm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    ex = (d__1 = xp[i__] - wk[i__], abs(d__1));
	    if (wk[i__] != 0.) {
		ex /= (d__1 = wk[i__], abs(d__1));
	    }
	    ey = (d__1 = yp[i__] - wk[nm1 + i__], abs(d__1));
	    if (wk[nm1 + i__] != 0.) {
		ey /= (d__1 = wk[nm1 + i__], abs(d__1));
	    }
/* Computing MAX */
	    d__1 = max(dyp,ex);
	    dyp = max(d__1,ey);
/* L4: */
	}
	if (icnt == 0 || dyp <= dyptol) {
	    goto L7;
	}
	if (! loop2) {

/*   Loop 1:  reinitialize SIGMA to zeros. */

	    i__2 = nm1;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		sigma[i__] = 0.;
/* L5: */
	    }
	}
/* L6: */
    }

/*   The loop failed to converge within MAXIT iterations. */

    iter = maxit;
    if (! loop2) {
	loop2 = FALSE_;
	goto L2;
    }

/* Store convergence parameters. */

L7:
    wk[1] = dyp;
    wk[2] = dsmax;
    if (loop2) {
	iter += maxit;
    }

/* No error encountered. */

L10:
    *ier = iter;
    return 0;

/* Invalid input parameter N, NCD, or IENDC. */

L11:
    *ier = -1;
    return 0;

/* LWK too small. */

L12:
    *ier = -2;
    return 0;

/* Adjacent duplicate data points encountered. */

L14:
    *ier = -4;
    return 0;

/* Invalid constraint encountered by SIGBP. */

L15:
    *ier = -5;
    return 0;

/* Error flag returned by YPC1, YPC1P, YPC2, or YPC2P: */
/*   T is not strictly increasing. */

L16:
    *ier = -6;
    return 0;
} /* tspbp_ */

/* Subroutine */ int tspsi_(integer *n, doublereal *x, doublereal *y, integer 
	*ncd, integer *iendc, logical *per, logical *unifrm, integer *lwk, 
	doublereal *wk, doublereal *yp, doublereal *sigma, integer *ier)
{
    /* Initialized data */

    static doublereal sbig = 85.;
    static doublereal stol = 0.;
    static integer maxit = 99;
    static doublereal dyptol = .01;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer icnt, ierr, iter;
    extern /* Subroutine */ int sigs_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , ypc1p_(integer *, doublereal *, doublereal *, doublereal *, 
	    integer *), ypc2p_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static doublereal e;
    static integer i__;
    static doublereal dsmax;
    static integer nn, nm1;
    static doublereal yp1, sig, dyp, ypn;
    extern /* Subroutine */ int ypc1_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *), ypc2_(integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/08/92 */

/*   This subroutine computes a set of parameter values which */
/* define a Hermite interpolatory tension spline H(x).  The */
/* parameters consist of knot derivative values YP computed */
/* by Subroutine YPC1, YPC1P, YPC2, or YPC2P, and tension */
/* factors SIGMA computed by Subroutine SIGS (unless UNIFRM = */
/* TRUE).  Alternative methods for computing SIGMA are pro- */
/* vided by Subroutine TSPBI and Functions SIG0, SIG1, and */
/* SIG2. */

/*   Refer to Subroutine TSPSS for a means of computing */
/* parameters which define a smoothing curve rather than an */
/* interpolatory curve. */

/*   The tension spline may be evaluated by Subroutine TSVAL1 */
/* or Functions HVAL (values), HPVAL (first derivatives), */
/* HPPVAL (second derivatives), and TSINTL (integrals). */

/* On input: */

/*       N = Number of data points.  N .GE. 2 and N .GE. 3 if */
/*           PER = TRUE. */

/*       X = Array of length N containing a strictly in- */
/*           creasing sequence of abscissae:  X(I) < X(I+1) */
/*           for I = 1,...,N-1. */

/*       Y = Array of length N containing data values asso- */
/*           ciated with the abscissae.  H(X(I)) = Y(I) for */
/*           I = 1,...,N.  If NCD = 1 and PER = TRUE, Y(N) */
/*           is set to Y(1). */

/*       NCD = Number of continuous derivatives at the knots. */
/*             NCD = 1 or NCD = 2.  If NCD = 1, the YP values */
/*             are computed by local monotonicity-constrained */
/*             quadratic fits.  Otherwise, a linear system is */
/*             solved for the derivative values which result */
/*             in second derivative continuity.  Unless */
/*             UNIFRM = TRUE, this requires iterating on */
/*             calls to YPC2 or YPC2P and calls to SIGS, and */
/*             generally results in more nonzero tension */
/*             factors (hence more expensive evaluation). */

/*       IENDC = End condition indicator for NCD = 2 and PER */
/*               = FALSE (or dummy parameter otherwise): */
/*               IENDC = 0 if YP(1) and YP(N) are to be com- */
/*                         puted by monotonicity-constrained */
/*                         parabolic fits to the first three */
/*                         and last three points, respective- */
/*                         ly.  This is identical to the */
/*                         values computed by YPC1. */
/*               IENDC = 1 if the first derivatives of H at */
/*                         X(1) and X(N) are user-specified */
/*                         in YP(1) and YP(N), respectively. */
/*               IENDC = 2 if the second derivatives of H at */
/*                         X(1) and X(N) are user-specified */
/*                         in YP(1) and YP(N), respectively. */
/*               IENDC = 3 if the end conditions are to be */
/*                         computed by Subroutine ENDSLP and */
/*                         vary with SIGMA(1) and SIGMA(N-1). */

/*       PER = Logical variable with value TRUE if and only */
/*             H(x) is to be a periodic function with period */
/*             X(N)-X(1).  It is assumed without a test that */
/*             Y(N) = Y(1) in this case.  On output, YP(N) = */
/*             YP(1).  If H(x) is one of the components of a */
/*             parametric curve, this option may be used to */
/*             obtained a closed curve. */

/*       UNIFRM = Logical variable with value TRUE if and */
/*                only if constant (uniform) tension is to be */
/*                used.  The tension factor must be input in */
/*                SIGMA(1) in this case and must be in the */
/*                range 0 to 85.  If SIGMA(1) = 0, H(x) is */
/*                piecewise cubic (a cubic spline if NCD = */
/*                2), and as SIGMA increases, H(x) approaches */
/*                the piecewise linear interpolant.  If */
/*                UNIFRM = FALSE, tension factors are chosen */
/*                (by SIGS) to preserve local monotonicity */
/*                and convexity of the data.  This often */
/*                improves the appearance of the curve over */
/*                the piecewise cubic fit. */

/*       LWK = Length of work space WK:  no work space is */
/*             needed if NCD = 1; at least N-1 locations */
/*             are required if NCD = 2; another N-1 locations */
/*             are required if PER = TRUE; and an additional */
/*             N-1 locations are required for the convergence */
/*             test if SIGS is called (UNIFRM = FALSE): */

/*             LWK GE 0    if NCD=1 */
/*             LWK GE N-1  if NCD=2, PER=FALSE, UNIFRM=TRUE */
/*             LWK GE 2N-2 if NCD=2, PER=TRUE,  UNIFRM=TRUE */
/*             LWK GE 2N-2 if NCD=2, PER=FALSE, UNIFRM=FALSE */
/*             LWK GE 3N-3 if NCD=2, PER=TRUE,  UNIFRM=FALSE */

/*   The above parameters, except possibly Y(N), are not */
/* altered by this routine. */

/*       WK = Array of length at least LWK to be used as */
/*            temporary work space. */

/*       YP = Array of length .GE. N containing end condition */
/*            values in positions 1 and N if NCD = 2 and */
/*            IENDC = 1 or IENDC = 2. */

/*       SIGMA = Array of length .GE. N-1 containing a ten- */
/*               sion factor (0 to 85) in the first position */
/*               if UNIFRM = TRUE. */

/* On output: */

/*       WK = Array containing convergence parameters in the */
/*            first two locations if IER > 0 (NCD = 2 and */
/*            UNIFRM = FALSE): */
/*            WK(1) = Maximum relative change in a component */
/*                    of YP on the last iteration. */
/*            WK(2) = Maximum relative change in a component */
/*                    of SIGMA on the last iteration. */

/*       YP = Array containing derivatives of H at the */
/*            abscissae.  YP is not altered if -4 < IER < 0, */
/*            and YP is only partially defined if IER = -4. */

/*       SIGMA = Array containing tension factors.  SIGMA(I) */
/*               is associated with interval (X(I),X(I+1)) */
/*               for I = 1,...,N-1.  SIGMA is not altered if */
/*               -4 < IER < 0 (unless IENDC is invalid), and */
/*               SIGMA is constant (not optimal) if IER = -4 */
/*               or IENDC (if used) is invalid. */

/*       IER = Error indicator or iteration count: */
/*             IER = IC .GE. 0 if no errors were encountered */
/*                      and IC calls to SIGS and IC+1 calls */
/*                      to YPC1, YPC1P, YPC2 or YPC2P were */
/*                      employed.  (IC = 0 if NCD = 1). */
/*             IER = -1 if N, NCD, or IENDC is outside its */
/*                      valid range. */
/*             IER = -2 if LWK is too small. */
/*             IER = -3 if UNIFRM = TRUE and SIGMA(1) is out- */
/*                      side its valid range. */
/*             IER = -4 if the abscissae X are not strictly */
/*                      increasing. */

/* Modules required by TSPSI:  ENDSLP, SIGS, SNHCSH, STORE, */
/*                               YPCOEF, YPC1, YPC1P, YPC2, */
/*                               YPC2P */

/* Intrinsic functions called by TSPSI:  ABS, MAX */

/* *********************************************************** */


    /* Parameter adjustments */
    --sigma;
    --yp;
    --y;
    --x;
    --wk;

    /* Function Body */

/* Convergence parameters: */

/*   STOL = Absolute tolerance for SIGS */
/*   MAXIT = Maximum number of YPC2/SIGS iterations */
/*   DYPTOL = Bound on the maximum relative change in a */
/*              component of YP defining convergence of */
/*              the YPC2/SIGS iteration when NCD = 2 and */
/*              UNIFRM = FALSE */


/* Test for invalid input parameters (other than X and */
/*   IENDC). */

    nn = *n;
    nm1 = nn - 1;
    if (nn < 2 || *per && nn < 3 || *ncd < 1 || *ncd > 2) {
	goto L11;
    }
    if (*unifrm) {
	if (*ncd == 2 && (*lwk < nm1 || *per && *lwk < nm1 << 1)) {
	    goto L12;
	}
	sig = sigma[1];
	if (sig < 0. || sig > sbig) {
	    goto L13;
	}
    } else {
	if (*ncd == 2 && (*lwk < nm1 << 1 || *per && *lwk < nm1 * 3)) {
	    goto L12;
	}
	sig = 0.;
    }

/* Initialize iteration count ITER, and store uniform */
/*   tension factors, or initialize SIGMA to zeros. */

    iter = 0;
    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sigma[i__] = sig;
/* L1: */
    }
    if (*ncd == 1) {

/* NCD = 1. */

	if (! (*per)) {
	    ypc1_(&nn, &x[1], &y[1], &yp[1], &ierr);
	} else {
	    ypc1p_(&nn, &x[1], &y[1], &yp[1], &ierr);
	}
	if (ierr != 0) {
	    goto L14;
	}
	if (! (*unifrm)) {

/*   Call SIGS for UNIFRM = FALSE. */

	    sigs_(&nn, &x[1], &y[1], &yp[1], &stol, &sigma[1], &dsmax, &ierr);
	}
	goto L10;
    }

/* NCD = 2. */

    if (! (*per)) {

/*   Nonperiodic case:  call YPC2 and test for IENDC or X */
/*     invalid. */

	yp1 = yp[1];
	ypn = yp[nn];
	ypc2_(&nn, &x[1], &y[1], &sigma[1], iendc, iendc, &yp1, &ypn, &wk[1], 
		&yp[1], &ierr);
	if (ierr == 1) {
	    goto L11;
	}
	if (ierr > 1) {
	    goto L14;
	}
    } else {

/*   Periodic fit:  call YPC2P. */

	ypc2p_(&nn, &x[1], &y[1], &sigma[1], &wk[1], &yp[1], &ierr);
	if (ierr > 1) {
	    goto L14;
	}
    }
    if (*unifrm) {
	goto L10;
    }

/*   Iterate on calls to SIGS and YPC2 (or YPC2P).  The first */
/*     N-1 WK locations are used to store the derivative */
/*     estimates YP from the previous iteration. */

/*   DYP is the maximum relative change in a component of YP. */
/*   ICNT is the number of tension factors which were */
/*        increased by SIGS. */
/*   DSMAX is the maximum relative change in a component of */
/*         SIGMA. */

    i__1 = maxit;
    for (iter = 1; iter <= i__1; ++iter) {
	dyp = 0.;
	i__2 = nm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    wk[i__] = yp[i__];
/* L2: */
	}
	sigs_(&nn, &x[1], &y[1], &yp[1], &stol, &sigma[1], &dsmax, &icnt);
	if (! (*per)) {
	    ypc2_(&nn, &x[1], &y[1], &sigma[1], iendc, iendc, &yp1, &ypn, &wk[
		    nn], &yp[1], &ierr);
	} else {
	    ypc2p_(&nn, &x[1], &y[1], &sigma[1], &wk[nn], &yp[1], &ierr);
	}
	i__2 = nm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    e = (d__1 = yp[i__] - wk[i__], abs(d__1));
	    if (wk[i__] != 0.) {
		e /= (d__1 = wk[i__], abs(d__1));
	    }
	    dyp = max(dyp,e);
/* L3: */
	}
	if (icnt == 0 || dyp <= dyptol) {
	    goto L5;
	}
/* L4: */
    }
    iter = maxit;

/* Store convergence parameters in WK. */

L5:
    wk[1] = dyp;
    wk[2] = dsmax;

/* No error encountered. */

L10:
    *ier = iter;
    return 0;

/* Invalid input parameter N, NCD, or IENDC. */

L11:
    *ier = -1;
    return 0;

/* LWK too small. */

L12:
    *ier = -2;
    return 0;

/* UNIFRM = TRUE and SIGMA(1) outside its valid range. */

L13:
    *ier = -3;
    return 0;

/* Abscissae are not strictly increasing. */

L14:
    *ier = -4;
    return 0;
} /* tspsi_ */

/* Subroutine */ int tspsp_(integer *n, integer *nd, doublereal *x, 
	doublereal *y, doublereal *z__, integer *ncd, integer *iendc, logical 
	*per, logical *unifrm, integer *lwk, doublereal *wk, doublereal *t, 
	doublereal *xp, doublereal *yp, doublereal *zp, doublereal *sigma, 
	integer *ier)
{
    /* Initialized data */

    static doublereal sbig = 85.;
    static doublereal stol = 0.;
    static integer maxit = 99;
    static doublereal dyptol = .01;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer icnt, ierr, iter;
    extern /* Subroutine */ int sigs_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    , ypc1p_(integer *, doublereal *, doublereal *, doublereal *, 
	    integer *), ypc2p_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static integer i__;
    static doublereal dsmax;
    static logical scurv;
    extern /* Subroutine */ int arcl2d_(integer *, doublereal *, doublereal *,
	     doublereal *, integer *), arcl3d_(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static integer nn;
    static doublereal ex, ey, ez;
    static integer nm1, iw1;
    static doublereal xp1, yp1, zp1;
    static integer n2m2;
    static doublereal sig, dyp, xpn, ypn, zpn;
    extern /* Subroutine */ int ypc1_(integer *, doublereal *, doublereal *, 
	    doublereal *, integer *), ypc2_(integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   07/08/92 */

/*   This subroutine computes a set of values which define a */
/* parametric planar curve C(t) = (H1(t),H2(t)) or space */
/* curve C(t) = (H1(t),H2(t),H3(t)) whose components are Her- */
/* mite interpolatory tension splines.  The output values */
/* consist of parameters (knots) T computed by ARCL2D or */
/* ARCL3D, knot derivative values XP, YP, (and ZP) computed */
/* by Subroutine YPC1, YPC1P, YPC2, or YPC2P, and tension */
/* factors SIGMA computed by Subroutine SIGS (unless UNIFRM = */
/* TRUE). */

/*   Refer to Subroutine TSPSP for an alternative method of */
/* computing tension factors in the case of a planar curve. */

/*   The tension splines may be evaluated by Subroutine */
/* TSVAL2 (or TSVAL3) or Functions HVAL (values), HPVAL */
/* (first derivatives), HPPVAL (second derivatives), and */
/* TSINTL (integrals). */

/* On input: */

/*       N = Number of knots and data points.  N .GE. 2 and */
/*           N .GE. 3 if PER = TRUE. */

/*       ND = Number of dimensions: */
/*            ND = 2 if a planar curve is to be constructed. */
/*            ND = 3 if a space curve is to be constructed. */

/*       X,Y,Z = Arrays of length N containing the Cartesian */
/*               coordinates of an ordered sequence of data */
/*               points C(I), I = 1 to N, such that C(I) .NE. */
/*               C(I+1).  C(t) is constrained to pass through */
/*               these points.  Z is an unused dummy parame- */
/*               ter if ND = 2.  In the case of a closed curve */
/*               (PER = TRUE), the first and last points */
/*               should coincide.  In this case, X(N), Y(N), */
/*               (and Z(N)) are set to X(1), Y(1), (and Z(1)) */
/*               if NCD = 1, but not altered if NCD = 2. */

/*       NCD = Number of continuous derivatives at the knots. */
/*             NCD = 1 or NCD = 2.  If NCD = 1, XP, YP, (and */
/*             ZP) are computed by local monotonicity- */
/*             constrained quadratic fits.  Otherwise, a */
/*             linear system is solved for the derivative */
/*             values which result in second derivative con- */
/*             tinuity.  Unless UNIFRM = FALSE, this requires */
/*             iterating on calls to YPC2 or YPC2P and calls */
/*             to SIGS, and generally results in more nonzero */
/*             tension factors (hence more expensive evalua- */
/*             tion). */

/*       IENDC = End condition indicator for NCD = 2 and PER */
/*               = FALSE (or dummy parameter otherwise): */
/*               IENDC = 0 if XP(1), XP(N), YP(1), YP(N) (and */
/*                         ZP(1) and ZP(N)) are to be com- */
/*                         puted by monotonicity-constrained */
/*                         parabolic fits (YPC1). */
/*               IENDC = 1 if the first derivatives of H1 at */
/*                         the left and right endpoints are */
/*                         user-specified in XP(1) and XP(N), */
/*                         respectively, the first deriva- */
/*                         tives of H2 at the ends are */
/*                         specified in YP(1) and YP(N), and, */
/*                         if ND = 3, the first derivatives */
/*                         of H3 are specified in ZP(1) and */
/*                         ZP(N). */
/*               IENDC = 2 if the second derivatives of H1, */
/*                         H2, (and H3) at the endpoints are */
/*                         user-specified in XP(1), XP(N), */
/*                         YP(1), YP(N), (ZP(1), and ZP(N)). */
/*               IENDC = 3 if the end conditions are to be */
/*                         computed by Subroutine ENDSLP and */
/*                         vary with SIGMA(1) and SIGMA(N-1). */

/*       PER = Logical variable with value TRUE if and only */
/*             a closed curve is to be constructed -- H1(t), */
/*             H2(t), (and H3(t)) are to be periodic func- */
/*             tions with period T(N)-T(1), where T(1) and */
/*             T(N) are the parameter values associated with */
/*             the first and last data points.  It is assumed */
/*             in this case that X(N) = X(1), Y(N) = Y(1) */
/*             and, if ND = 3, Z(N) = Z(1), and, on output, */
/*             XP(N) = XP(1), YP(N) = YP(1), (and ZP(N) = */
/*             ZP(1) if ND = 3). */

/*       UNIFRM = Logical variable with value TRUE if and */
/*                only if constant (uniform) tension is to be */
/*                used.  The tension factor must be input in */
/*                SIGMA(1) in this case and must be in the */
/*                range 0 to 85.  If SIGMA(1) = 0, H(t) is */
/*                piecewise cubic (a cubic spline if NCD = */
/*                2), and as SIGMA increases, H(t) approaches */
/*                the piecewise linear interpolant, where H */
/*                is H1, H2, or H3.  If UNIFRM = FALSE, */
/*                tension factors are chosen (by SIGS) to */
/*                preserve local monotonicity and convexity */
/*                of the data.  This often improves the */
/*                appearance of the curve over the piecewise */
/*                cubic fitting functions. */

/*       LWK = Length of work space WK:  no work space is */
/*             needed if NCD = 1; at least N-1 locations */
/*             are required if NCD = 2; another N-1 locations */
/*             are required if PER = TRUE; and an additional */
/*             ND*(N-1) locations are required for the con- */
/*             vergence test if SIGS is called (UNIFRM = */
/*             FALSE): */
/*               If NCD=1 then LWK = 0 (not tested). */
/*               If NCD=2 then */

/*             LWK GE N-1          if PER=FALSE, UNIFRM=TRUE */
/*             LWK GE 2N-2         if PER=TRUE,  UNIFRM=TRUE */
/*             LWK GE (ND+1)*(N-1) if PER=FALSE, UNIFRM=FALSE */
/*             LWK GE (ND+2)*(N-1) if PER=TRUE,  UNIFRM=FALSE */

/*   The above parameters, except possibly X(N), Y(N), and */
/* Z(N), are not altered by this routine. */

/*       WK = Array of length .GE. LWK to be used as tempor- */
/*            ary work space. */

/*       T = Array of length .GE. N. */

/*       XP = Array of length .GE. N containing end condition */
/*            values in positions 1 and N if NCD = 2 and */
/*            IENDC = 1 or IENDC = 2. */

/*       YP = Array of length .GE. N containing end condition */
/*            values in positions 1 and N if NCD = 2 and */
/*            IENDC = 1 or IENDC = 2. */

/*       ZP = Dummy argument if ND=2, or, if ND=3, array of */
/*            length .GE. N containing end condition values */
/*            in positions 1 and N if NCD = 2 and IENDC = 1 */
/*            or IENDC = 2. */

/*       SIGMA = Array of length .GE. N-1 containing a ten- */
/*               sion factor (0 to 85) in the first position */
/*               if UNIFRM = TRUE. */

/* On output: */

/*       WK = Array containing convergence parameters in the */
/*            first two locations if IER > 0 (NCD = 2 and */
/*            UNIFRM = FALSE): */
/*            WK(1) = Maximum relative change in a component */
/*                    of XP, YP, or ZP on the last iteration. */
/*            WK(2) = Maximum relative change in a component */
/*                    of SIGMA on the last iteration. */

/*       T = Array containing parameter values computed by */
/*           Subroutine ARCL2D or ARCL3D unless -4 < IER < 0. */
/*           T is only partially defined if IER = -4. */

/*       XP = Array containing derivatives of H1 at the */
/*            knots.  XP is not altered if -5 < IER < 0, */
/*            and XP is only partially defined if IER = -6. */

/*       YP = Array containing derivatives of H2 at the */
/*            knots.  YP is not altered if -5 < IER < 0, */
/*            and YP is only partially defined if IER = -6. */

/*       ZP = Array containing derivatives of H3 at the knots */
/*            if ND=3.  ZP is not altered if -5 < IER < 0, */
/*            and ZP is only partially defined if IER = -6. */

/*       SIGMA = Array containing tension factors.  SIGMA(I) */
/*               is associated with interval (T(I),T(I+1)) */
/*               for I = 1,...,N-1.  SIGMA is not altered if */
/*               -5 < IER < 0 (unless IENDC is invalid), and */
/*               SIGMA is constant (not optimal) if IER = -6 */
/*               or IENDC (if used) is invalid. */

/*       IER = Error indicator or iteration count: */
/*             IER = IC .GE. 0 if no errors were encountered */
/*                      and IC calls to SIGS and IC+1 calls */
/*                      to YPC1, YPC1P, YPC2 or YPC2P were */
/*                      employed.  (IC = 0 if NCD = 1). */
/*             IER = -1 if N, NCD, or IENDC is outside its */
/*                      valid range. */
/*             IER = -2 if LWK is too small. */
/*             IER = -3 if UNIFRM = TRUE and SIGMA(1) is out- */
/*                      side its valid range. */
/*             IER = -4 if a pair of adjacent data points */
/*                      coincide:  X(I) = X(I+1), Y(I) = */
/*                      Y(I+1), (and Z(I) = Z(I+1)) for some */
/*                      I in the range 1 to N-1. */
/*             IER = -6 if invalid knots T were returned by */
/*                      ARCL2D or ARCL3D.  This should not */
/*                      occur. */

/* Modules required by TSPSP:  ARCL2D, ARCL3D, ENDSLP, SIGS, */
/*                               SNHCSH, STORE, YPCOEF, YPC1, */
/*                               YPC1P, YPC2, YPC2P */

/* Intrinsic functions called by TSPSP:  ABS, MAX */

/* *********************************************************** */


    /* Parameter adjustments */
    --sigma;
    --zp;
    --yp;
    --xp;
    --t;
    --z__;
    --y;
    --x;
    --wk;

    /* Function Body */

/* Convergence parameters: */

/*   STOL = Absolute tolerance for SIGS */
/*   MAXIT = Maximum number of YPC2/SIGS iterations */
/*   DYPTOL = Bound on the maximum relative change in a com- */
/*              ponent of XP, YP, or ZP defining convergence */
/*              of the YPC2/SIGS iteration when NCD = 2 and */
/*              UNIFRM = FALSE */


/* Test for invalid input parameters N, NCD, or LWK. */

    nn = *n;
    nm1 = nn - 1;
    n2m2 = nm1 << 1;
    if (nn < 2 || *per && nn < 3 || *ncd < 1 || *ncd > 2) {
	goto L11;
    }
    if (*unifrm) {
	if (*ncd == 2 && (*lwk < nm1 || *per && *lwk < n2m2)) {
	    goto L12;
	}
	sig = sigma[1];
	if (sig < 0. || sig > sbig) {
	    goto L13;
	}
    } else {
	if (*ncd == 2 && (*lwk < (*nd + 1) * nm1 || *per && *lwk < (*nd + 2) *
		 nm1)) {
	    goto L12;
	}
	sig = 0.;
    }

/* Compute the sequence of parameter values T. */

    scurv = *nd == 3;
    if (! scurv) {
	arcl2d_(&nn, &x[1], &y[1], &t[1], &ierr);
    } else {
	arcl3d_(&nn, &x[1], &y[1], &z__[1], &t[1], &ierr);
    }
    if (ierr > 0) {
	goto L14;
    }

/* Initialize iteration count ITER, and store uniform */
/*   tension factors, or initialize SIGMA to zeros. */

    iter = 0;
    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sigma[i__] = sig;
/* L1: */
    }
    if (*ncd == 1) {

/* NCD = 1. */

	if (! (*per)) {
	    ypc1_(&nn, &t[1], &x[1], &xp[1], &ierr);
	    ypc1_(&nn, &t[1], &y[1], &yp[1], &ierr);
	    if (scurv) {
		ypc1_(&nn, &t[1], &z__[1], &zp[1], &ierr);
	    }
	} else {
	    ypc1p_(&nn, &t[1], &x[1], &xp[1], &ierr);
	    ypc1p_(&nn, &t[1], &y[1], &yp[1], &ierr);
	    if (scurv) {
		ypc1p_(&nn, &t[1], &z__[1], &zp[1], &ierr);
	    }
	}
	if (ierr != 0) {
	    goto L16;
	}
	if (! (*unifrm)) {

/*   Call SIGS for UNIFRM = FALSE. */

	    sigs_(&nn, &t[1], &x[1], &xp[1], &stol, &sigma[1], &dsmax, &ierr);
	    sigs_(&nn, &t[1], &y[1], &yp[1], &stol, &sigma[1], &dsmax, &ierr);
	    if (scurv) {
		sigs_(&nn, &t[1], &z__[1], &zp[1], &stol, &sigma[1], &dsmax, &
			ierr);
	    }
	}
	goto L10;
    }

/* NCD = 2. */

    if (! (*per)) {

/*   Nonperiodic case:  call YPC2 and test for IENDC invalid. */

	xp1 = xp[1];
	xpn = xp[nn];
	ypc2_(&nn, &t[1], &x[1], &sigma[1], iendc, iendc, &xp1, &xpn, &wk[1], 
		&xp[1], &ierr);
	yp1 = yp[1];
	ypn = yp[nn];
	ypc2_(&nn, &t[1], &y[1], &sigma[1], iendc, iendc, &yp1, &ypn, &wk[1], 
		&yp[1], &ierr);
	if (scurv) {
	    zp1 = zp[1];
	    zpn = zp[nn];
	    ypc2_(&nn, &t[1], &z__[1], &sigma[1], iendc, iendc, &zp1, &zpn, &
		    wk[1], &zp[1], &ierr);
	}
	if (ierr == 1) {
	    goto L11;
	}
	if (ierr > 1) {
	    goto L16;
	}
    } else {

/*   Periodic fit:  call YPC2P. */

	ypc2p_(&nn, &t[1], &x[1], &sigma[1], &wk[1], &xp[1], &ierr);
	ypc2p_(&nn, &t[1], &y[1], &sigma[1], &wk[1], &yp[1], &ierr);
	if (scurv) {
	    ypc2p_(&nn, &t[1], &z__[1], &sigma[1], &wk[1], &zp[1], &ierr);
	}
	if (ierr != 0) {
	    goto L16;
	}
    }
    if (*unifrm) {
	goto L10;
    }

/*   Iterate on calls to SIGS and YPC2 (or YPC2P).  The */
/*     first ND*(N-1) WK locations are used to store the */
/*     derivative estimates XP, YP, (and ZP) from the */
/*     previous iteration.  IW1 is the first free WK location */
/*     following the stored derivatives. */

/*   DYP is the maximum relative change in a component of XP, */
/*       YP, or ZP. */
/*   ICNT is the number of tension factors which were */
/*        increased by SIGS. */
/*   DSMAX is the maximum relative change in a component of */
/*         SIGMA. */

    iw1 = *nd * nm1 + 1;
    i__1 = maxit;
    for (iter = 1; iter <= i__1; ++iter) {
	dyp = 0.;
	i__2 = nm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    wk[i__] = xp[i__];
	    wk[nm1 + i__] = yp[i__];
/* L2: */
	}
	if (scurv) {
	    i__2 = nm1;
	    for (i__ = 2; i__ <= i__2; ++i__) {
		wk[n2m2 + i__] = zp[i__];
/* L3: */
	    }
	}
	sigs_(&nn, &t[1], &x[1], &xp[1], &stol, &sigma[1], &dsmax, &icnt);
	sigs_(&nn, &t[1], &y[1], &yp[1], &stol, &sigma[1], &dsmax, &icnt);
	if (scurv) {
	    sigs_(&nn, &t[1], &z__[1], &zp[1], &stol, &sigma[1], &dsmax, &
		    icnt);
	}
	if (! (*per)) {
	    ypc2_(&nn, &t[1], &x[1], &sigma[1], iendc, iendc, &xp1, &xpn, &wk[
		    iw1], &xp[1], &ierr);
	    ypc2_(&nn, &t[1], &y[1], &sigma[1], iendc, iendc, &yp1, &ypn, &wk[
		    iw1], &yp[1], &ierr);
	    if (scurv) {
		ypc2_(&nn, &t[1], &z__[1], &sigma[1], iendc, iendc, &zp1, &
			zpn, &wk[iw1], &zp[1], &ierr);
	    }
	} else {
	    ypc2p_(&nn, &t[1], &x[1], &sigma[1], &wk[iw1], &xp[1], &ierr);
	    ypc2p_(&nn, &t[1], &y[1], &sigma[1], &wk[iw1], &yp[1], &ierr);
	    if (scurv) {
		ypc2p_(&nn, &t[1], &z__[1], &sigma[1], &wk[iw1], &zp[1], &
			ierr);
	    }
	}
	ez = 0.;
	i__2 = nm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    ex = (d__1 = xp[i__] - wk[i__], abs(d__1));
	    if (wk[i__] != 0.) {
		ex /= (d__1 = wk[i__], abs(d__1));
	    }
	    ey = (d__1 = yp[i__] - wk[nm1 + i__], abs(d__1));
	    if (wk[nm1 + i__] != 0.) {
		ey /= (d__1 = wk[nm1 + i__], abs(d__1));
	    }
	    if (scurv) {
		ez = (d__1 = zp[i__] - wk[n2m2 + i__], abs(d__1));
		if (wk[n2m2 + i__] != 0.) {
		    ez /= (d__1 = wk[n2m2 + i__], abs(d__1));
		}
	    }
/* Computing MAX */
	    d__1 = max(dyp,ex), d__1 = max(d__1,ey);
	    dyp = max(d__1,ez);
/* L4: */
	}
	if (icnt == 0 || dyp <= dyptol) {
	    goto L6;
	}
/* L5: */
    }
    iter = maxit;

/* Store convergence parameters in WK. */

L6:
    wk[1] = dyp;
    wk[2] = dsmax;

/* No error encountered. */

L10:
    *ier = iter;
    return 0;

/* Invalid input parameter N, NCD, or IENDC. */

L11:
    *ier = -1;
    return 0;

/* LWK too small. */

L12:
    *ier = -2;
    return 0;

/* UNIFRM = TRUE and SIGMA(1) outside its valid range. */

L13:
    *ier = -3;
    return 0;

/* Adjacent duplicate data points encountered. */

L14:
    *ier = -4;
    return 0;

/* Error flag returned by YPC1, YPC1P, YPC2, or YPC2P: */
/*   T is not strictly increasing. */

L16:
    *ier = -6;
    return 0;
} /* tspsp_ */

/* Subroutine */ int tspss_(integer *n, doublereal *x, doublereal *y, logical 
	*per, logical *unifrm, doublereal *w, doublereal *sm, doublereal *
	smtol, integer *lwk, doublereal *wk, doublereal *sigma, doublereal *
	ys, doublereal *yp, integer *nit, integer *ier)
{
    /* Initialized data */

    static doublereal sbig = 85.;
    static doublereal stol = 0.;
    static integer maxit = 99;
    static doublereal dystol = .01;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer icnt, ierr, iter;
    extern /* Subroutine */ int sigs_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *)
	    ;
    static doublereal e;
    static integer i__;
    static doublereal dsmax;
    extern /* Subroutine */ int smcrv_(integer *, doublereal *, doublereal *, 
	    doublereal *, logical *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, integer *);
    static integer nn, nm1;
    static doublereal sig, dys;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/10/92 */

/*   This subroutine computes a set of parameter values which */
/* define a smoothing tension spline H(x).  The parameters */
/* consist of knot values YS and derivatives YP computed */
/* by Subroutine SMCRV, and tension factors SIGMA computed by */
/* Subroutine SIGS (unless UNIFRM = TRUE).  The Hermite */
/* interpolatory tension spline H(x) defined by the knot */
/* values and derivatives has two continuous derivatives and */
/* satisfies either natural or periodic end conditions. */

/*   The tension spline may be evaluated by Subroutine TSVAL1 */
/* or Functions HVAL (values), HPVAL (first derivatives), */
/* HPPVAL (second derivatives), and TSINTL (integrals). */

/* On input: */

/*       N = Number of data points.  N .GE. 2 and N .GE. 3 if */
/*           PER = TRUE. */

/*       X = Array of length N containing a strictly in- */
/*           creasing sequence of abscissae:  X(I) < X(I+1) */
/*           for I = 1,...,N-1. */

/*       Y = Array of length N containing data values asso- */
/*           ciated with the abscissae.  If PER = TRUE, it is */
/*           assumed that Y(N) = Y(1). */

/*       PER = Logical variable with value TRUE if and only */
/*             H(x) is to be a periodic function with period */
/*             X(N)-X(1).  It is assumed without a test that */
/*             Y(N) = Y(1) in this case.  On output, YP(N) = */
/*             YP(1) and, more generally, the values and */
/*             first two derivatives of H at X(1) agree with */
/*             those at X(N).  If H(x) is one of the compo- */
/*             nents of a parametric curve, this option may */
/*             be used to obtained a closed curve.  If PER = */
/*             FALSE, H satisfies natural end conditions: */
/*             zero second derivatives at X(1) and X(N). */

/*       UNIFRM = Logical variable with value TRUE if and */
/*                only if constant (uniform) tension is to be */
/*                used.  The tension factor must be input in */
/*                SIGMA(1) in this case and must be in the */
/*                range 0 to 85.  If SIGMA(1) = 0, H(x) is */
/*                a cubic spline, and as SIGMA increases, */
/*                H(x) approaches piecewise linear.  If */
/*                UNIFRM = FALSE, tension factors are chosen */
/*                (by SIGS) to preserve local monotonicity */
/*                and convexity of the data.  This may re- */
/*                sult in a better fit than the case of */
/*                uniform tension, but requires an iteration */
/*                on calls to SMCRV and SIGS. */

/*       W = Array of length N containing positive weights */
/*           associated with the data values.  The recommend- */
/*           ed value of W(I) is 1/DY**2, where DY is the */
/*           standard deviation associated with Y(I).  If */
/*           nothing is known about the errors in Y, a con- */
/*           stant (estimated value) should be used for DY. */
/*           If PER = TRUE, it is assumed that W(N) = W(1). */

/*       SM = Positive parameter specifying an upper bound on */
/*            Q2(YS), where Q2(YS) is the weighted sum of */
/*            squares of deviations from the data (differ- */
/*            ences between YS and Y).  H(x) is linear (and */
/*            Q2 is minimized) if SM is sufficiently large */
/*            that the constraint is not active.  It is */
/*            recommended that SM satisfy N-SQRT(2N) .LE. SM */
/*            .LE. N+SQRT(2N) and SM = N is reasonable if */
/*            W(I) = 1/DY**2. */

/*       SMTOL = Parameter in the range (0,1) specifying the */
/*               relative error allowed in satisfying the */
/*               constraint:  the constraint is assumed to */
/*               be satisfied if SM*(1-SMTOL) .LE. Q2 .LE. */
/*               SM*(1+SMTOL).  A reasonable value for SMTOL */
/*               is SQRT(2/N). */

/*       LWK = Length of work space WK: */
/*             LWK .GE. 6N   if PER=FALSE  and  UNIFRM=TRUE */
/*             LWK .GE. 7N   if PER=FALSE  and  UNIFRM=FALSE */
/*             LWK .GE. 10N  if PER=TRUE   and  UNIFRM=TRUE */
/*             LWK .GE. 11N  if PER=TRUE   and  UNIFRM=FALSE */

/* The above parameters are not altered by this routine. */

/*       WK = Array of length at least LWK to be used as */
/*            temporary work space. */

/*       SIGMA = Array of length .GE. N-1 containing a ten- */
/*               sion factor (0 to 85) in the first position */
/*               if UNIFRM = TRUE. */

/*       YS = Array of length .GE. N. */

/*       YP = Array of length .GE. N. */

/* On output: */

/*       WK = Array containing convergence parameters in the */
/*            first two locations if NIT > 0: */
/*            WK(1) = Maximum relative change in a component */
/*                    of YS on the last iteration. */
/*            WK(2) = Maximum relative change in a component */
/*                    of SIGMA on the last iteration. */

/*       SIGMA = Array containing tension factors.  SIGMA(I) */
/*               is associated with interval (X(I),X(I+1)) */
/*               for I = 1,...,N-1.  SIGMA is not altered if */
/*               N is invalid or -4 < IER < -1, and SIGMA is */
/*               constant if IER = -1 (and N is valid) or */
/*               IER = -4. */

/*       YS = Array of length N containing values of H at the */
/*            abscissae.  YS(N) = YS(1) if PER = TRUE.  YS is */
/*            not altered if IER < 0. */

/*       YP = Array of length N containing first derivative */
/*            values of H at the abscissae.  YP(N) = YP(1) */
/*            if PER = TRUE.  YP is not altered if IER < 0. */

/*       NIT = Number of iterations (calls to SIGS).  NIT = 0 */
/*             if IER < 0 or UNIFRM = TRUE.  If NIT > 0, */
/*             NIT+1 calls to SMCRV were employed. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered and the */
/*                     constraint is active:  Q2(YS) is ap- */
/*                     proximately equal to SM. */
/*             IER = 1 if no errors were encountered but the */
/*                     constraint is not active:  YS and YP */
/*                     are the values and derivatives of the */
/*                     linear function (constant function if */
/*                     PERIOD = TRUE) which minimizes Q2, and */
/*                     Q1 = 0 (refer to SMCRV). */
/*             IER = -1 if N, W, SM, or SMTOL is outside its */
/*                      valid range. */
/*             IER = -2 if LWK is too small. */
/*             IER = -3 if UNIFRM = TRUE and SIGMA(1) is out- */
/*                      side its valid range. */
/*             IER = -4 if the abscissae X are not strictly */
/*                      increasing. */

/* Modules required by TSPSS:  B2TRI or B2TRIP, SIGS, SMCRV, */
/*                               SNHCSH, STORE, YPCOEF */

/* Intrinsic functions called by TSPSS:  ABS, MAX */

/* *********************************************************** */


    /* Parameter adjustments */
    --yp;
    --ys;
    --sigma;
    --w;
    --y;
    --x;
    --wk;

    /* Function Body */

/* Convergence parameters: */

/*   STOL = Absolute tolerance for SIGS */
/*   MAXIT = Maximum number of SMCRV/SIGS iterations */
/*   DYSTOL = Bound on the maximum relative change in a */
/*              component of YS defining convergence of */
/*              the SMCRV/SIGS iteration when UNIFRM = FALSE */


/* Initialize NIT, and test for invalid input parameters LWK */
/*   and SIGMA(1). */

    *nit = 0;
    nn = *n;
    nm1 = nn - 1;
    if (nn < 2 || *per && nn < 3) {
	goto L11;
    }
    if (*unifrm) {
	if (*lwk < nn * 6 || *per && *lwk < nn * 10) {
	    goto L12;
	}
	sig = sigma[1];
	if (sig < 0. || sig > sbig) {
	    goto L13;
	}
    } else {
	if (*lwk < nn * 7 || *per && *lwk < nn * 11) {
	    goto L12;
	}
	sig = 0.;
    }

/* Store uniform tension factors, or initialize SIGMA to */
/*   zeros. */

    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sigma[i__] = sig;
/* L1: */
    }

/* Compute smoothing curve for uniform tension. */

    smcrv_(&nn, &x[1], &y[1], &sigma[1], per, &w[1], sm, smtol, &wk[1], &ys[1]
	    , &yp[1], ier);
    if (*ier <= -2) {
	*ier = -4;
    }
    if (*ier < 0 || *unifrm) {
	return 0;
    }

/*   Iterate on calls to SIGS and SMCRV.  The first N-1 WK */
/*     locations are used to store the function values YS */
/*     from the previous iteration. */

/*   DYS is the maximum relative change in a component of YS. */
/*   ICNT is the number of tension factors which were */
/*        increased by SIGS. */
/*   DSMAX is the maximum relative change in a component of */
/*         SIGMA. */

    i__1 = maxit;
    for (iter = 1; iter <= i__1; ++iter) {
	dys = 0.;
	i__2 = nm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    wk[i__] = ys[i__];
/* L2: */
	}
	sigs_(&nn, &x[1], &y[1], &yp[1], &stol, &sigma[1], &dsmax, &icnt);
	smcrv_(&nn, &x[1], &y[1], &sigma[1], per, &w[1], sm, smtol, &wk[nn], &
		ys[1], &yp[1], &ierr);
	i__2 = nm1;
	for (i__ = 2; i__ <= i__2; ++i__) {
	    e = (d__1 = ys[i__] - wk[i__], abs(d__1));
	    if (wk[i__] != 0.) {
		e /= (d__1 = wk[i__], abs(d__1));
	    }
	    dys = max(dys,e);
/* L3: */
	}
	if (icnt == 0 || dys <= dystol) {
	    goto L5;
	}
/* L4: */
    }
    iter = maxit;

/* No error encountered. */

L5:
    wk[1] = dys;
    wk[2] = dsmax;
    *nit = iter;
    *ier = ierr;
    return 0;

/* Invalid input parameter N, W, SM, or SMTOL. */

L11:
    *ier = -1;
    return 0;

/* LWK too small. */

L12:
    *ier = -2;
    return 0;

/* UNIFRM = TRUE and SIGMA(1) outside its valid range. */

L13:
    *ier = -3;
    return 0;
} /* tspss_ */

/* Subroutine */ int tsval1_(integer *n, doublereal *x, doublereal *y, 
	doublereal *yp, doublereal *sigma, integer *iflag, integer *ne, 
	doublereal *te, doublereal *v, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer iflg;
    extern doublereal hval_(doublereal *, integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, integer *);
    static integer nval, ierr, i__;
    extern doublereal hpval_(doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static integer nx;
    extern doublereal hppval_(doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/10/92 */

/*   This subroutine evaluates a Hermite interpolatory ten- */
/* sion spline H or its first or second derivative at a set */
/* of points TE. */

/*   Note that a large tension factor in SIGMA may cause */
/* underflow.  The result is assumed to be zero.  If not the */
/* default, this may be specified by either a compiler option */
/* or operating system option. */

/* On input: */

/*       N = Number of data points.  N .GE. 2. */

/*       X = Array of length N containing the abscissae. */
/*           These must be in strictly increasing order: */
/*           X(I) < X(I+1) for I = 1,...,N-1. */

/*       Y = Array of length N containing data values or */
/*           function values returned by Subroutine SMCRV. */
/*           Y(I) = H(X(I)) for I = 1,...,N. */

/*       YP = Array of length N containing first deriva- */
/*            tives.  YP(I) = HP(X(I)) for I = 1,...,N, where */
/*            HP denotes the derivative of H. */

/*       SIGMA = Array of length N-1 containing tension fac- */
/*               tors whose absolute values determine the */
/*               balance between cubic and linear in each */
/*               interval.  SIGMA(I) is associated with int- */
/*               erval (I,I+1) for I = 1,...,N-1. */

/*       IFLAG = Output option indicator: */
/*               IFLAG = 0 if values of H are to be computed. */
/*               IFLAG = 1 if first derivative values are to */
/*                         be computed. */
/*               IFLAG = 2 if second derivative values are to */
/*                         be computed. */

/*       NE = Number of evaluation points.  NE > 0. */

/*       TE = Array of length NE containing the evaluation */
/*            points.  The sequence should be strictly in- */
/*            creasing for maximum efficiency.  Extrapolation */
/*            is performed if a point is not in the interval */
/*            [X(1),X(N)]. */

/* The above parameters are not altered by this routine. */

/*       V = Array of length at least NE. */

/* On output: */

/*       V = Array of function, first derivative, or second */
/*           derivative values at the evaluation points un- */
/*           less IER < 0.  If IER = -1, V is not altered. */
/*           If IER = -2, V may be only partially defined. */

/*       IER = Error indicator: */
/*             IER = 0  if no errors were encountered and */
/*                      no extrapolation occurred. */
/*             IER > 0  if no errors were encountered but */
/*                      extrapolation was required at IER */
/*                      points. */
/*             IER = -1 if N < 2, IFLAG < 0, IFLAG > 2, or */
/*                      NE < 1. */
/*             IER = -2 if the abscissae are not in strictly */
/*                      increasing order.  (This error will */
/*                      not necessarily be detected.) */

/* Modules required by TSVAL1:  HPPVAL, HPVAL, HVAL, INTRVL, */
/*                                SNHCSH */

/* *********************************************************** */


    /* Parameter adjustments */
    --sigma;
    --yp;
    --y;
    --x;
    --v;
    --te;

    /* Function Body */
    iflg = *iflag;
    nval = *ne;

/* Test for invalid input. */

    if (*n < 2 || iflg < 0 || iflg > 2 || nval < 1) {
	goto L2;
    }

/* Initialize the number of extrapolation points NX and */
/*   loop on evaluation points. */

    nx = 0;
    i__1 = nval;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iflg == 0) {
	    v[i__] = hval_(&te[i__], n, &x[1], &y[1], &yp[1], &sigma[1], &
		    ierr);
	} else if (iflg == 1) {
	    v[i__] = hpval_(&te[i__], n, &x[1], &y[1], &yp[1], &sigma[1], &
		    ierr);
	} else {
	    v[i__] = hppval_(&te[i__], n, &x[1], &y[1], &yp[1], &sigma[1], &
		    ierr);
	}
	if (ierr < 0) {
	    goto L3;
	}
	nx += ierr;
/* L1: */
    }

/* No errors encountered. */

    *ier = nx;
    return 0;

/* N, IFLAG, or NE is outside its valid range. */

L2:
    *ier = -1;
    return 0;

/* X is not strictly increasing. */

L3:
    *ier = -2;
    return 0;
} /* tsval1_ */

/* Subroutine */ int tsval2_(integer *n, doublereal *t, doublereal *x, 
	doublereal *y, doublereal *xp, doublereal *yp, doublereal *sigma, 
	integer *iflag, integer *ne, doublereal *te, doublereal *vx, 
	doublereal *vy, integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer iflg;
    extern doublereal hval_(doublereal *, integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, integer *);
    static integer nval, ierr, i__;
    extern doublereal hpval_(doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static integer nx;
    extern doublereal hppval_(doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/10/92 */

/*   This subroutine returns values or derivatives of a pair */
/* of Hermite interpolatory tension splines H1 and H2 which */
/* form the components of a parametric planar curve C(t) = */
/* (H1(t),H2(t)).  Refer to Subroutines TSPBP and TSPSP. */

/*   Note that a large tension factor in SIGMA may cause */
/* underflow.  The result is assumed to be zero.  If not the */
/* default, this may be specified by either a compiler option */
/* or operating system option. */

/* On input: */

/*       N = Number of data points.  N .GE. 2. */

/*       T = Array of length N containing a strictly in- */
/*           creasing sequence of abscissae (parameter */
/*           values).  Refer to Subroutine ARCL2D. */

/*       X = Array of length N containing data values or */
/*           function values returned by Subroutine SMCRV. */
/*           X(I) = H1(T(I)) for I = 1,...,N. */

/*       Y = Array of length N containing data values or */
/*           function values returned by Subroutine SMCRV. */
/*           Y(I) = H2(T(I)) for I = 1,...,N. */

/*       XP = Array of length N containing first deriva- */
/*            tives.  XP(I) = H1P(T(I)) for I = 1,...,N, */
/*            where H1P denotes the derivative of H1. */

/*       YP = Array of length N containing first deriva- */
/*            tives.  YP(I) = H2P(T(I)) for I = 1,...,N, */
/*            where H2P denotes the derivative of H2. */

/*   Note that C(T(I)) = (X(I),Y(I)) and CP(T(I)) = (XP(I), */
/* YP(I)), I = 1,...,N, are data (control) points and deriva- */
/* tive (velocity) vectors, respectively. */

/*       SIGMA = Array of length N-1 containing tension fac- */
/*               tors whose absolute values determine the */
/*               balance between cubic and linear in each */
/*               interval.  SIGMA(I) is associated with int- */
/*               erval (I,I+1) for I = 1,...,N-1. */

/*       IFLAG = Output option indicator: */
/*               IFLAG = 0 if values of H1 and H2 (points on */
/*                         the curve) are to be computed. */
/*               IFLAG = 1 if first derivative vectors are to */
/*                         be computed.  Unit tangent vectors */
/*                         can be obtained by normalizing */
/*                         these to unit vectors. */
/*               IFLAG = 2 if second derivative (accelera- */
/*                         tion) vectors are to be computed. */
/*                         Given a unit tangent vector U and */
/*                         a second derivative vector V, the */
/*                         corresponding curvature vector */
/*                         can be computed as the cross */
/*                         product U X V X U. */

/*       NE = Number of evaluation points.  NE > 0. */

/*       TE = Array of length NE containing the evaluation */
/*            points.  The sequence should be strictly in- */
/*            creasing for maximum efficiency.  Extrapolation */
/*            is performed if a point is not in the interval */
/*            [T(1),T(N)]. */

/* The above parameters are not altered by this routine. */

/*       VX,VY = Arrays of length at least NE. */

/* On output: */

/*       VX,VY = Arrays containing values, first derivatives, */
/*               or second derivatives of H1 and H2, respec- */
/*               tively, at the evaluation points (unless */
/*               IER < 0).  If IER = -1, VX and VY are not */
/*               altered.  If IER = -2, VX and VY may be only */
/*               partially defined. */

/*       IER = Error indicator: */
/*             IER = 0  if no errors were encountered and */
/*                      no extrapolation occurred. */
/*             IER > 0  if no errors were encountered but */
/*                      extrapolation was required at IER */
/*                      points. */
/*             IER = -1 if N < 2, IFLAG < 0, IFLAG > 2, or */
/*                      NE < 1. */
/*             IER = -2 if the abscissae are not in strictly */
/*                      increasing order.  (This error will */
/*                      not necessarily be detected.) */

/* Modules required by TSVAL2:  HPPVAL, HPVAL, HVAL, INTRVL, */
/*                                SNHCSH */

/* *********************************************************** */


    /* Parameter adjustments */
    --sigma;
    --yp;
    --xp;
    --y;
    --x;
    --t;
    --vy;
    --vx;
    --te;

    /* Function Body */
    iflg = *iflag;
    nval = *ne;

/* Test for invalid input. */

    if (*n < 2 || iflg < 0 || iflg > 2 || nval < 1) {
	goto L2;
    }

/* Initialize the number of extrapolation points NX and */
/*   loop on evaluation points. */

    nx = 0;
    i__1 = nval;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iflg == 0) {
	    vx[i__] = hval_(&te[i__], n, &t[1], &x[1], &xp[1], &sigma[1], &
		    ierr);
	    vy[i__] = hval_(&te[i__], n, &t[1], &y[1], &yp[1], &sigma[1], &
		    ierr);
	} else if (iflg == 1) {
	    vx[i__] = hpval_(&te[i__], n, &t[1], &x[1], &xp[1], &sigma[1], &
		    ierr);
	    vy[i__] = hpval_(&te[i__], n, &t[1], &y[1], &yp[1], &sigma[1], &
		    ierr);
	} else {
	    vx[i__] = hppval_(&te[i__], n, &t[1], &x[1], &xp[1], &sigma[1], &
		    ierr);
	    vy[i__] = hppval_(&te[i__], n, &t[1], &y[1], &yp[1], &sigma[1], &
		    ierr);
	}
	if (ierr < 0) {
	    goto L3;
	}
	nx += ierr;
/* L1: */
    }

/* No errors encountered. */

    *ier = nx;
    return 0;

/* N, IFLAG, or NE is outside its valid range. */

L2:
    *ier = -1;
    return 0;

/* T is not strictly increasing. */

L3:
    *ier = -2;
    return 0;
} /* tsval2_ */

/* Subroutine */ int tsval3_(integer *n, doublereal *t, doublereal *x, 
	doublereal *y, doublereal *z__, doublereal *xp, doublereal *yp, 
	doublereal *zp, doublereal *sigma, integer *iflag, integer *ne, 
	doublereal *te, doublereal *vx, doublereal *vy, doublereal *vz, 
	integer *ier)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer iflg;
    extern doublereal hval_(doublereal *, integer *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, integer *);
    static integer nval, ierr, i__;
    extern doublereal hpval_(doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static integer nx;
    extern doublereal hppval_(doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/10/92 */

/*   This subroutine returns values or derivatives of three */
/* Hermite interpolatory tension splines H1, H2, and H3 which */
/* form the components of a parametric space curve C(t) = */
/* (H1(t),H2(t),H3(t)).  Refer to Subroutines TSPBP and */
/* TSPSP. */

/*   Note that a large tension factor in SIGMA may cause */
/* underflow.  The result is assumed to be zero.  If not the */
/* default, this may be specified by either a compiler option */
/* or operating system option. */

/* On input: */

/*       N = Number of data points.  N .GE. 2. */

/*       T = Array of length N containing a strictly in- */
/*           creasing sequence of abscissae (parameter */
/*           values).  Refer to Subroutine ARCL3D. */

/*       X = Array of length N containing data values or */
/*           function values returned by Subroutine SMCRV. */
/*           X(I) = H1(T(I)) for I = 1,...,N. */

/*       Y = Array of length N containing data values or */
/*           function values returned by Subroutine SMCRV. */
/*           Y(I) = H2(T(I)) for I = 1,...,N. */

/*       Z = Array of length N containing data values or */
/*           function values returned by Subroutine SMCRV. */
/*           Z(I) = H3(T(I)) for I = 1,...,N. */

/*       XP = Array of length N containing first deriva- */
/*            tives.  XP(I) = H1P(T(I)) for I = 1,...,N, */
/*            where H1P denotes the derivative of H1. */

/*       YP = Array of length N containing first deriva- */
/*            tives.  YP(I) = H2P(T(I)) for I = 1,...,N, */
/*            where H2P denotes the derivative of H2. */

/*       ZP = Array of length N containing first deriva- */
/*            tives.  ZP(I) = H3P(T(I)) for I = 1,...,N, */
/*            where H3P denotes the derivative of H3. */

/*   Note that C(T(I)) = (X(I),Y(I),Z(I)) and CP(T(I)) = */
/* (XP(I),YP(I),ZP(I)), I = 1,...,N, are data (control) */
/* points and derivative (velocity) vectors, respectively. */

/*       SIGMA = Array of length N-1 containing tension fac- */
/*               tors whose absolute values determine the */
/*               balance between cubic and linear in each */
/*               interval.  SIGMA(I) is associated with int- */
/*               erval (I,I+1) for I = 1,...,N-1. */

/*       IFLAG = Output option indicator: */
/*               IFLAG = 0 if values of H1, H2, and H3 */
/*                         (points on the curve) are to be */
/*                         computed. */
/*               IFLAG = 1 if first derivative vectors are to */
/*                         be computed.  Unit tangent vectors */
/*                         can be obtained by normalizing */
/*                         these to unit vectors. */
/*               IFLAG = 2 if second derivative (accelera- */
/*                         tion) vectors are to be computed. */
/*                         Given a unit tangent vector U and */
/*                         a second derivative vector V, the */
/*                         corresponding curvature vector */
/*                         can be computed as the cross */
/*                         product U X V X U. */

/*       NE = Number of evaluation points.  NE > 0. */

/*       TE = Array of length NE containing the evaluation */
/*            points.  The sequence should be strictly in- */
/*            creasing for maximum efficiency.  Extrapolation */
/*            is performed if a point is not in the interval */
/*            [T(1),T(N)]. */

/* The above parameters are not altered by this routine. */

/*       VX,VY,VZ = Arrays of length at least NE. */

/* On output: */

/*       VX,VY,VZ = Arrays containing values, first deriva- */
/*                  tives, or second derivatives of H1, H2, */
/*                  and H3, respectively, at the evaluation */
/*                  points (unless IER < 0).  If IER = -1, */
/*                  VX, VY, and VZ are not altered.  If IER */
/*                  = -2, VX, VY, and VZ may be only partial- */
/*                  ly defined. */

/*       IER = Error indicator: */
/*             IER = 0  if no errors were encountered and */
/*                      no extrapolation occurred. */
/*             IER > 0  if no errors were encountered but */
/*                      extrapolation was required at IER */
/*                      points. */
/*             IER = -1 if N < 2, IFLAG < 0, IFLAG > 2, or */
/*                      NE < 1. */
/*             IER = -2 if the abscissae are not in strictly */
/*                      increasing order.  (This error will */
/*                      not necessarily be detected.) */

/* Modules required by TSVAL3:  HPPVAL, HPVAL, HVAL, INTRVL, */
/*                                SNHCSH */

/* *********************************************************** */


    /* Parameter adjustments */
    --sigma;
    --zp;
    --yp;
    --xp;
    --z__;
    --y;
    --x;
    --t;
    --vz;
    --vy;
    --vx;
    --te;

    /* Function Body */
    iflg = *iflag;
    nval = *ne;

/* Test for invalid input. */

    if (*n < 2 || iflg < 0 || iflg > 2 || nval < 1) {
	goto L2;
    }

/* Initialize the number of extrapolation points NX and */
/*   loop on evaluation points. */

    nx = 0;
    i__1 = nval;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iflg == 0) {
	    vx[i__] = hval_(&te[i__], n, &t[1], &x[1], &xp[1], &sigma[1], &
		    ierr);
	    vy[i__] = hval_(&te[i__], n, &t[1], &y[1], &yp[1], &sigma[1], &
		    ierr);
	    vz[i__] = hval_(&te[i__], n, &t[1], &z__[1], &zp[1], &sigma[1], &
		    ierr);
	} else if (iflg == 1) {
	    vx[i__] = hpval_(&te[i__], n, &t[1], &x[1], &xp[1], &sigma[1], &
		    ierr);
	    vy[i__] = hpval_(&te[i__], n, &t[1], &y[1], &yp[1], &sigma[1], &
		    ierr);
	    vz[i__] = hpval_(&te[i__], n, &t[1], &z__[1], &zp[1], &sigma[1], &
		    ierr);
	} else {
	    vx[i__] = hppval_(&te[i__], n, &t[1], &x[1], &xp[1], &sigma[1], &
		    ierr);
	    vy[i__] = hppval_(&te[i__], n, &t[1], &y[1], &yp[1], &sigma[1], &
		    ierr);
	    vz[i__] = hppval_(&te[i__], n, &t[1], &z__[1], &zp[1], &sigma[1], 
		    &ierr);
	}
	if (ierr < 0) {
	    goto L3;
	}
	nx += ierr;
/* L1: */
    }

/* No errors encountered. */

    *ier = nx;
    return 0;

/* N, IFLAG, or NE is outside its valid range. */

L2:
    *ier = -1;
    return 0;

/* T is not strictly increasing. */

L3:
    *ier = -2;
    return 0;
} /* tsval3_ */

/* Subroutine */ int ypc1_(integer *n, doublereal *x, doublereal *y, 
	doublereal *yp, integer *ier)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal asim1, dxim1;
    static integer i__;
    static doublereal t, s2, si;
    static integer nm1;
    static doublereal dx2, asi, dxi, sgn, sim1;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/10/92 */

/*   This subroutine employs a three-point quadratic interpo- */
/* lation method to compute local derivative estimates YP */
/* associated with a set of data points.  The interpolation */
/* formula is the monotonicity-constrained parabolic method */
/* described in the reference cited below.  A Hermite int- */
/* erpolant of the data values and derivative estimates pre- */
/* serves monotonicity of the data.  Linear interpolation is */
/* used if N = 2.  The method is invariant under a linear */
/* scaling of the coordinates but is not additive. */

/* On input: */

/*       N = Number of data points.  N .GE. 2. */

/*       X = Array of length N containing a strictly in- */
/*           creasing sequence of abscissae:  X(I) < X(I+1) */
/*           for I = 1,...,N-1. */

/*       Y = Array of length N containing data values asso- */
/*           ciated with the abscissae. */

/* Input parameters are not altered by this routine. */

/* On output: */

/*       YP = Array of length N containing estimated deriv- */
/*            atives at the abscissae unless IER .NE. 0. */
/*            YP is not altered if IER = 1, and is only par- */
/*            tially defined if IER > 1. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if N < 2. */
/*             IER = I if X(I) .LE. X(I-1) for some I in the */
/*                     range 2,...,N. */

/* Reference:  J. M. Hyman, "Accurate Monotonicity-preserving */
/*               Cubic Interpolation",  LA-8796-MS, Los */
/*               Alamos National Lab, Feb. 1982. */

/* Modules required by YPC1:  None */

/* Intrinsic functions called by YPC1:  ABS, MAX, MIN, SIGN */

/* *********************************************************** */


    /* Parameter adjustments */
    --yp;
    --y;
    --x;

    /* Function Body */
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L2;
    }
    i__ = 1;
    dxi = x[2] - x[1];
    if (dxi <= 0.) {
	goto L3;
    }
    si = (y[2] - y[1]) / dxi;
    if (nm1 == 1) {

/* Use linear interpolation for N = 2. */

	yp[1] = si;
	yp[2] = si;
	*ier = 0;
	return 0;
    }

/* N .GE. 3.  YP(1) = S1 + DX1*(S1-S2)/(DX1+DX2) unless this */
/*   results in YP(1)*S1 .LE. 0 or abs(YP(1)) > 3*abs(S1). */

    i__ = 2;
    dx2 = x[3] - x[2];
    if (dx2 <= 0.) {
	goto L3;
    }
    s2 = (y[3] - y[2]) / dx2;
    t = si + dxi * (si - s2) / (dxi + dx2);
    if (si >= 0.) {
/* Computing MIN */
	d__1 = max(0.,t), d__2 = si * 3.;
	yp[1] = min(d__1,d__2);
    } else {
/* Computing MAX */
	d__1 = min(0.,t), d__2 = si * 3.;
	yp[1] = max(d__1,d__2);
    }

/* YP(I) = (DXIM1*SI+DXI*SIM1)/(DXIM1+DXI) subject to the */
/*   constraint that YP(I) has the sign of either SIM1 or */
/*   SI, whichever has larger magnitude, and abs(YP(I)) .LE. */
/*   3*min(abs(SIM1),abs(SI)). */

    i__1 = nm1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	dxim1 = dxi;
	dxi = x[i__ + 1] - x[i__];
	if (dxi <= 0.) {
	    goto L3;
	}
	sim1 = si;
	si = (y[i__ + 1] - y[i__]) / dxi;
	t = (dxim1 * si + dxi * sim1) / (dxim1 + dxi);
	asim1 = abs(sim1);
	asi = abs(si);
	sgn = d_sign(&c_b82, &si);
	if (asim1 > asi) {
	    sgn = d_sign(&c_b82, &sim1);
	}
	if (sgn > 0.) {
/* Computing MIN */
	    d__1 = max(0.,t), d__2 = min(asim1,asi) * 3.;
	    yp[i__] = min(d__1,d__2);
	} else {
/* Computing MAX */
	    d__1 = min(0.,t), d__2 = min(asim1,asi) * -3.;
	    yp[i__] = max(d__1,d__2);
	}
/* L1: */
    }

/* YP(N) = SNM1 + DXNM1*(SNM1-SNM2)/(DXNM2+DXNM1) subject to */
/*   the constraint that YP(N) has the sign of SNM1 and */
/*   abs(YP(N)) .LE. 3*abs(SNM1).  Note that DXI = DXNM1 and */
/*   SI = SNM1. */

    t = si + dxi * (si - sim1) / (dxim1 + dxi);
    if (si >= 0.) {
/* Computing MIN */
	d__1 = max(0.,t), d__2 = si * 3.;
	yp[*n] = min(d__1,d__2);
    } else {
/* Computing MAX */
	d__1 = min(0.,t), d__2 = si * 3.;
	yp[*n] = max(d__1,d__2);
    }
    *ier = 0;
    return 0;

/* N is outside its valid range. */

L2:
    *ier = 1;
    return 0;

/* X(I+1) .LE. X(I). */

L3:
    *ier = i__ + 1;
    return 0;
} /* ypc1_ */

/* Subroutine */ int ypc1p_(integer *n, doublereal *x, doublereal *y, 
	doublereal *yp, integer *ier)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal asim1, dxim1;
    static integer i__;
    static doublereal t, si;
    static integer nm1;
    static doublereal asi, dxi, sgn, sim1;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/10/92 */

/*   This subroutine employs a three-point quadratic interpo- */
/* lation method to compute local derivative estimates YP */
/* associated with a set of N data points (X(I),Y(I)).  It */
/* is assumed that Y(N) = Y(1), and YP(N) = YP(1) on output. */
/* Thus, a Hermite interpolant H(x) defined by the data */
/* points and derivative estimates is periodic with period */
/* X(N)-X(1).  The derivative-estimation formula is the */
/* monotonicity-constrained parabolic fit described in the */
/* reference cited below:  H(x) is monotonic in intervals in */
/* which the data is monotonic.  The method is invariant */
/* under a linear scaling of the coordinates but is not */
/* additive. */

/* On input: */

/*       N = Number of data points.  N .GE. 3. */

/*       X = Array of length N containing a strictly in- */
/*           creasing sequence of abscissae:  X(I) < X(I+1) */
/*           for I = 1,...,N-1. */

/*       Y = Array of length N containing data values asso- */
/*           ciated with the abscissae.  Y(N) is set to Y(1) */
/*           on output unless IER = 1. */

/*   Input parameters, other than Y(N), are not altered by */
/* this routine. */

/* On output: */

/*       YP = Array of length N containing estimated deriv- */
/*            atives at the abscissae unless IER .NE. 0. */
/*            YP is not altered if IER = 1, and is only par- */
/*            tially defined if IER > 1. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if N < 3. */
/*             IER = I if X(I) .LE. X(I-1) for some I in the */
/*                     range 2,...,N. */

/* Reference:  J. M. Hyman, "Accurate Monotonicity-preserving */
/*               Cubic Interpolation",  LA-8796-MS, Los */
/*               Alamos National Lab, Feb. 1982. */

/* Modules required by YPC1P:  None */

/* Intrinsic functions called by YPC1P:  ABS, MAX, MIN, SIGN */

/* *********************************************************** */


    /* Parameter adjustments */
    --yp;
    --y;
    --x;

    /* Function Body */
    nm1 = *n - 1;
    if (nm1 < 2) {
	goto L2;
    }
    y[*n] = y[1];

/* Initialize for loop on interior points. */

    i__ = 1;
    dxi = x[2] - x[1];
    if (dxi <= 0.) {
	goto L3;
    }
    si = (y[2] - y[1]) / dxi;

/* YP(I) = (DXIM1*SI+DXI*SIM1)/(DXIM1+DXI) subject to the */
/*   constraint that YP(I) has the sign of either SIM1 or */
/*   SI, whichever has larger magnitude, and abs(YP(I)) .LE. */
/*   3*min(abs(SIM1),abs(SI)). */

    i__1 = nm1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	dxim1 = dxi;
	dxi = x[i__ + 1] - x[i__];
	if (dxi <= 0.) {
	    goto L3;
	}
	sim1 = si;
	si = (y[i__ + 1] - y[i__]) / dxi;
	t = (dxim1 * si + dxi * sim1) / (dxim1 + dxi);
	asim1 = abs(sim1);
	asi = abs(si);
	sgn = d_sign(&c_b82, &si);
	if (asim1 > asi) {
	    sgn = d_sign(&c_b82, &sim1);
	}
	if (sgn > 0.) {
/* Computing MIN */
	    d__1 = max(0.,t), d__2 = min(asim1,asi) * 3.;
	    yp[i__] = min(d__1,d__2);
	} else {
/* Computing MAX */
	    d__1 = min(0.,t), d__2 = min(asim1,asi) * -3.;
	    yp[i__] = max(d__1,d__2);
	}
/* L1: */
    }

/* YP(N) = YP(1), I = 1, and IM1 = N-1. */

    dxim1 = dxi;
    dxi = x[2] - x[1];
    sim1 = si;
    si = (y[2] - y[1]) / dxi;
    t = (dxim1 * si + dxi * sim1) / (dxim1 + dxi);
    asim1 = abs(sim1);
    asi = abs(si);
    sgn = d_sign(&c_b82, &si);
    if (asim1 > asi) {
	sgn = d_sign(&c_b82, &sim1);
    }
    if (sgn > 0.) {
/* Computing MIN */
	d__1 = max(0.,t), d__2 = min(asim1,asi) * 3.;
	yp[1] = min(d__1,d__2);
    } else {
/* Computing MAX */
	d__1 = min(0.,t), d__2 = min(asim1,asi) * -3.;
	yp[1] = max(d__1,d__2);
    }
    yp[*n] = yp[1];

/* No error encountered. */

    *ier = 0;
    return 0;

/* N is outside its valid range. */

L2:
    *ier = 1;
    return 0;

/* X(I+1) .LE. X(I). */

L3:
    *ier = i__ + 1;
    return 0;
} /* ypc1p_ */

/* Subroutine */ int ypc2_(integer *n, doublereal *x, doublereal *y, 
	doublereal *sigma, integer *isl1, integer *isln, doublereal *bv1, 
	doublereal *bvn, doublereal *wk, doublereal *yp, integer *ier)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal d__;
    static integer i__;
    static doublereal s, d1, d2, r1, r2;
    static integer nn;
    static doublereal dx;
    extern doublereal endslp_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int ypcoef_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal sd1, sd2;
    static integer nm1;
    static doublereal yp1, sig, ypn;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/10/92 */

/*   This subroutine solves a linear system for a set of */
/* first derivatives YP associated with a Hermite interpola- */
/* tory tension spline H(x).  The derivatives are chosen so */
/* that H(x) has two continuous derivatives for all x and H */
/* satisfies user-specified end conditions. */

/* On input: */

/*       N = Number of data points.  N .GE. 2. */

/*       X = Array of length N containing a strictly in- */
/*           creasing sequence of abscissae:  X(I) < X(I+1) */
/*           for I = 1,...,N-1. */

/*       Y = Array of length N containing data values asso- */
/*           ciated with the abscissae.  H(X(I)) = Y(I) for */
/*           I = 1,...,N. */

/*       SIGMA = Array of length N-1 containing tension */
/*               factors.  SIGMA(I) is associated with inter- */
/*               val (X(I),X(I+1)) for I = 1,...,N-1.  If */
/*               SIGMA(I) = 0, H is the Hermite cubic interp- */
/*               olant of the data values and computed deriv- */
/*               atives at X(I) and X(I+1), and if all */
/*               tension factors are zero, H is the C-2 cubic */
/*               spline interpolant which satisfies the end */
/*               conditions. */

/*       ISL1 = Option indicator for the condition at X(1): */
/*              ISL1 = 0 if YP(1) is to be estimated inter- */
/*                       nally by a constrained parabolic */
/*                       fit to the first three points. */
/*                       This is identical to the method used */
/*                       by Subroutine YPC1.  BV1 is not used */
/*                       in this case. */
/*              ISL1 = 1 if the first derivative of H at X(1) */
/*                       is specified by BV1. */
/*              ISL1 = 2 if the second derivative of H at */
/*                       X(1) is specified by BV1. */
/*              ISL1 = 3 if YP(1) is to be estimated inter- */
/*                       nally from the derivative of the */
/*                       tension spline (using SIGMA(1)) */
/*                       which interpolates the first three */
/*                       data points and has third derivative */
/*                       equal to zero at X(1).  Refer to */
/*                       ENDSLP.  BV1 is not used in this */
/*                       case. */

/*       ISLN = Option indicator for the condition at X(N): */
/*              ISLN = 0 if YP(N) is to be estimated inter- */
/*                       nally by a constrained parabolic */
/*                       fit to the last three data points. */
/*                       This is identical to the method used */
/*                       by Subroutine YPC1.  BVN is not used */
/*                       in this case. */
/*              ISLN = 1 if the first derivative of H at X(N) */
/*                       is specified by BVN. */
/*              ISLN = 2 if the second derivative of H at */
/*                       X(N) is specified by BVN. */
/*              ISLN = 3 if YP(N) is to be estimated inter- */
/*                       nally from the derivative of the */
/*                       tension spline (using SIGMA(N-1)) */
/*                       which interpolates the last three */
/*                       data points and has third derivative */
/*                       equal to zero at X(N).  Refer to */
/*                       ENDSLP.  BVN is not used in this */
/*                       case. */

/*       BV1,BVN = Boundary values or dummy parameters as */
/*                 defined by ISL1 and ISLN. */

/* The above parameters are not altered by this routine. */

/*       WK = Array of length at least N-1 to be used as */
/*            temporary work space. */

/*       YP = Array of length .GE. N. */

/* On output: */

/*       YP = Array containing derivatives of H at the */
/*            abscissae.  YP is not defined if IER .NE. 0. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if N, ISL1, or ISLN is outside its */
/*                     valid range. */
/*             IER = I if X(I) .LE. X(I-1) for some I in the */
/*                     range 2,...,N. */

/* Modules required by YPC2:  ENDSLP, SNHCSH, YPCOEF */

/* Intrinsic function called by YPC2:  ABS */

/* *********************************************************** */


    /* Parameter adjustments */
    --yp;
    --wk;
    --sigma;
    --y;
    --x;

    /* Function Body */
    nn = *n;
    if (nn < 2 || *isl1 < 0 || *isl1 > 3 || *isln < 0 || *isln > 3) {
	goto L3;
    }
    nm1 = nn - 1;

/* Set YP1 and YPN to the endpoint values. */

    if (*isl1 == 0) {
	if (nn > 2) {
	    yp1 = endslp_(&x[1], &x[2], &x[3], &y[1], &y[2], &y[3], &c_b192);
	}
    } else if (*isl1 != 3) {
	yp1 = *bv1;
    } else {
	if (nn > 2) {
	    yp1 = endslp_(&x[1], &x[2], &x[3], &y[1], &y[2], &y[3], &sigma[1])
		    ;
	}
    }
    if (*isln == 0) {
	if (nn > 2) {
	    ypn = endslp_(&x[nn], &x[nm1], &x[nn - 2], &y[nn], &y[nm1], &y[nn 
		    - 2], &c_b192);
	}
    } else if (*isln != 3) {
	ypn = *bvn;
    } else {
	if (nn > 2) {
	    ypn = endslp_(&x[nn], &x[nm1], &x[nn - 2], &y[nn], &y[nm1], &y[nn 
		    - 2], &sigma[nm1]);
	}
    }

/* Solve the symmetric positive-definite tridiagonal linear */
/*   system.  The forward elimination step consists of div- */
/*   iding each row by its diagonal entry, then introducing a */
/*   zero below the diagonal.  This requires saving only the */
/*   superdiagonal (in WK) and the right hand side (in YP). */

    i__ = 1;
    dx = x[2] - x[1];
    if (dx <= 0.) {
	goto L4;
    }
    s = (y[2] - y[1]) / dx;
    if (nn == 2) {
	if (*isl1 == 0 || *isl1 == 3) {
	    yp1 = s;
	}
	if (*isln == 0 || *isln == 3) {
	    ypn = s;
	}
    }

/* Begin forward elimination. */

    sig = abs(sigma[1]);
    ypcoef_(&sig, &dx, &d1, &sd1);
    r1 = (sd1 + d1) * s;
    wk[1] = 0.;
    yp[1] = yp1;
    if (*isl1 == 2) {
	wk[1] = sd1 / d1;
	yp[1] = (r1 - yp1) / d1;
    }
    i__1 = nm1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	dx = x[i__ + 1] - x[i__];
	if (dx <= 0.) {
	    goto L4;
	}
	s = (y[i__ + 1] - y[i__]) / dx;
	sig = (d__1 = sigma[i__], abs(d__1));
	ypcoef_(&sig, &dx, &d2, &sd2);
	r2 = (sd2 + d2) * s;
	d__ = d1 + d2 - sd1 * wk[i__ - 1];
	wk[i__] = sd2 / d__;
	yp[i__] = (r1 + r2 - sd1 * yp[i__ - 1]) / d__;
	d1 = d2;
	sd1 = sd2;
	r1 = r2;
/* L1: */
    }
    d__ = d1 - sd1 * wk[nm1];
    yp[nn] = ypn;
    if (*isln == 2) {
	yp[nn] = (r1 + ypn - sd1 * yp[nm1]) / d__;
    }

/* Back substitution: */

    for (i__ = nm1; i__ >= 1; --i__) {
	yp[i__] -= wk[i__] * yp[i__ + 1];
/* L2: */
    }
    *ier = 0;
    return 0;

/* Invalid integer input parameter. */

L3:
    *ier = 1;
    return 0;

/* Abscissae out of order or duplicate points. */

L4:
    *ier = i__ + 1;
    return 0;
} /* ypc2_ */

/* Subroutine */ int ypc2p_(integer *n, doublereal *x, doublereal *y, 
	doublereal *sigma, doublereal *wk, doublereal *yp, integer *ier)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static doublereal sdnm1, d__, ypnm1;
    static integer i__;
    static doublereal s, d1, d2, r1, r2;
    static integer nn;
    static doublereal dx;
    extern /* Subroutine */ int ypcoef_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal sd1, sd2;
    static integer nm1, nm2, nm3, np1;
    static doublereal din, sig;
    static integer npi;
    static doublereal dnm1, rnm1;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   06/10/92 */

/*   This subroutine solves a linear system for a set of */
/* first derivatives YP associated with a Hermite interpola- */
/* tory tension spline H(x).  The derivatives are chosen so */
/* that H(x) has two continuous derivatives for all x, and H */
/* satisfies periodic end conditions:  first and second der- */
/* ivatives of H at X(1) agree with those at X(N), and thus */
/* the length of a period is X(N) - X(1).  It is assumed that */
/* Y(N) = Y(1), and Y(N) is not referenced. */

/* On input: */

/*       N = Number of data points.  N .GE. 3. */

/*       X = Array of length N containing a strictly in- */
/*           creasing sequence of abscissae:  X(I) < X(I+1) */
/*           for I = 1,...,N-1. */

/*       Y = Array of length N containing data values asso- */
/*           ciated with the abscissae.  H(X(I)) = Y(I) for */
/*           I = 1,...,N. */

/*       SIGMA = Array of length N-1 containing tension */
/*               factors.  SIGMA(I) is associated with inter- */
/*               val (X(I),X(I+1)) for I = 1,...,N-1.  If */
/*               SIGMA(I) = 0, H is the Hermite cubic interp- */
/*               olant of the data values and computed deriv- */
/*               atives at X(I) and X(I+1), and if all */
/*               tension factors are zero, H is the C-2 cubic */
/*               spline interpolant which satisfies the end */
/*               conditions. */

/* The above parameters are not altered by this routine. */

/*       WK = Array of length at least 2N-2 to be used as */
/*            temporary work space. */

/*       YP = Array of length .GE. N. */

/* On output: */

/*       YP = Array containing derivatives of H at the */
/*            abscissae.  YP is not defined if IER .NE. 0. */

/*       IER = Error indicator: */
/*             IER = 0 if no errors were encountered. */
/*             IER = 1 if N is outside its valid range. */
/*             IER = I if X(I) .LE. X(I-1) for some I in the */
/*                     range 2,...,N. */

/* Modules required by YPC2P:  SNHCSH, YPCOEF */

/* Intrinsic function called by YPC2P:  ABS */

/* *********************************************************** */


    /* Parameter adjustments */
    --yp;
    --sigma;
    --y;
    --x;
    --wk;

    /* Function Body */
    nn = *n;
    if (nn < 3) {
	goto L4;
    }
    nm1 = nn - 1;
    nm2 = nn - 2;
    nm3 = nn - 3;
    np1 = nn + 1;

/* The system is order N-1, symmetric, positive-definite, and */
/*   tridiagonal except for nonzero elements in the upper */
/*   right and lower left corners.  The forward elimination */
/*   step zeros the subdiagonal and divides each row by its */
/*   diagonal entry for the first N-2 rows.  The superdiago- */
/*   nal is stored in WK(I), the negative of the last column */
/*   (fill-in) in WK(N+I), and the right hand side in YP(I) */
/*   for I = 1,...,N-2. */

    i__ = nm1;
    dx = x[nn] - x[nm1];
    if (dx <= 0.) {
	goto L5;
    }
    s = (y[1] - y[nm1]) / dx;
    sig = (d__1 = sigma[nm1], abs(d__1));
    ypcoef_(&sig, &dx, &dnm1, &sdnm1);
    rnm1 = (sdnm1 + dnm1) * s;
    i__ = 1;
    dx = x[2] - x[1];
    if (dx <= 0.) {
	goto L5;
    }
    s = (y[2] - y[1]) / dx;
    sig = abs(sigma[1]);
    ypcoef_(&sig, &dx, &d1, &sd1);
    r1 = (sd1 + d1) * s;
    d__ = dnm1 + d1;
    wk[1] = sd1 / d__;
    wk[np1] = -sdnm1 / d__;
    yp[1] = (rnm1 + r1) / d__;
    i__1 = nm2;
    for (i__ = 2; i__ <= i__1; ++i__) {
	dx = x[i__ + 1] - x[i__];
	if (dx <= 0.) {
	    goto L5;
	}
	s = (y[i__ + 1] - y[i__]) / dx;
	sig = (d__1 = sigma[i__], abs(d__1));
	ypcoef_(&sig, &dx, &d2, &sd2);
	r2 = (sd2 + d2) * s;
	d__ = d1 + d2 - sd1 * wk[i__ - 1];
	din = 1. / d__;
	wk[i__] = sd2 * din;
	npi = nn + i__;
	wk[npi] = -sd1 * wk[npi - 1] * din;
	yp[i__] = (r1 + r2 - sd1 * yp[i__ - 1]) * din;
	sd1 = sd2;
	d1 = d2;
	r1 = r2;
/* L1: */
    }

/* The backward elimination step zeros the superdiagonal */
/*   (first N-3 elements).  WK(I) and YP(I) are overwritten */
/*   with the negative of the last column and the new right */
/*   hand side, respectively, for I = N-2, N-3, ..., 1. */

    npi = nn + nm2;
    wk[nm2] = wk[npi] - wk[nm2];
    for (i__ = nm3; i__ >= 1; --i__) {
	yp[i__] -= wk[i__] * yp[i__ + 1];
	npi = nn + i__;
	wk[i__] = wk[npi] - wk[i__] * wk[i__ + 1];
/* L2: */
    }

/* Solve the last equation for YP(N-1). */

    ypnm1 = (r1 + rnm1 - sdnm1 * yp[1] - sd1 * yp[nm2]) / (d1 + dnm1 + sdnm1 *
	     wk[1] + sd1 * wk[nm2]);

/* Back substitute for the remainder of the solution */
/*   components. */

    yp[nm1] = ypnm1;
    i__1 = nm2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yp[i__] += wk[i__] * ypnm1;
/* L3: */
    }

/* YP(N) = YP(1). */

    yp[*n] = yp[1];
    *ier = 0;
    return 0;

/* N is outside its valid range. */

L4:
    *ier = 1;
    return 0;

/* Abscissae out of order or duplicate points. */

L5:
    *ier = i__ + 1;
    return 0;
} /* ypc2p_ */

/* Subroutine */ int ypcoef_(doublereal *sigma, doublereal *dx, doublereal *
	d__, doublereal *sd)
{
    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static doublereal e, coshm, sinhm, ssinh, coshmm;
    extern /* Subroutine */ int snhcsh_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal scm, sig, ems, ssm;


/* *********************************************************** */

/*                                                From TSPACK */
/*                                            Robert J. Renka */
/*                                  Dept. of Computer Science */
/*                                       Univ. of North Texas */
/*                                           renka@cs.unt.edu */
/*                                                   11/17/96 */

/*   This subroutine computes the coefficients of the deriva- */
/* tives in the symmetric diagonally dominant tridiagonal */
/* system associated with the C-2 derivative estimation pro- */
/* cedure for a Hermite interpolatory tension spline. */

/* On input: */

/*       SIGMA = Nonnegative tension factor associated with */
/*               an interval. */

/*       DX = Positive interval width. */

/* Input parameters are not altered by this routine. */

/* On output: */

/*       D = Component of the diagonal term associated with */
/*           the interval.  D = SIG*(SIG*COSHM(SIG) - */
/*           SINHM(SIG))/(DX*E), where SIG = SIGMA and E = */
/*           SIG*SINH(SIG) - 2*COSHM(SIG). */

/*       SD = Subdiagonal (superdiagonal) term.  SD = SIG* */
/*            SINHM(SIG)/E. */

/* Module required by YPCOEF:  SNHCSH */

/* Intrinsic function called by YPCOEF:  EXP */

/* *********************************************************** */


    sig = *sigma;
    if (sig < 1e-9) {

/* SIG = 0:  cubic interpolant. */

	*d__ = 4. / *dx;
	*sd = 2. / *dx;
    } else if (sig <= .5) {

/* 0 .LT. SIG .LE. .5:  use approximations designed to avoid */
/*                      cancellation error in the hyperbolic */
/*                      functions when SIGMA is small. */

	snhcsh_(&sig, &sinhm, &coshm, &coshmm);
	e = (sig * sinhm - coshmm - coshmm) * *dx;
	*d__ = sig * (sig * coshm - sinhm) / e;
	*sd = sig * sinhm / e;
    } else {

/* SIG > .5:  scale SINHM and COSHM by 2*EXP(-SIG) in order */
/*            to avoid overflow when SIGMA is large. */

	ems = exp(-sig);
	ssinh = 1. - ems * ems;
	ssm = ssinh - sig * 2. * ems;
	scm = (1. - ems) * (1. - ems);
	e = (sig * ssinh - scm - scm) * *dx;
	*d__ = sig * (sig * scm - ssm) / e;
	*sd = sig * ssm / e;
    }
    return 0;
} /* ypcoef_ */
