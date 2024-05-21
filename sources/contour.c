/***
 * dcontour.c
 *
 * This is a neat routine for tracing contours.
 *
 * Makes a contourplot with z[x,y] with x[nx1..nx2] and y[ny1..ny2].
 *
 * This is a very clever program based on a routine by
 * B. R. Heap, Nat. Phys. Labs in Teddington, England.
 * The whole idea is to trace out the contours with a convention that HIGH
 * GROUND is always on your RIGHT.
 *
 * The procedure has three basic steps:
 *
 * STEP #1: Search through the array for those points representing
 *  the high side of a contour when moving in a right-handed
 *  manner. The array UNUSED holds the result of this.
 *
 * STEP #2: Trace out all OPEN CONTOURS. If you trace out one of
 *  the UNUSED locations, set it FALSE.
 *
 * STEP #3: Trace out all of the CLOSED CONTOURS.
 *
 *  This differs from ContourM in that it does not include the "midpoint"
 *  of a grid cell to help determine the path of the contour.	This one
 *  is faster and uses fewer points. 	 (Nicer if smoothed by MacDraw...)
 *
 * (c) M. E. Mauel -- Columbia University, New York City, March 16, 1992
 ***/

#include "nrutil.h"
#include "contour.h"

/* These six "cases" describe how a contour crosses a rectangular grid */
#define CASE_A 	1				/* Sharp right */
#define CASE_B 	2				/* Through top */
#define CASE_C 	3				/* Sharp left */
#define CASE_D	4
#define CASE_E	5
#define CASE_F	6
#define CASE_BAD 0

  typedef struct delsquare {
	  int           ixa;
	  int           iya;
	  int           ixd;
	  int           iyd;
	  int           ix;
	  int           iy;
  }
DELSQUARE;

/*** Shared variables within this file ***/

int         **unused;
int           ixmin, ixmax;
int           iymin, iymax;
double       *x, *y, **z;
double        height;

double        midpoint(double, double, double, double, double);
void          copydel(DELSQUARE *, DELSQUARE *);
void          Trace1(DELSQUARE *, DELSQUARE *, int, void (*)(double, double, double, int));
void          Trace2(DELSQUARE *, DELSQUARE *, int, void (*)(double, double, double, int));
void          trace_open_contour(DELSQUARE *, int, int, void (*)(double, double, double, int));
void          trace_closed_contour(DELSQUARE *, int, int, void (*)(double, double, double, int));
void          trace_unknown_contour(DELSQUARE *, int, int, void (*)(double, double, double, int));

/*** Functions and Subroutines ***/

double        midpoint(double r1, double z1, double r2, double z2, double h)
{
	double        ans;
	ans = r1 + (h - z1) * (r2 - r1) / (z2 - z1);
	return ans;
}

void          copydel(DELSQUARE * from, DELSQUARE * to)
{
	to->ix = from->ix;
	to->iy = from->iy;
	to->ixa = from->ixa;
	to->iya = from->iya;
	to->ixd = from->ixd;
	to->iyd = from->iyd;
}

/***
 * This procedure traces a contour across one square of the grid.
 * It presumes that the point has been set previously and we only need to
 * "LINETO" to the various paths of the contour.
 * The values (ixa,iya) represent the increment from (ix,iy) to the other
 * location of the ENTERING SIDE. (ixd,iyd) represent the increment to oposite
 * side. Notice that this can represent all four ways in which a contour can
 * enter a square:
 *
 * ENTER BOTTOM: ixa = -1 iya = 0  ixd = 0  iyd = 1
 * ENTER RIGHT:  ixa = 0  iya = -1 ixd = -1 iyd = 0
 * ENTER TOP:    ixa = 1  iya = 0  ixd = 0  iyd = -1
 * ENTER LEFT:   ixa = 0  iya = 1  ixd = 1  iyd = 0
 *
 * Once we've set up our aliases for our "standard square", we note that there
 * are six possible ways of tracing out the contour depending on the extent
 * of the high ground including the value of z(x,y) at the center of the square
 * which we represent by a simple average of the corners.
 *
 * D-------C
 * |\	  /|
 * | \ 	 / |
 * |  \ /  |
 * |   M   |
 * |  / \  |
 * | / 	 \ |
 * |/	  \|
 * A----^--B
 * 	    |
 * enter here!
 *
 ***/
void          Trace1(DELSQUARE * in, DELSQUARE * out, int ccode, void (*doAT) (double, double, double, int))
{
	double        za, zb, zc, zd;
	double        xa, xb, xc, xd;
	double        ya, yb, yc, yd;
	double        xx, yy;
	int           case_type = CASE_BAD;

	xa = x[in->ix + in->ixa];
	xb = x[in->ix];
	xc = x[in->ix + in->ixd];
	xd = x[in->ix + in->ixa + in->ixd];
	ya = y[in->iy + in->iya];
	yb = y[in->iy];
	yc = y[in->iy + in->iyd];
	yd = y[in->iy + in->iya + in->iyd];
	za = z[in->ix + in->ixa][in->iy + in->iya];
	zb = z[in->ix][in->iy];
	zc = z[in->ix + in->ixd][in->iy + in->iyd];
	zd = z[in->ix + in->ixa + in->ixd][in->iy + in->iya + in->iyd];

	/*In the following, we know that za < height and zb >= height.*/
	if (zc < height)
		case_type = CASE_A;
	if ((zd >= height) && (zc >= height))
		case_type = CASE_C;
	if ((zd < height) && (zc >= height))
		case_type = CASE_B;

	if (ccode == CONTOUR_DRAW) {
		switch (case_type) {
		  case CASE_A:			/* Sharp right */
			  xx = midpoint(xb, zb, xc, zc, height);
			  yy = midpoint(yb, zb, yc, zc, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  break;
		  case CASE_B:			/* Through top */
			  xx = midpoint(xc, zc, xd, zd, height);
			  yy = midpoint(yc, zc, yd, zd, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  break;
		  case CASE_C:			/* Sharp left */
			  xx = midpoint(xd, zd, xa, za, height);
			  yy = midpoint(yd, zd, ya, za, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  break;
                case CASE_BAD : nrerror("Program error in contour.");

		}
	}
	switch (case_type) {
	  case CASE_A:				/* Sharp right */
		  /* Rotate our orientation left.*/
		  out->ix = in->ix;
		  out->iy = in->iy;
		  out->ixa = in->iya;
		  out->iya = -in->ixa;
		  out->ixd = in->iyd;
		  out->iyd = -in->ixd;
		  break;
	  case CASE_B:				/* Through top */
		  /* Translate our orientation upwards*/
		  out->ix = in->ix + in->ixd;
		  out->iy = in->iy + in->iyd;
		  out->ixa = in->ixa;
		  out->iya = in->iya;
		  out->ixd = in->ixd;
		  out->iyd = in->iyd;
		  break;
	  case CASE_C:				/* Sharp left */
		  /* Rotate our orientation right*/
		  out->ix = in->ix + in->ixa + in->ixd;
		  out->iy = in->iy + in->iya + in->iyd;
		  out->ixa = -in->iya;
		  out->iya = in->ixa;
		  out->ixd = -in->iyd;
		  out->iyd = in->ixd;
		  break;
            case CASE_BAD : nrerror("Program error in contour.");
	}
}

/***
 * MidPoint Contouring
 *
 * This routine also takes the average of the four corner grid points, creating
 * a fifth "midpoint".  This adds more detail to the contours.
 *
 ***/

void          Trace2(DELSQUARE * in, DELSQUARE * out, int ccode, void (*doAT) (double, double, double, int))
{
	double        za, zb, zc, zd;
	double        xa, xb, xc, xd;
	double        ya, yb, yc, yd;
	double        xx, yy;
	double        xm, ym, zm;
	int           case_type = CASE_BAD;

	xa = x[in->ix + in->ixa];
	xb = x[in->ix];
	xc = x[in->ix + in->ixd];
	xd = x[in->ix + in->ixa + in->ixd];
	ya = y[in->iy + in->iya];
	yb = y[in->iy];
	yc = y[in->iy + in->iyd];
	yd = y[in->iy + in->iya + in->iyd];
	za = z[in->ix + in->ixa][in->iy + in->iya];
	zb = z[in->ix][in->iy];
	zc = z[in->ix + in->ixd][in->iy + in->iyd];
	zd = z[in->ix + in->ixa + in->ixd][in->iy + in->iya + in->iyd];

	xm = 0.5 * (xb + xd);		/* Look at diagonals and always find midpoint */
	ym = 0.5 * (yd + yb);
	zm = 0.25 * (za + zb + zc + zd);

	/* In the following, we know that za < height and zb >= height.*/
	if ((zc < height) && (zm < height))
		case_type = CASE_A;
	if ((zd >= height) && (zm >= height))
		case_type = CASE_D;
	if ((zd < height) && (zc < height) && (zm >= height))
		case_type = CASE_F;
	if ((zd >= height) && (zc >= height) && (zm < height))
		case_type = CASE_C;
	if ((zd < height) && (zc >= height) && (zm < height))
		case_type = CASE_B;
	if ((zd < height) && (zc >= height) && (zm >= height))
		case_type = CASE_E;

	if (ccode == CONTOUR_DRAW) {
		switch (case_type) {
		  case CASE_A:			/* Sharp right */
			  xx = midpoint(xb, zb, xm, zm, height);
			  yy = midpoint(yb, zb, ym, zm, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  xx = midpoint(xb, zb, xc, zc, height);
			  yy = midpoint(yb, zb, yc, zc, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  break;
		  case CASE_B:			/* Through top, right of M */
			  xx = midpoint(xb, zb, xm, zm, height);
			  yy = midpoint(yb, zb, ym, zm, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  xx = midpoint(xc, zc, xm, zm, height);
			  yy = midpoint(yc, zc, ym, zm, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  xx = midpoint(xc, zc, xd, zd, height);
			  yy = midpoint(yc, zc, yd, zd, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  break;
		  case CASE_C:			/* Through left, go around M */
			  xx = midpoint(xb, zb, xm, zm, height);
			  yy = midpoint(yb, zb, ym, zm, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  xx = midpoint(xc, zc, xm, zm, height);
			  yy = midpoint(yc, zc, ym, zm, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  xx = midpoint(xd, zd, xm, zm, height);
			  yy = midpoint(yd, zd, ym, zm, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  xx = midpoint(xd, zd, xa, za, height);
			  yy = midpoint(yd, zd, ya, za, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  break;
		  case CASE_D:			/* Sharp left */
			  xx = midpoint(xa, za, xm, zm, height);
			  yy = midpoint(ya, za, ym, zm, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  xx = midpoint(xa, za, xd, zd, height);
			  yy = midpoint(ya, za, yd, zd, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  break;
		  case CASE_E:			/* Through top, left of M */
			  xx = midpoint(xa, za, xm, zm, height);
			  yy = midpoint(ya, za, ym, zm, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  xx = midpoint(xd, zd, xm, zm, height);
			  yy = midpoint(yd, zd, ym, zm, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  xx = midpoint(xc, zc, xd, zd, height);
			  yy = midpoint(yc, zc, yd, zd, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  break;
		  case CASE_F:			/* Through right, go around M */
			  xx = midpoint(xa, za, xm, zm, height);
			  yy = midpoint(ya, za, ym, zm, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  xx = midpoint(xd, zd, xm, zm, height);
			  yy = midpoint(yd, zd, ym, zm, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  xx = midpoint(xc, zc, xm, zm, height);
			  yy = midpoint(yc, zc, ym, zm, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  xx = midpoint(xc, zc, xb, zb, height);
			  yy = midpoint(yc, zc, yb, zb, height);
			  (*doAT) (xx, yy, height, CONTOUR_TRACE);
			  break;
                    case CASE_BAD : nrerror("Program error in contour.");
		}
	}
	switch (case_type) {
	  case CASE_A:
	  case CASE_F:				/* Sharp right */
		  /* Rotate our orientation left.*/
		  out->ix = in->ix;
		  out->iy = in->iy;
		  out->ixa = in->iya;
		  out->iya = -in->ixa;
		  out->ixd = in->iyd;
		  out->iyd = -in->ixd;
		  break;
	  case CASE_B:
	  case CASE_E:				/* Through top */
		  /* Translate our orientation upwards*/
		  out->ix = in->ix + in->ixd;
		  out->iy = in->iy + in->iyd;
		  out->ixa = in->ixa;
		  out->iya = in->iya;
		  out->ixd = in->ixd;
		  out->iyd = in->iyd;
		  break;
	  case CASE_C:
	  case CASE_D:				/* Sharp left */
		  /* Rotate our orientation right*/
		  out->ix = in->ix + in->ixa + in->ixd;
		  out->iy = in->iy + in->iya + in->iyd;
		  out->ixa = -in->iya;
		  out->iya = in->ixa;
		  out->ixd = -in->iyd;
		  out->iyd = in->ixd;
		  break;
            case CASE_BAD : nrerror("Program error in contour.");

	}
}

void          trace_open_contour(DELSQUARE * in1, int ccode, int midpt, void (*doAT) (double, double, double, int))
{
	int           edge = 0;
	DELSQUARE     in, out;
	double        xx, yy;
	int           cdraw;

	cdraw = (((ccode == CONTOUR_ALL) || (ccode == CONTOUR_ONLY_OPEN)) ?
			 CONTOUR_DRAW : CONTOUR_DONT_DRAW);

	copydel(in1, &in);
	xx = midpoint(x[in.ix], z[in.ix][in.iy], x[in.ix + in.ixa], z[in.ix + in.ixa][in.iy + in.iya], height);
	yy = midpoint(y[in.iy], z[in.ix][in.iy], y[in.iy + in.iya], z[in.ix + in.ixa][in.iy + in.iya], height);
	if (cdraw == CONTOUR_DRAW)
		(*doAT) (xx, yy, height, CONTOUR_START);

	while (!edge) {
		if (midpt == CONTOUR_NO_MIDPOINT)
			Trace1(&in, &out, cdraw, doAT);
		else
			Trace2(&in, &out, cdraw, doAT);
		copydel(&out, &in);
		if ((in.ix == ixmin) && (in.ixd == -1))	/*If entering right*/
			edge = 1;
		if ((in.iy == iymin) && (in.iyd == -1))	/*If entering top*/
			edge = 1;
		if ((in.ix == ixmax) && (in.ixd == 1))	/*If entering left*/
			edge = 1;
		if ((in.iy == iymax) && (in.iyd == 1))	/*If entering bottom*/
			edge = 1;
		unused[in.ix][in.iy] = 0;
	}
	if (cdraw == CONTOUR_DRAW)
		(*doAT) (xx, yy, height, CONTOUR_STOP);	/* dummy arguments during STOP */
}

void          trace_closed_contour(DELSQUARE * in1, int dummy, int midpt, void (*doAT) (double, double, double, int))
{
#pragma unused(dummy)
	DELSQUARE     in, out;
	int           finished = 0;
	double        xx, yy;

	/* Each open contour starts from the bottom by construction of UNUSED. */
	copydel(in1, &in);

	xx = midpoint(x[in.ix], z[in.ix][in.iy], x[in.ix + in.ixa], z[in.ix + in.ixa][in.iy + in.iya], height);
	yy = midpoint(y[in.iy], z[in.ix][in.iy], y[in.iy + in.iya], z[in.ix + in.ixa][in.iy + in.iya], height);
	(*doAT) (xx, yy, height, CONTOUR_START);

	while (!finished) {
		if (midpt == CONTOUR_NO_MIDPOINT)
			Trace1(&in, &out, CONTOUR_DRAW, doAT);
		else
			Trace2(&in, &out, CONTOUR_DRAW, doAT);
		copydel(&out, &in);
		if ((in.ix == in1->ix) && (in.iy == in1->iy) && (in.ixa == -1) && (in.iya == 0))
			finished = 1;
		unused[in.ix][in.iy] = 0;
	}
	(*doAT) (xx, yy, height, CONTOUR_STOP);	/* dummy arguments during STOP */
}

void  contour_point(double *xarg, double *yarg, double **zarg,
	  int nx1, int nx2, int ny1, int ny2, double xx, double yy, int midpt,
					  void          (*doAT) (double, double, double, int))
/* traces a contour (either open or closed) from a starting point xx, yy */
/* currently this has a bug... it will only do half of an open contour */
{
	DELSQUARE   in, out, in1;
	int 		done = 0;
	int         ix, iy;
	double      zz,t,u;

	/* Copy to file-global variables */
	x = xarg;
	y = yarg;
	z = zarg;
	ixmin = nx1;
	ixmax = nx2;
	iymin = ny1;
	iymax = ny2;

    /* find border grid points */
	for (ix = ixmin; ix < ixmax; ix++)
		if (((x[ix] <= xx) && (xx < x[ix+1])) ||
		    ((x[ix] >= xx) && (xx > x[ix+1]))    ) break;
	if (ix == ixmax) nrerror("Error in contour_point:  bad xx input\n");

	for (iy = iymin; iy < iymax; iy++)
		if (((y[iy] <= yy) && (yy < y[iy+1])) ||
		    ((y[iy] >= yy) && (yy > y[iy+1]))    ) break;
	if (iy == iymax) nrerror("Error in contour_point:  bad yy input\n");

    /* find z[xx,yy] by bilinear interpolation */

	t = (xx - x[ix])/(x[ix+1]-x[ix]);
	u = (yy - y[iy])/(y[iy+1]-y[iy]);

	zz = height = (1.0 - t) * (1.0 - u) * z[ix][iy] +
	                 t * (1.0 - u) * z[ix + 1][iy] +
			         t * u         * z[ix + 1][iy + 1] +
			 (1.0 - t) * u         * z[ix][iy + 1];

	/* now lets find the orientation of the grid */
	if        ( (z[ix][iy] <= zz) && (z[ix+1][iy] > zz) ) {     /* enter bottom */
		in.ix = ix+1;
		in.iy = iy;
		in.ixa = -1;
		in.iya = 0;
		in.ixd = 0;
		in.iyd = 1;

	} else if ( (z[ix+1][iy] <= zz) && (z[ix+1][iy+1] > zz) ) { /* enter right */
		in.ix = ix+1;
		in.iy = iy+1;
		in.ixa = 0;
		in.iya = -1;
		in.ixd = -1;
		in.iyd = 0;

	} else if ( (z[ix+1][iy+1] <= zz) && (z[ix][iy+1] > zz) ) { /* enter top */
		in.ix = ix;
		in.iy = iy+1;
		in.ixa = 1;
		in.iya = 0;
		in.ixd = 0;
		in.iyd = -1;
	} else if ( (z[ix][iy+1] <= zz) && (z[ix][iy] > zz) ) {     /* enter left */
		in.ix = ix;
		in.iy = iy;
		in.ixa = 0;
		in.iya = 1;
		in.ixd = 1;
		in.iyd = 0;
	} else
		nrerror("Error in contour_point: Internal error . \n");

	copydel(&in,&in1);

	/* tell the drawing routine to start! */
	(*doAT) (xx, yy, height, CONTOUR_START);

	/* now lets find the appropriat */

	while (done == 0) {
		if (midpt == CONTOUR_NO_MIDPOINT)
			Trace1(&in, &out, CONTOUR_DRAW, doAT);
		else
			Trace2(&in, &out, CONTOUR_DRAW, doAT);
		copydel(&out, &in);
		if ((in.ix == in1.ix) && (in.iy == in1.iy) && (in.ixa == in1.ixa) && (in.iya == in1.iya))
			done = 1;
		if ((in.ix == ixmin) && (in.ixd == -1))	/*If entering right*/
			done = 2;
		if ((in.iy == iymin) && (in.iyd == -1))	/*If entering top*/
			done = 2;
		if ((in.ix == ixmax) && (in.ixd == 1))	/*If entering left*/
			done = 2;
		if ((in.iy == iymax) && (in.iyd == 1))	/*If entering bottom*/
			done = 2;
	}

	(*doAT) (xx, yy, height, CONTOUR_STOP);
}

void          contour(double *xarg, double *yarg, double **zarg,
	  int nx1, int nx2, int ny1, int ny2, double h, int ccode, int midpt,
					  void          (*doAT) (double, double, double, int))
{
	DELSQUARE     in;
	int           ix, iy;

	/* Copy to file-global variables */
	x = xarg;
	y = yarg;
	z = zarg;
	ixmin = nx1;
	ixmax = nx2;
	iymin = ny1;
	iymax = ny2;
	height = h;

	/*** STEP #1 ***/
	/* First, look at each cell and decide whether or not a contour passes with */
	/* high ground toward increasing ix. 			*/
	unused = imatrix(ixmin, ixmax, iymin, iymax);
	for (ix = ixmin + 1; ix <= ixmax; ix++)
		for (iy = iymin + 1; iy < iymax; iy++)
			unused[ix][iy] = (z[ix - 1][iy] < height) && (z[ix][iy] >= height);

	/*** STEP #2 ***/
	/* Now, we trace all open contours...*/

	/*** ENTER Bottom ***/
	for (ix = ixmin + 1; ix <= ixmax; ix++) {
		in.ix = ix;
		in.iy = iymin;
		in.ixa = -1;
		in.iya = 0;
		in.ixd = 0;
		in.iyd = 1;
		if ((z[ix - 1][iymin] < height) && (z[ix][iymin] >= height))
			trace_open_contour(&in, ccode, midpt, doAT);
	}
	/*** ENTER Right ***/
	for (iy = iymin + 1; iy <= iymax; iy++) {
		in.ix = ixmax;
		in.iy = iy;
		in.ixa = 0;
		in.iya = -1;
		in.ixd = -1;
		in.iyd = 0;
		if ((z[ixmax][iy - 1] < height) && (z[ixmax][iy] >= height))
			trace_open_contour(&in, ccode, midpt, doAT);
	}
	/*** ENTER Top ***/
	for (ix = ixmax - 1; ix >= ixmin; ix--) {
		in.ix = ix;
		in.iy = iymax;
		in.ixa = 1;
		in.iya = 0;
		in.ixd = 0;
		in.iyd = -1;
		if ((z[ix + 1][iymax] < height) && (z[ix][iymax] >= height))
			trace_open_contour(&in, ccode, midpt, doAT);
	}
	/*** ENTER Left ***/
	for (iy = iymax - 1; iy >= iymin; iy--) {
		in.ix = ixmin;
		in.iy = iy;
		in.ixa = 0;
		in.iya = 1;
		in.ixd = 1;
		in.iyd = 0;
		if ((z[ixmin][iy + 1] < height) && (z[ixmin][iy] >= height))
			trace_open_contour(&in, ccode, midpt, doAT);
	}

	if ((ccode == CONTOUR_ALL) || (ccode == CONTOUR_ONLY_CLOSED)) {
		/*** STEP #3 ***/
		/* Finally, we loop within the interior for closed contours.  */

		for (ix = ixmin + 1; ix <= ixmax; ix++)
			for (iy = iymin + 1; iy < iymax; iy++)
				if (unused[ix][iy] == 1) {
					in.ix = ix;
					in.iy = iy;
					in.ixa = -1;
					in.iya = 0;
					in.ixd = 0;
					in.iyd = 1;
					trace_closed_contour(&in, ccode, midpt, doAT);
				}
	}
	free_imatrix(unused, ixmin, ixmax, iymin, iymax);
}
