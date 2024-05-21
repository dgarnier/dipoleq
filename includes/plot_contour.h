/*
**
**    plot_contour.h
**
**    12/12/98
**
**    make a contour plot using the ClibPDF routines and
**    contour.c from the tokamac code.
**
**    (c) D. Garnier -- Columbia University
*/

#ifndef _PLOT_CONTOUR_H_
#define _PLOT_CONTOUR_H_ 1

#ifdef __cplusplus
extern "C" {
#endif

#define NOFILL  0
#define FILL    1
#define NOMIDPT 0
#define MIDPT   2

void  plot_contour(CPDFdoc *pdf, int fill, double *xarg, double *yarg, double **zarg, double level,
	  int nx1, int nx2, int ny1, int ny2);

#ifdef __cplusplus
}
#endif



#endif /* _PLOT_CONTOUR_H_ */
