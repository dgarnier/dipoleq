/*
**	stdio_dmatrix.c
**
**	A utility to read and write dmatrix created
**	by nrutil.c
**
**
** File:		stdio_dmatrix.c
** Date:		March 24, 1993
**
** Modifications:
**
**		October 12, 1993		Fixed bug with non-zero NR_END
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdio.h>
#include "stdio_dmatrix.h"

/*
**	The following should be the same as found in nrutil.c
*/
#define NR_END 			1

size_t        fread_dmatrix(double **dm, long nrl, long nrh, long ncl, long nch, FILE * fi)
{
	size_t        c;			/* read count */
	long          nrow = nrh - nrl + 1, ncol = nch - ncl + 1;

	c = fread(dm[nrl] - NR_END + ncl, sizeof(double), nrow * ncol + NR_END, fi);

	return c - NR_END;
}

size_t        fwrite_dmatrix(double **dm, long nrl, long nrh, long ncl, long nch, FILE * fi)
{
	size_t        c;			/* write count */
	long          nrow = nrh - nrl + 1, ncol = nch - ncl + 1;

	c = fwrite(dm[nrl] - NR_END + ncl, sizeof(double), nrow * ncol + NR_END, fi);

	return c - NR_END;
}
