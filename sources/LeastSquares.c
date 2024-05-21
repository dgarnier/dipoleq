/*
**	TokaMac v2.0
**
**	LeastSquares.c
**
**	This file  contains the least-squares fitting routines which
**	adjust the plasma parameters in order to fit the
**	experimental measurements.
**
**
**	File:		LeastSquares.c
**	Date:		March 20, 1993
**
**	Revisions:
**
**		August 6, 1993		Added meas->Now
**		Oct. 3, 1993		Added covariance calculation
**		Nov. 16, 1993		Added singluar values to SVD output
**
**
**	(c) L. Bai and M. Mauel -- Columbia University
*/

#include <math.h>
#include <stdlib.h>
#include "VAX_Alloc.h"
#include <stdio.h>
#include <string.h>
#include "nrutil.h"
#include "SVDFit.h"
#include "plasma.h"
#include "psigrid.h"
#include "coil.h"
#include "measurement.h"
#include "tokamak.h"
#include "Find_dJdy.h"
#include "FindMeasFit.h"
#include "LeastSquares.h"

extern double ***meas_dJdy;		/* the pointer used for magnetic L array */

extern FILE  *LogFile;


void          SVDFitOutput(TOKAMAK * td, double **L, double *w, double *uold, double *unew,
						   double chisq1, double chisq2);


/*
**	SVDFitOutput
**
**	This writes to a file the L array and the SVDFit info.
*/
void          SVDFitOutput(TOKAMAK * td, double **L, double *w, double *uold, double *unew,
						   double chisq1, double chisq2)
{
	FILE         *fi;
	int           im, iu;
	char          fname[32] = "";
	double        **Covar;

	strncat(fname, td->Oname, 20);	/* take 1st 20 characters */
	strcat(fname, "_SVDFit.out");

	fi = fopen(fname, "w");
	if (!fi)
		nrerror("ERROR:	Could not open file for writing in SVDFitOutput.");

	Covar = td->Covar;

	/*  H E A D E R */
	fprintf(fi, "TokaMac SVDFit Output. From Input FileName: %s\n", td->Iname);
	fprintf(fi, "    RunName: %s. Info: %s\n", td->Name, td->Info);
	fprintf(fi, "    Run started at %s\n", td->Start);
	fprintf(fi, "    Chisq1 = %g, Chisq2 = %g\n", chisq1, chisq2);
	fprintf(fi, "\n");

	/* U N K N O W N S */
	/*          1234567890123456789012345678901234567890123456789012345678901234567890 */
	fprintf(fi, "Unknowns     uOld       uNew      +-err  +-PerCent   sing val\n");
	for (iu = 1; iu <= td->NumUnkns; iu++)
		fprintf(fi, " %3d   %10.3g %10.3g %10.3g %10.3g %10.3g\n",
				iu, uold[iu], unew[iu], sqrt(Covar[iu][iu]),
				100.0 * sqrt(Covar[iu][iu]) / fabs(unew[iu]), w[iu]);
	fprintf(fi, "\n");

	/* C O V A R I A N C E */
	/*          1234567890123456789012345678901234567890123456789012345678901234567890 */
	fprintf(fi, "Covariance Martix\n");
	fprintf(fi, "Unkn |  Unknowns ->   \n");
	fprintf(fi, "     ");
	for (iu = 1; iu <= td->NumUnkns; iu++)
		fprintf(fi, "        %2d ", iu);
	fprintf(fi, "\n");
	for (im = 1; im <= td->NumUnkns; im++) {
		fprintf(fi, " %3d ", im);
		for (iu = 1; iu <= td->NumUnkns; iu++)
			fprintf(fi, "%10.2g ", Covar[im][iu]);
		fprintf(fi, "\n");
	}
	fprintf(fi, "\n");

	/* M E A S U R E M E N T S */
	/*          1234567890123456789012345678901234567890123456789012345678901234567890 */
	fprintf(fi, "L Array\n");
	fprintf(fi, "Meas |  Unknowns ->   \n");
	fprintf(fi, "     ");
	for (iu = 1; iu <= td->NumUnkns; iu++)
		fprintf(fi, "        %2d ", iu);
	fprintf(fi, "\n");
	for (im = 1; im <= td->NumMeasures; im++) {
		fprintf(fi, " %3d ", im);
		for (iu = 1; iu <= td->NumUnkns; iu++)
			fprintf(fi, "%10.2g ", L[im][iu]);
		fprintf(fi, "\n");
	}
	fprintf(fi, "\n");

	fclose(fi);
}

/*
** CopyUnknowns
**
*/
void	CopyUnknowns(TOKAMAK *td, double *UnknOld)
{
	PLASMA	*pl;
	int		idx = 0,i;

	pl = td->Plasma;

	/* F I L L   P R E S E N T   U N K N O W N S */
	switch (pl->ModelType) {
	  case Plasma_Std:
		  UnknOld[1] = pl->G2p[1];
		  UnknOld[2] = pl->Pp[1];
		  idx = 2;
		  break;
	  case Plasma_IsoNoFlow:
		  for (i = 1; i <= pl->G2pTerms; i++)
			  UnknOld[i + idx] = pl->G2p[i];
		  idx = pl->G2pTerms;
		  for (i = 1; i <= pl->PpTerms; i++)
			  UnknOld[i + idx] = pl->Pp[i];
		  idx += pl->PpTerms;
		  break;
	  case Plasma_IsoFlow:
		  for (i = 1; i <= pl->G2pTerms; i++)
			  UnknOld[i + idx] = pl->G2p[i];
		  idx = pl->G2pTerms;
		  for (i = 1; i <= pl->HTerms; i++)
			  UnknOld[i + idx] = pl->H[i];
		  idx += pl->HTerms;
		  for (i = 1; i <= pl->RotTerms; i++)
			  UnknOld[i + idx] = pl->Rot[i];
		  idx += pl->RotTerms;
		  for (i = 1; i <= pl->SisoTerms; i++)
			  UnknOld[i + idx] = pl->Siso[i];
		  idx += pl->SisoTerms;
		  break;
	  case Plasma_AnisoNoFlow:
		  for (i = 1; i <= pl->G2pTerms; i++)
			  UnknOld[i + idx] = pl->G2p[i];
		  idx = pl->G2pTerms;
		  for (i = 1; i <= pl->HTerms; i++)
			  UnknOld[i + idx] = pl->H[i];
		  idx += pl->HTerms;
		  for (i = 1; i <= pl->SparTerms; i++)
			  UnknOld[i + idx] = pl->Spar[i];
		  idx += pl->SparTerms;
		  for (i = 1; i <= pl->SperTerms; i++)
			  UnknOld[i + idx] = pl->Sper[i];
		  idx += pl->SperTerms;
		  break;
	  case Plasma_AnisoFlow:
		  for (i = 1; i <= pl->G2pTerms; i++)
			  UnknOld[i + idx] = pl->G2p[i];
		  idx = pl->G2pTerms;
		  for (i = 1; i <= pl->HTerms; i++)
			  UnknOld[i + idx] = pl->H[i];
		  idx += pl->HTerms;
		  for (i = 1; i <= pl->RotTerms; i++)
			  UnknOld[i + idx] = pl->Rot[i];
		  idx += pl->RotTerms;
		  for (i = 1; i <= pl->SparTerms; i++)
			  UnknOld[i + idx] = pl->Spar[i];
		  idx += pl->SparTerms;
		  for (i = 1; i <= pl->SperTerms; i++)
			  UnknOld[i + idx] = pl->Sper[i];
		  idx += pl->SperTerms;
		  break;
	}
	for (i = 0; i < td->NumCoils; i++)
		UnknOld[i + 1 + idx] = td->Coils[i]->CoilCurrent;
}

/*
**	RewriteUnknowns
**
*/
void		RewriteUnknowns(TOKAMAK *td, double *UnknNew)
{
	PLASMA	*pl;
	int		idx = 0,i;

	pl = td->Plasma;

	/* U P D A T E    P L A S M A    V A L U E S */
	switch (pl->ModelType) {
	  case Plasma_Std:
		  pl->G2p[1] = UnknNew[1];
		  pl->Pp[1] = UnknNew[2];
		  idx = 2;
		  break;
	  case Plasma_IsoNoFlow:
		  for (i = 1; i <= pl->G2pTerms; i++)
			  pl->G2p[i] = UnknNew[i + idx];
		  idx = pl->G2pTerms;
		  for (i = 1; i <= pl->PpTerms; i++)
			  pl->Pp[i] = UnknNew[i + idx];
		  idx += pl->PpTerms;
		  break;
	  case Plasma_IsoFlow:
		  for (i = 1; i <= pl->G2pTerms; i++)
			  pl->G2p[i] = UnknNew[i + idx];
		  idx = pl->G2pTerms;
		  for (i = 1; i <= pl->HTerms; i++)
			  pl->H[i] = UnknNew[i + idx];
		  idx += pl->HTerms;
		  for (i = 1; i <= pl->RotTerms; i++)
			  pl->Rot[i] = UnknNew[i + idx];
		  idx += pl->RotTerms;
		  for (i = 1; i <= pl->SisoTerms; i++)
			  pl->Siso[i] = UnknNew[i + idx];
		  idx += pl->SisoTerms;
		  break;
	  case Plasma_AnisoNoFlow:
		  for (i = 1; i <= pl->G2pTerms; i++)
			  pl->G2p[i] = UnknNew[i + idx];
		  idx = pl->G2pTerms;
		  for (i = 1; i <= pl->HTerms; i++)
			  pl->H[i] = UnknNew[i + idx];
		  idx += pl->HTerms;
		  for (i = 1; i <= pl->SparTerms; i++)
			  pl->Spar[i] = UnknNew[i + idx];
		  idx += pl->SparTerms;
		  for (i = 1; i <= pl->SperTerms; i++)
			  pl->Sper[i] = UnknNew[i + idx];
		  idx += pl->SperTerms;
		  break;
	  case Plasma_AnisoFlow:
		  for (i = 1; i <= pl->G2pTerms; i++)
			  pl->G2p[i] = UnknNew[i + idx];
		  idx = pl->G2pTerms;
		  for (i = 1; i <= pl->HTerms; i++)
			  pl->H[i] = UnknNew[i + idx];
		  idx += pl->HTerms;
		  for (i = 1; i <= pl->RotTerms; i++)
			  pl->Rot[i] = UnknNew[i + idx];
		  idx += pl->RotTerms;
		  for (i = 1; i <= pl->SparTerms; i++)
			  pl->Spar[i] = UnknNew[i + idx];
		  idx += pl->SparTerms;
		  for (i = 1; i <= pl->SperTerms; i++)
			  pl->Sper[i] = UnknNew[i + idx];
		  idx += pl->SperTerms;
		  break;
	}
	for (i = 0; i < td->NumCoils; i++)
		td->Coils[i]->CoilCurrent = UnknNew[i + 1 + idx];
}

/*
**	LeastSquares
**
**
**	We use the Numerical Recipies routines for SVD.
**	These dimension the arrays from [1..NumMeasures][1..NumUnkns].
*/

void          LeastSquares(TOKAMAK * td, int IsFirstTime)
{
	PSIGRID      *pg;
	PLASMA       *pl;
	MEAS         *aMeas;
	int           im, iu, i;
	double       *m, *mNow;		/* measurements */
	double       *s;			/* StdDev */
	double       *UnknOld, *UnknNew;
	double       *w;
	double      **L, **u, **V;
	double        chisq1, chisq2 = 0.0, dot;

	pg = td->PsiGrid;
	pl = td->Plasma;

	/* H O W   M A N Y   U N K N O W N S  ? */
	switch (pl->ModelType) {
	  case Plasma_Std:
		  td->NumUnkns = 2 + td->NumCoils;
		  break;
	  case Plasma_IsoNoFlow:
		  td->NumUnkns = pl->PpTerms + pl->G2pTerms + td->NumCoils;
		  break;
	  case Plasma_IsoFlow:
		  td->NumUnkns = pl->HTerms + pl->G2pTerms +
			  pl->SisoTerms + pl->RotTerms + td->NumCoils;
		  break;
	  case Plasma_AnisoNoFlow:
		  td->NumUnkns = pl->HTerms + pl->G2pTerms +
			  pl->SperTerms + pl->SparTerms + td->NumCoils;
		  break;
	  case Plasma_AnisoFlow:
		  td->NumUnkns = pl->HTerms + pl->G2pTerms +
			  pl->SperTerms + pl->SparTerms +
			  pl->RotTerms + td->NumCoils;
		  break;
	}

	printf("INFO:	LeastSquares with %d measurements and %d unknowns.\n",
		   td->NumMeasures, td->NumUnkns);
	fprintf(LogFile, "INFO:	LeastSquares with %d measurements and %d unknowns.\n",
			td->NumMeasures, td->NumUnkns);

	if (td->NumMeasures < td->NumUnkns)
		nrerror("ERROR:	Too few measurements for these unknowns.");

	/* F I N D   M E A S U R E M E N T S   W I T H   P R E S E N T   U N K N O W N S */
	FindMeasNow(td);

	/* F I L L   L   A R R A Y */
	L = dmatrix(1, td->NumMeasures, 1, td->NumUnkns);
	meas_dJdy = new_dJdy(td->NumUnkns, pg->Nsize);
	Find_dJdy(td, meas_dJdy);
	for (im = 0; im < td->NumMeasures; im++) {
		aMeas = td->Measures[im];
		(*(aMeas->FindL)) (td, aMeas, L[im + 1]);
	}
	free_dJdy(meas_dJdy, td->NumUnkns, pg->Nsize);

	/*  A L L O C A T E   O T H E R   A R R A Y S */
	u = dmatrix(1, td->NumMeasures, 1, td->NumUnkns);
	V = dmatrix(1, td->NumUnkns, 1, td->NumUnkns);
	w = dvector(1, td->NumUnkns);
	UnknOld = dvector(1, td->NumUnkns);
	UnknNew = dvector0(1, td->NumUnkns);
	m = dvector(1, td->NumMeasures);
	mNow = dvector(1, td->NumMeasures);
	s = dvector(1, td->NumMeasures);
	if (IsFirstTime) {
		td->Covar = dmatrix(1, td->NumUnkns, 1, td->NumUnkns);
		td->UnknVectors = dmatrix(1, td->NumUnkns, 1, td->NumUnkns);
		td->SValues = dvector(1, td->NumUnkns);
	}

	/* F I L L   M   &   S   A R R A Y S */
	for (im = 0; im < td->NumMeasures; im++) {
		aMeas = td->Measures[im];
		m[im + 1] = aMeas->Value;
		mNow[im + 1] = aMeas->Now;
		s[im + 1] = aMeas->StdDev;
	}

	/* F I L L   P R E S E N T   U N K N O W N S */
	CopyUnknowns(td,UnknOld);

	/* Singular Value Decomposition based on Chapter 14 of NUMERICAL RECIPES. */

	/* N O R M A L I Z E   B Y   S T D D E V
	**
	** For a linear model (dot - mNow) == 0.0
	**
	** When IsFirstTime, UnknOld is usually zero.
	*/
	for (im = 1; im <= td->NumMeasures; im++) {
		dot = 0.0;
		for (iu = 1; iu <= td->NumUnkns; iu++) {
			dot += L[im][iu] * UnknOld[iu];
			L[im][iu] /= s[im];
		}
		dot -= mNow[im];
		m[im] = m[im] + dot;
		m[im] /= s[im];
	}

	/* S V D F I T */
	svdfit(L, u, w, V, UnknNew, m, td->NumMeasures, td->NumUnkns, &chisq1);

	/* C O V A R I A N C E   M A T R I X */
	svdvar(V, td->NumUnkns, w, td->Covar);

	/* C o p y   U n k n V e c t o r s */
	for (i=1; i<= td->NumUnkns; i++) {
		td->SValues[i] = w[i];
		for (iu = 1; iu <= td->NumUnkns; iu++)
			td->UnknVectors[i][iu] = V[i][iu];
	}

	/***
	* It is important to realize that considerable information is contained
	* in the SVD results.
	* Probably the most important of these are the principle axis of the
	* chisq confidence ellipses in the space of the unknowns.
	***/

	/* U N D E R  R E L A X   U N K N O W N S */
	if (!IsFirstTime)
		for (iu=1; i<= td->NumUnkns; iu++)
			UnknNew[iu] = (1.0 - pg->UnderRelax2)*UnknNew[iu] + pg->UnderRelax2*UnknOld[iu];

	/* U P D A T E    P L A S M A    V A L U E S */
	RewriteUnknowns(td, UnknNew);

	/* F I N D   M E A S   F I T */
	FindMeasFit(td);
	for (im = 0; im < td->NumMeasures; im++) {
		aMeas = td->Measures[im];
		chisq2 += DSQR((aMeas->Value - aMeas->Fit) / aMeas->StdDev);
	}

	pl->ChiSqr = chisq2;		/* this is more accurate than chisq1 */

	printf("		[ChiSq1 = %g, Chisq2 = %g]\n", chisq1, chisq2);
	fprintf(LogFile, "		[ChiSq1 = %g, Chisq2 = %g]\n", chisq1, chisq2);

	SVDFitOutput(td, L, w, UnknOld, UnknNew, chisq1, chisq2);

	free_dmatrix(L, 1, td->NumMeasures, 1, td->NumUnkns);

	free_dmatrix(u, 1, td->NumMeasures, 1, td->NumUnkns);
	free_dmatrix(V, 1, td->NumUnkns, 1, td->NumUnkns);
	free_dvector(w, 1, td->NumUnkns);
	free_dvector(UnknOld, 1, td->NumUnkns);
	free_dvector(UnknNew, 1, td->NumUnkns);
	free_dvector(m, 1, td->NumMeasures);
	free_dvector(mNow, 1, td->NumMeasures);
	free_dvector(s, 1, td->NumMeasures);

}
