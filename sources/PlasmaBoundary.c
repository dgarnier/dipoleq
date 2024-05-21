/*
** TokaMac v2.0
**
** PlasmaBoundary.c
**
**
**
** File:		PlasmaBoundary.c
** Date:		March 22, 1993
**
** Modification history:
**
**	April 7, 1993		Changed FindAllFlats to include assiging one
**						magnetic axis to the minimum of Psi.
**
**	Oct. 25, 1993		Added SetEdgeCurrent.
**
** August 6, 1996 Fixed typo in SetEdgeCurrent.
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "psigrid.h"
#include "limiter.h"
#include "plasma.h"
#include "separatrix.h"
#include "tokamak.h"
#include "PlasmaBoundary.h"
#include "multitask.h"

#define FALSE 0
#define TRUE 1
#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

#define SQUARE(x)	((x) * (x))
#define MAX(A,B) 	((A)>(B) ? (A) : (B))
#define MIN(A,B) 	((A)<(B) ? (A) : (B))

/*
**	The following definitions are used in FindAllFlats
*/
#define	GET_DX(ix,iz)	tdx*(Psi[(ix) + 1][iz] - Psi[(ix) - 1][iz])
#define	GET_DZ(ix,iz)	tdz*(Psi[ix][(iz) + 1] - Psi[ix][(iz) - 1])
#define GET_DX_4(ix,iz)  (- Psi[(ix)+2][iz] + Psi[(ix)-2][iz] + 8.0*(Psi[(ix)+1][iz] - Psi[(ix)-1][iz]))/12.0/dx
#define GET_DZ_4(ix,iz)  (- Psi[ix][(iz)+2] + Psi[ix][(iz)-2] + 8.0*(Psi[ix][(iz)+1] - Psi[ix][(iz)-1]))/12.0/dz

/*
**	The following definitions are used to make sure we do not
**	call pow(0.0, 0) when we compute the polynomial expansions
**	of the flux functions.
*/
#define	SMALL_POSITIVE_FACTOR	1.000001
#define	SMALL_NEGATIVE_FACTOR	0.999999
#define	SMALL_SMALL_FACTOR		-0.000001

  enum sides {
	  TOP, BOT, IN, OUT
  };

#define FLAT_NONE		0
#define FLAT_SEP		1
#define	FLAT_AXIS		2
#define MAX_FLATNUM		14		/* Maximum number of FLATS to locate */

  typedef struct flat {
	  int           fType;		/* FLAT_SEP or FLAT_AXIS */
	  double        X, Z;		/* Location where |grad(Psi)|^2 = 0 */
	  double        Psi;		/* Value of Psi at (X,Z) */
  }
FLAT;
#ifndef IDL
extern FILE  *LogFile;
#endif
/*
**
**	Private function prototypes
**
*/

void          FindAllFlats(PSIGRID *, FLAT *);
void          FindValidFlats(TOKAMAK * td, FLAT * fl);
void          MarkDivertors(TOKAMAK *, int **);
int           CheckIsDivertor(PSIGRID *, int **, double, double);
void          FindMinPsiLimiter(PSIGRID *, LIMITER *, int **);
void          FindMaxPsiLimiter(PSIGRID *, LIMITER *, int **);
double        get_dxdz(double **Psi, int ix, int iz);
int           IsTruePlasma(PSIGRID * pg, int ixp, int izp, int ixa, int iza);
void          FinalCheckIsPlasma(PSIGRID * pg);
void          SetEdgeCurrent(TOKAMAK *td);


double        get_dxdz(double **Psi, int ix, int iz)
{
	double        sum = 0.0;
	int           ixp, izp, dx, dz;

	for (dx = -1; dx <= 1; dx = dx + 2)
		for (dz = -1; dz <= 1; dz = dz + 2) {
			ixp = ix + dx;
			izp = iz + dz;
			sum += dx * dz * (64.0 * Psi[ixp][izp] -
						8.0 * (Psi[ixp + dx][izp] + Psi[ixp][izp + dz]) +
							  Psi[ixp + dx][izp + dz]);
		}

	return sum / 144.0;
}

/*
**
**	FindAllFlats
**
**	This function creates an array filled with information
**	concerning the location of all separatrices and all
**	magnetic axes.
**
**	A separatrix or magnetic axis is defined when
**		|�p/�x|^2 + |�p/�z|^2 < |�p/�x|^2 + |�p/�z|^2 on all sides
**
**	The minimum is a magnetic axis when
**		�2p/�x2 * �2p/�z2 - (�2p/�z�x)^2 > 0
**
**	The minimum is a separatrix when
**		�2p/�x2 * �2p/�z2 - (�2p/�z�x)^2 < 0
**
**	Interpolation is used to find the actual coordinates for the
**	where |grad(Psi)|^2 vanishes.
**
**	INPUTS:
**		PSIGRID		*pg		==> a pointer to the PSIGRID data
**		FLAT		*fl		==> array of FLATs to hold results
**
**
**	April 7, 1993	We also add a loop to find the minimum of
**					Psi as one of the estimates for the magnetic
**					axis.
**
**	April 7, 1993	We increased the order of finite-difference
**					formula used for �2p/�x2, �2p/�z2, and �2p/�z�x.
**					We now use 5 point formula.
**
*/
void          FindAllFlats(PSIGRID * pg, FLAT * fl)
{
	int           ix, iz, i, s, nmax;
	double      **Psi;
	double        dPsiX4[4], dPsiZ4[4];
	double        dPsiZ, dPsiX, dPsiXX, dPsiZZ, dPsiXZ, grad, det;
	int           IsMagAxis, IsSeparatrix, IsMinimum;
	int           FoundAxis = FALSE;
	double        xFlat, zFlat, PsiFlat;
	double        dx, dz, tdx, tdz;

	Psi = pg->Psi;
	nmax = pg->Nsize;
	dx = pg->dx;
	dz = pg->dz;
	tdx = 0.5 / dx;
	tdz = 0.5 / dz;
MULTI;


	for (i = 0; i < MAX_FLATNUM; i++)
		fl[i].fType = FLAT_NONE;
	i = 0;

	/* F I N D   M I N   P S I   F O R   A  M A G N E T I C   A X I S */
	PsiFlat = Psi[2][2];
	for (ix = 2; ix < nmax - 2; ix++)
		for (iz = 2; iz < nmax - 2; iz++)
			if ((Psi[ix][iz] < PsiFlat) && (Psi[ix][iz] < Psi[ix][iz + 1]) &&
				(Psi[ix][iz] < Psi[ix][iz - 1]) && (Psi[ix][iz] < Psi[ix + 1][iz]) &&
				(Psi[ix][iz] < Psi[ix - 1][iz])) {
				PsiFlat = Psi[ix][iz];
				xFlat = pg->X[ix];
				zFlat = pg->Z[iz];
				FoundAxis = TRUE;
			}
	if (FoundAxis) {
		fl[i].fType = FLAT_AXIS;
		fl[i].X = xFlat;
		fl[i].Z = zFlat;
		if (PsiFlat < 0.0)
			fl[i].Psi = PsiFlat * SMALL_POSITIVE_FACTOR;
		else if (PsiFlat > 0.0)
			fl[i].Psi = PsiFlat * SMALL_NEGATIVE_FACTOR;
		else					/* PsiFlat = 0.0 */
			fl[i].Psi = SMALL_SMALL_FACTOR;
		i++;
		#ifndef IDL
		printf("		[PsiMin found at (X = %.3f, Z = %.3f) Psi = %g]\n", xFlat, zFlat, PsiFlat);
		fprintf(LogFile, "		[PsiMin found at (X = %.3f, Z = %.3f) Psi = %g]\n", xFlat, zFlat, PsiFlat);
		#endif
	}
	/* F I N D   A L L   A X E S   A N D   S E P S */
	for (ix = 2; ix < nmax - 2; ix++)
		for (iz = 2; iz < nmax - 2; iz++) {
			IsMagAxis = FALSE;
			IsSeparatrix = FALSE;
			dPsiX4[TOP] = GET_DX(ix, iz + 1);
			dPsiX4[BOT] = GET_DX(ix, iz - 1);
			dPsiX4[IN] = GET_DX(ix - 1, iz);
			dPsiX4[OUT] = GET_DX(ix + 1, iz);
			dPsiZ4[TOP] = GET_DZ(ix, iz + 1);
			dPsiZ4[BOT] = GET_DZ(ix, iz - 1);
			dPsiZ4[IN] = GET_DZ(ix - 1, iz);
			dPsiZ4[OUT] = GET_DZ(ix + 1, iz);
			dPsiX = GET_DX(ix, iz);
			dPsiZ = GET_DZ(ix, iz);
			IsMinimum = TRUE;
			grad = SQUARE(dPsiX) + SQUARE(dPsiZ);
			for (s = TOP; s <= OUT; s++) {
				if (grad > SQUARE(dPsiX4[s]) + SQUARE(dPsiZ4[s]))
					IsMinimum = FALSE;
			}
			if (!IsMinimum)
				continue;

MULTI;


			/* fourth order finite differences */
			dPsiX = GET_DX_4(ix, iz);
			dPsiZ = GET_DZ_4(ix, iz);
			dPsiXX = (-Psi[ix + 2][iz] + 16.0 * (Psi[ix + 1][iz] + Psi[ix - 1][iz])
					  - 30.0 * Psi[ix][iz] - Psi[ix - 2][iz]) / 12.0;
			dPsiXX /= SQUARE(dx);
			dPsiZZ = (-Psi[ix][iz + 2] + 16.0 * (Psi[ix][iz + 1] + Psi[ix][iz - 1])
					  - 30.0 * Psi[ix][iz] - Psi[ix][iz - 2]) / 12.0;
			dPsiZZ /= SQUARE(dz);
			dPsiXZ = get_dxdz(Psi, ix, iz) / dx / dz;
			det = dPsiXX * dPsiZZ - SQUARE(dPsiXZ);
			if (det > 0.0)
				IsMagAxis = TRUE;
			if (det < 0.0)
				IsSeparatrix = TRUE;
			if (IsMagAxis || IsSeparatrix) {
				xFlat = (dPsiXZ * dPsiZ - dPsiZZ * dPsiX) / det;
				zFlat = (dPsiXZ * dPsiX - dPsiXX * dPsiZ) / det;
				PsiFlat = Psi[ix][iz] + xFlat * (dPsiX + 0.5 * xFlat * dPsiXX + zFlat * dPsiXZ) +
					zFlat * (dPsiZ + 0.5 * zFlat * dPsiZZ);
				if ((fabs(xFlat) > dx) || (fabs(zFlat) > dz)) {
					IsMagAxis = FALSE;
					IsSeparatrix = FALSE;
				}
				xFlat += pg->X[ix];
				zFlat += pg->Z[iz];
			}
			if ((IsMagAxis || IsSeparatrix) && (i == MAX_FLATNUM)) {
				printf("ERROR:	Too many separatrices and axes!\n");
				fprintf(LogFile, "ERROR:	Too many separatrices and axes!\n");
				break;
			}
			if (IsMagAxis) {
				FoundAxis = TRUE;
				printf("		[Axis found at (X = %.3f, Z = %.3f) Psi = %g]\n", xFlat, zFlat, PsiFlat);
				fprintf(LogFile, "		[Axis found at (X = %.3f, Z = %.3f) Psi = %g]\n", xFlat, zFlat, PsiFlat);
			}
			if (IsSeparatrix) {
				fl[i].fType = FLAT_SEP;
				printf("		[Sep  found at (X = %.3f, Z = %.3f) Psi = %g]\n", xFlat, zFlat, PsiFlat);
				fprintf(LogFile, "		[Sep  found at (X = %.3f, Z = %.3f) Psi = %g]\n", xFlat, zFlat, PsiFlat);
			}
			if (IsMagAxis || IsSeparatrix) {
				fl[i].Psi = PsiFlat;
				fl[i].X = xFlat;
				fl[i].Z = zFlat;
				i++;
			}
		}
	if (!FoundAxis)
		nrerror("ERROR:	No magnetic axis found.");
}

/*
**
**	FindValidFlats
**
**	This function creates an array filled with information
**	concerning the location of all separatrices and all
**	magnetic axes.
**
**	A separatrix or magnetic axis is defined when
**		|�p/�x|^2 + |�p/�z|^2 < |�p/�x|^2 + |�p/�z|^2 on all sides
**
**	The minimum is a magnetic axis when
**		�2p/�x2 * �2p/�z2 - (�2p/�z�x)^2 > 0
**
**	The minimum is a separatrix when
**		�2p/�x2 * �2p/�z2 - (�2p/�z�x)^2 < 0
**
**	Interpolation is used to find the actual coordinates for the
**	where |grad(Psi)|^2 vanishes.
**
**	INPUTS:
**		PSIGRID		*pg		==> a pointer to the PSIGRID data
**		FLAT		*fl		==> array of FLATs to hold results
**
**
**	DTG 12/9/98... I don't really care about magnetic axis and I
**  only want to find seperatricies in valid regions.
*/
void          FindValidFlats(TOKAMAK * td, FLAT * fl)
{
    int           ix, iz, i, s, is, nmax;
    double      **Psi;
    double        dPsiX4[4], dPsiZ4[4];
    double        dPsiZ, dPsiX, dPsiXX, dPsiZZ, dPsiXZ, grad, det;
    int           IsMagAxis, IsSeparatrix, IsMinimum;
    int           FoundAxis = FALSE;
    double        xFlat, zFlat, PsiFlat;
    double        dx, dz, tdx, tdz;
    PSIGRID		*pg;

    pg = td->PsiGrid;
    Psi = pg->Psi;
    nmax = pg->Nsize;
    dx = pg->dx;
    dz = pg->dz;
    tdx = 0.5 / dx;
    tdz = 0.5 / dz;

    for (i = 0; i < MAX_FLATNUM; i++)
        fl[i].fType = FLAT_NONE;
    i = 0;

    /* F I N D   M I N   P S I   F O R   A  M A G N E T I C   A X I S */
    PsiFlat = Psi[2][2];
    for (ix = 2; ix < nmax - 2; ix++) {
        for (iz = 2; iz < nmax - 2; iz++) {
            if ((Psi[ix][iz] < PsiFlat) && (Psi[ix][iz] < Psi[ix][iz + 1]) &&
                (Psi[ix][iz] < Psi[ix][iz - 1]) && (Psi[ix][iz] < Psi[ix + 1][iz]) &&
                (Psi[ix][iz] < Psi[ix - 1][iz])) {
                PsiFlat = Psi[ix][iz];
                xFlat = pg->X[ix];
                zFlat = pg->Z[iz];
                FoundAxis = TRUE;
            }
        }
    }
    if (FoundAxis) {
        fl[0].fType = FLAT_AXIS;
        fl[0].X = xFlat;
        fl[0].Z = zFlat;
        if (PsiFlat < 0.0)
            fl[0].Psi = PsiFlat * SMALL_POSITIVE_FACTOR;
        else if (PsiFlat > 0.0)
            fl[0].Psi = PsiFlat * SMALL_NEGATIVE_FACTOR;
        else					/* PsiFlat = 0.0 */
            fl[0].Psi = SMALL_SMALL_FACTOR;

        printf("		[PsiMin found at (X = %.3f, Z = %.3f) Psi = %g]\n", xFlat, zFlat, PsiFlat);
        fprintf(LogFile, "		[PsiMin found at (X = %.3f, Z = %.3f) Psi = %g]\n", xFlat, zFlat, PsiFlat);

    }
    i=1;
    MULTI;

    /* F I N D   A L L   A X E S   A N D   S E P S */
    for (ix = 2; ix < nmax - 2; ix++) {
        for (iz = 2; iz < nmax - 2; iz++) {
            IsMagAxis = FALSE;
            IsSeparatrix = FALSE;
            dPsiX4[TOP] = GET_DX(ix, iz + 1);
            dPsiX4[BOT] = GET_DX(ix, iz - 1);
            dPsiX4[IN] = GET_DX(ix - 1, iz);
            dPsiX4[OUT] = GET_DX(ix + 1, iz);
            dPsiZ4[TOP] = GET_DZ(ix, iz + 1);
            dPsiZ4[BOT] = GET_DZ(ix, iz - 1);
            dPsiZ4[IN] = GET_DZ(ix - 1, iz);
            dPsiZ4[OUT] = GET_DZ(ix + 1, iz);
            dPsiX = GET_DX(ix, iz);
            dPsiZ = GET_DZ(ix, iz);
            IsMinimum = TRUE;
            grad = SQUARE(dPsiX) + SQUARE(dPsiZ);
            for (s = TOP; s <= OUT; s++) {
                if (grad > SQUARE(dPsiX4[s]) + SQUARE(dPsiZ4[s]))
                    IsMinimum = FALSE;
            }
            if (!IsMinimum)
                continue;
            /* fourth order finite differences */
            dPsiX = GET_DX_4(ix, iz);
            dPsiZ = GET_DZ_4(ix, iz);
            dPsiXX = (-Psi[ix + 2][iz] + 16.0 * (Psi[ix + 1][iz] + Psi[ix - 1][iz])
                      - 30.0 * Psi[ix][iz] - Psi[ix - 2][iz]) / 12.0;
            dPsiXX /= SQUARE(dx);
            dPsiZZ = (-Psi[ix][iz + 2] + 16.0 * (Psi[ix][iz + 1] + Psi[ix][iz - 1])
                      - 30.0 * Psi[ix][iz] - Psi[ix][iz - 2]) / 12.0;
            dPsiZZ /= SQUARE(dz);
            dPsiXZ = get_dxdz(Psi, ix, iz) / dx / dz;
            det = dPsiXX * dPsiZZ - SQUARE(dPsiXZ);
            if (det > 0.0)
                IsMagAxis = TRUE;

            if (det < 0.0)
                IsSeparatrix = TRUE;
            if (IsMagAxis || IsSeparatrix) {
                xFlat = (dPsiXZ * dPsiZ - dPsiZZ * dPsiX) / det;
                zFlat = (dPsiXZ * dPsiX - dPsiXX * dPsiZ) / det;
                PsiFlat = Psi[ix][iz] + xFlat * (dPsiX + 0.5 * xFlat * dPsiXX + zFlat * dPsiXZ) +
                    zFlat * (dPsiZ + 0.5 * zFlat * dPsiZZ);
                if ((fabs(xFlat) > dx) || (fabs(zFlat) > dz)) {
                    IsMagAxis = FALSE;
                    IsSeparatrix = FALSE;
                }
                xFlat += pg->X[ix];
                zFlat += pg->Z[iz];
            }
            if ((IsMagAxis || IsSeparatrix) && (i == MAX_FLATNUM)) {
                printf("ERROR:	Too many separatrices and axes!\n");
                fprintf(LogFile, "ERROR:	Too many separatrices and axes!\n");
                break;
            }
            if (IsMagAxis) {
                if ((!FoundAxis) || (PsiFlat < fl[0].Psi)) { /* found new or better axis! */
                    FoundAxis = TRUE;
                    fl[0].Psi = PsiFlat;
                    fl[0].X = xFlat;
                    fl[0].Z = zFlat;
                    fl[0].fType=FLAT_AXIS;
                }
            }
            if (IsSeparatrix) {
                for (is = 0; is < td->NumSeps; is++) {
                    SEPARATRIX *sp = td->Seps[is];
                    if (!sp->Enabled) continue;
                    if (IsValidSeparatrix(sp,xFlat,zFlat)) {
                        fl[i].fType = FLAT_SEP;
                        fl[i].Psi = PsiFlat;
                        fl[i].X = xFlat;
                        fl[i].Z = zFlat;

                        sp->IsSeparatrix = TRUE;
                        if (PsiFlat < sp->PsiSep ) {
                            sp->PsiSep 	= PsiFlat;
                            sp->Xs 	= xFlat;
                            sp->Zs 	= zFlat;
                        }

                        i++;

                        printf("		[Sep  found at (X = %.3f, Z = %.3f) Psi = %g]\n",
                               xFlat, zFlat, PsiFlat);
                        fprintf(LogFile, "		[Sep  found at (X = %.3f, Z = %.3f) Psi = %g]\n",
                                xFlat, zFlat, PsiFlat);

                        break;   // only find first one...
                    }
                }
                if (is == td->NumSeps) {
                    printf("		[Unmatched Sep  found at (X = %.3f, Z = %.3f) Psi = %g]\n",
                           xFlat, zFlat, PsiFlat);
                    fprintf(LogFile, "		[Unmatched Sep  found at (X = %.3f, Z = %.3f) Psi = %g]\n",
                            xFlat, zFlat, PsiFlat);

                }
            }
        }
    }
    if (FoundAxis) {
        printf("		[Axis found at (X = %.3f, Z = %.3f) Psi = %g]\n",fl[0].X, fl[0].Z, fl[0].Psi);
        fprintf(LogFile, "		[Axis found at (X = %.3f, Z = %.3f) Psi = %g]\n", fl[0].X, fl[0].Z, fl[0].Psi);
    } else
        nrerror("ERROR:	No magnetic axis found.");
}

/*
**	MarkDivertors
**
**
**	This function examines all valid separatrices and sets
**	a flag in the divertor array if the location is within a
**	divertor region of a separatrix.
*/
void          MarkDivertors(TOKAMAK * td, int **div)
{
	int           i, ix, iz, nmax;
	PSIGRID      *pg;
	SEPARATRIX   *sp;
	double       *X, *Z;
	double        xa, za;
MULTI;


	pg = td->PsiGrid;
	nmax = pg->Nsize;
	X = pg->X;
	Z = pg->Z;
	xa = pg->XMagAxis;
	za = pg->ZMagAxis;

	for (ix = 0; ix <= nmax; ix++) {
MULTI;


		for (iz = 0; iz <= nmax; iz++) {
			div[ix][iz] = FALSE;
			for (i = 0; i < td->NumSeps; i++) {
				sp = td->Seps[i];
				if (sp->IsSeparatrix && IsPtDivertor(sp, X[ix], Z[iz], xa, za))
					div[ix][iz] = TRUE;
			}
		}
	}

#ifdef DEBUG_DIVERTORS
        printf("IsDiverter\n");
        for (iz = nmax; iz >= 0; iz--) {
            for (ix = 0; ix <= nmax; ix++)
                printf("%1d", div[ix][iz]);
            printf( "\n");
        }
        printf("\n");
#endif
}

/*
**	IsTruePlasma
**
**
*/
int           IsTruePlasma(PSIGRID * pg, int ixp, int izp, int ixa, int iza)
{
	double      **Psi, PsiLim;
	int           i, ix, iz, nmax, delx, delz, del, isTrue;

	nmax = pg->Nsize;
	Psi = pg->Psi;
	PsiLim = pg->PsiLim;
	delx = ixa - ixp;
	delz = iza - izp;
	del = MAX(delx, delz);

	isTrue = TRUE;
	for (i = 1; i <= del; i++) {
		ix = ixp + (int) ceil(i * (double) delx / del);
		iz = izp + (int) ceil(i * (double) delz / del);
		if (Psi[ix][iz] > PsiLim)
			isTrue = FALSE;
	}
	return isTrue;
}

/*
**	FinalCheckIsPlasma
**
**
*/
void          FinalCheckIsPlasma(PSIGRID * pg)
{
	int           ix, iz, nmax, ixa, iza;
	int         **ip;
    MULTI;


	nmax = pg->Nsize;
	ip = pg->IsPlasma;
	ixa = (int) floor((pg->XMagAxis - pg->Xmin) / (pg->dx));
	iza = (int) floor((pg->ZMagAxis - pg->Zmin) / (pg->dz));

	for (ix = 1; ix < nmax; ix++) {

MULTI;

		for (iz = 1; iz < nmax; iz++) {
			if (ip[ix][iz])
				ip[ix][iz] = IsTruePlasma(pg, ix, iz, ixa, iza);
		}
	}
}

/*
**	CheckIsDivertor
**
*/
int           CheckIsDivertor(PSIGRID * pg, int **div, double x, double z)
{
	double        xr, zr;
	int           ix0, ix1, iz0, iz1;

	xr = (x - pg->Xmin) / pg->dx;
	zr = (z - pg->Zmin) / pg->dz;
	ix0 = (int) floor(xr);
	ix1 = (int) ceil(xr);
	iz0 = (int) floor(zr);
	iz1 = (int) ceil(zr);

	return (div[ix0][iz0] || div[ix1][iz0] || div[ix0][iz1] || div[ix1][iz1]);
}

/*
**	FindMinPsiLimiter
**
** 	Finds minimum Psi along the LIMITER and
**	stores it in PsiMin IF NOT DIVERTOR.
**
** 	The idea is to exclude regions that are within a separtrix.
** 	This requires that we choose the minimum point along the LIMITER
** 	that actually has Psi decreasing along a line moving towards the
** 	magnetic axis, (xa,za).


DTG the above idea is pretty good, but not implemented.  Instead
we implement a rule that says that the plasma must be to the left when
traveling from the first limiter point to the second

*/
void          FindMinPsiLimiter(PSIGRID * pg, LIMITER * lm, int **div)
{
	int           i;
	double        xL, zL, Psi;
	double        x1, z1, x2, z2;
	double        dxl, dzl;
	int           imax, nmax;	/* The number of points along a LIMITER to check.*/

	nmax = pg->Nsize;

	x1 = lm->X1;
	z1 = lm->Z1;
	x2 = lm->X2;
	z2 = lm->Z2;

	imax = (int) ceil(sqrt(SQUARE((x2 - x1) / pg->dx) + SQUARE((z2 - z1) / pg->dz)));
	dxl = (x2 - x1) / imax;
	dzl = (z2 - z1) / imax;

	lm->PsiMin = pg->PsiLim;	/* initialize to PsiMax */
	lm->Xmin = x1;
	lm->Zmin = z1;
	/* DTG 6/17/98... this looks like a bug, should start with 0 */
	//for (i = 1; i <= imax; i++) {
	for (i = 0; i <= imax; i++) {
		xL = x1 + i * dxl;
		zL = z1 + i * dzl;
		Psi = GetPsi(pg, xL, zL);
		if ((Psi < lm->PsiMin) &&
		    !CheckIsDivertor(pg, div, xL, zL) &&
		    ((lm->Enabled < 2) || (Psi > GetPsi(pg,xL-dzl,zL+dxl)))) {
			lm->PsiMin = Psi;
			lm->Xmin = xL;
			lm->Zmin = zL;
		}
	}
}

#ifdef DIPOLE
/*
**	FindMaxPsiLimiter
**
** 	Finds maximum Psi along the LIMITER and
**	stores it in PsiMin
**
*/
void          FindMaxPsiLimiter(PSIGRID * pg, LIMITER * lm, int **dummy)
{
	int           i;
	double        xL, zL, Psi;
	double        x1, z1, x2, z2;
	double        dxl, dzl;
	int           imax, nmax;	/* The number of points along a LIMITER to check.*/

	nmax = pg->Nsize;

	x1 = lm->X1;
	z1 = lm->Z1;
	x2 = lm->X2;
	z2 = lm->Z2;

	imax = (int) ceil(sqrt(SQUARE((x2 - x1) / pg->dx) + SQUARE((z2 - z1) / pg->dz)));
	dxl = (x2 - x1) / imax;
	dzl = (z2 - z1) / imax;

	lm->PsiMin = pg->PsiAxis;	/* initialize to PsiMin */
	lm->Xmin = x1;
	lm->Zmin = z1;
	for (i = 0; i <= imax; i++) {
		xL = x1 + i * dxl;
		zL = z1 + i * dzl;
		Psi = GetPsi(pg, xL, zL);
		if (Psi > lm->PsiMin) {
			lm->PsiMin = Psi;
			lm->Xmin = xL;
			lm->Zmin = zL;
		}
	}
}

#endif


/*
**	SetEdgeCurrent
**
**	This will set G2p[0] according to an approximate formula.
**
**		J(edge) = - (Pi/R)*(R B0)^2/DelPsi/MU0 * (dG2/dPsi)
**
**	Note: Jedge is stored as MU0*Jedge.
**
** From J_IsoNoFlow.c.....
**
** 		J(x,z) = -2pi*(X �p/�Psi + ((R0*B0)^2/2X)*�g^2/�Psi)
**
** Therefore, at X = 1,
**
**		J(edge) = - Pi*(R B0)^2 * (dG2/dPsi)
**
*/
void SetEdgeCurrent(TOKAMAK *td)
{
	PLASMA	*pl;
	PSIGRID *pg;

	pg = td->PsiGrid;
	pl = td->Plasma;

	switch (pl->ModelType) {
	  case Plasma_Std:
		  /* no edge current yet */
		  break;
	  case Plasma_IsoNoFlow:
	  	  /* a constant... */
/*		  pl->G2p[0] = - pl->Jedge*(MU0*pl->Ip0/pl->B0/pl->B0/TWOPI); */
		  pl->G2p[0] = - pl->Jedge/(PI*DSQR(pl->B0R0));
		  break;
	  case Plasma_IsoFlow:
		  /* no edge current yet */
		  break;
	  case Plasma_AnisoNoFlow:
		  /* no edge current yet */
		  break;
	  case Plasma_AnisoFlow:
		  /* no edge current yet */
		  break;
	}
}

/*
**	PlasmaBoundary
**
** 	Finds Psi on the limiters and the separtricies.
**
** 	Added S. Jardin's suggestions on Saturday, January 27, 1990 11:06:46 AM.
**
** FIRST, we find any and all magnetic axes and separatrices.
**
** SECOND, we check the limiters...chosing PsiLim to be the smallest
** Psi value that is not within a divertor region.  We do not check
** that there may be a SEPARATRIX within a SEPARATRIX, although the user
** will be able to control where one looks for separatrices.
**
** FINALLY, we set IsPlasma = TRUE whenever Psi < PsiLim and it is
** not within a LIMITER region.
**
** Also, we seem to have trouble finding Seps on the axis.  (Why?)
** For a fix...let's just see if there is a maximum of Psi between the MagAxis
** and the walls.
**
** We can set up a few rules on how the magnetic axis and PLASMA
** boundary are set.
**
** MAGNETIC AXIS:	The magnetic axis is the minimum value of Psi such that
** 					Psi at a LIMITER is still larger.  Moral of story: Don't
** 					put a LIMITER near a coil with a low value of Psi.
**
** PLASMA BOUNDARY:	The PLASMA boundary is the minimum value of a SEPARATRIX
** 					or of a LIMITER...BUT we are careful that the minimum
** 					value of a LIMITER does not occur within a divertor
** 					region.
**
** Steve Jardin's suggestions:
**
** Separatrices:	For a SEPARATRIX, eliminate all PLASMA away from magnetic axis
** 					add to the outside of a line perpendicular to the line from the
** 					axis to the separatirx.
**
** Limiters: 		Make sure that Psi decreases away from the minimum point.  This
** 					insures that the point is not within the divertor region.
**
** These give the STEPS:
**
** 	STEP 1: Initialize
** 	STEP 2: Find all magnetic axes and separatrices
** 	STEP 3: Find minimum Psi along a LIMITER not within a divertor region
** 	STEP 4: Erase PLASMA from divertor regions.
**	STEP 5: Recheck all points for the "bean-ing problem".
**
*/
void          PlasmaBoundary(TOKAMAK * td)
{
	PSIGRID      *pg;
	SEPARATRIX   *sp;
	LIMITER      *lm;
	FLAT          flats[MAX_FLATNUM];
	int           ix, iz, i, nmax;
	double        dx, dz;
	double      **Psi;
	int         **IsDivertor;
	double        PsiMin, PsiMax;
	double        PsiMinL;		/*  PsiMin along a limiter */
	double        PsiMinS;		/* 	PsiMin at a separatrix */
	double        PsiMaxD;   /*  PsiMax on the Dipole coil */
	double 		  PsiAxis;

	/* I N I T I A L I Z E   L O C A L   V A R I A B L E S  */
	pg = td->PsiGrid;
	nmax = pg->Nsize;
	Psi = pg->Psi;
	dx = pg->dx;
	dz = pg->dz;

#ifndef IDL
	printf("INFO:	Plasma Boundary.\n");
	fprintf(LogFile, "INFO:	Plasma Boundary.\n");
#endif

	/* S T E P   1 :    I N I T I A L I Z A T I O N  */
	for (ix = 0; ix <= nmax; ix++)
		for (iz = 0; iz <= nmax; iz++)
			pg->IsPlasma[ix][iz] = FALSE;

	/*                  F I N D   P S I   M I N   &   M A X  */
	PsiMax = PsiMin = Psi[1][1];
        for (ix = 1; ix < nmax; ix++) {
            for (iz = 1; iz < nmax; iz++) {
                if (PsiMax < Psi[ix][iz])
                    PsiMax = Psi[ix][iz];
                if (PsiMin > Psi[ix][iz])
                    PsiMin = Psi[ix][iz];

//			PsiMax = MAX(Psi[ix][iz], PsiMax);
//			PsiMin = MIN(Psi[ix][iz], PsiMin);
            }
        }
	PsiAxis=pg->PsiAxis = pg->PsiLim = PsiMax;

	/*                  I N I T I A L I Z E   S E P S  */
	PsiMinL = PsiMinS = PsiMax;
	for (i = 0; i < td->NumSeps; i++) {
		sp = td->Seps[i];
		sp->PsiSep = PsiMax;
		sp->IsSeparatrix = FALSE;
	}

	/* S T E P   2 :   Check everywhere and see if there is a MagAxis or a Sep */

/*
        FindAllFlats(pg, flats);

	//                 Assign magnetic axes  //
	for (j = 0; j < MAX_FLATNUM; j++)
		if ((flats[j].fType == FLAT_AXIS) && (flats[j].Psi < pg->PsiAxis)) {
			PsiAxis = pg->PsiAxis = flats[j].Psi;
			pg->XMagAxis = flats[j].X;
			pg->ZMagAxis = flats[j].Z;
		}
	td->Plasma->XMagAxis = pg->XMagAxis;
	td->Plasma->ZMagAxis = pg->ZMagAxis;

	//                 Assign separarticies  /
	for (i = 0; i < td->NumSeps; i++) {
		sp = td->Seps[i];
		if (!sp->Enabled)
			continue;
		for (j = 0; j < MAX_FLATNUM; j++)
			if ((flats[j].fType == FLAT_SEP) &&
				IsValidSeparatrix(sp, flats[j].X, flats[j].Z) &&
				(flats[j].Psi < sp->PsiSep)) {
				sp->IsSeparatrix = TRUE;
				sp->PsiSep = flats[j].Psi;
				sp->Xs = flats[j].X;
				sp->Zs = flats[j].Z;
			}
	}

*/
        FindValidFlats(td, flats);

        // first flat is the real magnetic axis
        PsiAxis = pg->PsiAxis = flats[0].Psi;
        td->Plasma->XMagAxis = pg->XMagAxis = flats[0].X;
        td->Plasma->ZMagAxis = pg->ZMagAxis = flats[0].Z;
		td->Plasma->PsiMagAxis = pg->PsiMagAxis = PsiAxis;


        // all other flats have been assigned to separatrices...
        for (i = 1; i < MAX_FLATNUM; i++) {
            if (flats[i].fType == FLAT_SEP) {
                if (PsiMinS > flats[i].Psi)
                    PsiMinS = flats[i].Psi;
            }
	}

	/*                 Mark divertor regions  */
	IsDivertor = imatrix(0, nmax, 0, nmax);
	MarkDivertors(td, IsDivertor);

	/* S T E P   3 :   Find minimum Psi on limiters            			*/
	/*					...being careful to avoid divertor regions.		*/
	/*																	*/
	/* This uses PsiMin as the magnetic axis, and therefore can			*/
	/* eliminate limiters that are within divertor regions.				*/
	PsiMinL = PsiMinS;
	for (i = 0; i < td->NumLimiters; i++) {
		lm = td->Limiters[i];
		if (lm->Enabled > 0) {

MULTI;

			FindMinPsiLimiter(pg, lm, IsDivertor);
			PsiMinL = MIN(PsiMinL, lm->PsiMin);
		}
	}

#ifdef DIPOLE
	/* here's where we have to make a new choice and consider the inner
	** limiters of the floating coil.  I've saved it till last to
	** keep the old checks there... now I want to exclude the innermost region.
	*/

	PsiMaxD  = PsiAxis;

   for (i = 0; i < td->NumLimiters; i++) {
		lm = td->Limiters[i];
		if (lm->Enabled < 0) {

MULTI;


			FindMaxPsiLimiter(pg, lm, IsDivertor);
			PsiMaxD = MAX(PsiMaxD, lm->PsiMin);
		}
	}
#endif

	/* S T E P   4 :   F i l l   I s P l a s m a  */
	/*
	** Note: at this point, PsiMinL contains the PLASMA boundary.
	**
	** As an afterthought, should we follow the prescription in Johnson, et al.,
	** and set the actual PsiLim to be the 0.999 flux surface?
	*/

  /* DTG 6/98 I think this is a bug in the code... It should be the .999
  ** normalized flux surface.  This actually does the opposite of what one
  ** desires...

	// pg->PsiLim = 0.999 * PsiMinL;  */

	pg->PsiLim = PsiMinL;

#ifdef DIPOLE
	pg->PsiAxis = PsiMaxD;
#endif /* DIPOLE */

	pg->DelPsi = pg->PsiLim - pg->PsiAxis;
	td->Plasma->PsiLim = pg->PsiLim;
	td->Plasma->PsiAxis = pg->PsiAxis;

	if (pg->DelPsi <= 0.0)
MULTI;


	for (ix = 0; ix < nmax; ix++)
		for (iz = 0; iz < nmax; iz++)
			pg->IsPlasma[ix][iz] = (Psi[ix][iz] < pg->PsiLim) && !IsDivertor[ix][iz] &&
								   (Psi[ix][iz] > pg->PsiAxis);


	free_imatrix(IsDivertor, 0, nmax, 0, nmax);

	/* F I N A L   C H E C K */
	/*
	** Here, we check every point where IsPlasma, and make sure that
	** Psi doesn't read PsiLim as one moves towards the magnetic axis.
	**
	** This would help, for example, with bean-shaped tokamaks.
	*/

    MULTI;

	FinalCheckIsPlasma(pg);


	#ifndef IDL
	printf("		[PsiLim = %f, PsiAxis = %f]\n", pg->PsiLim, pg->PsiAxis);
	fprintf(LogFile, "		[PsiLim = %f, PsiAxis = %f]\n", pg->PsiLim, pg->PsiAxis);
  #endif
MULTI;

	SetEdgeCurrent(td);
}
