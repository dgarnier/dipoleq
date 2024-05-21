/*
** TokaMac v2.0
**
** Find_ShellCurrent.c
**
** This file contains a subroutine which calculates the currents
** which flow in a perfectly conducting shell given the plasma
** and coil currents.
**
** It uses a method of minimization of magnetic energy. That is
** a variational approach.
**
** In MKS, the magnetic energy of a shell is given by
**
**		Wm = SUM (Psi[i] Cur[i] / 2) + lam1 SUM(Cur[i]) + lam2 SUM(Cur[i]) + ...
**
** where "lam1" is a Lagrange multiplier for the variational problem
** used to insure that SUM(Cur[i]) = 0 for the 1st shell, "lam2" is for
** the 2nd shell, etc.
**
** Psi[i] represents the poloidal flux at the ith SUBSHELL. The value
** of the self-inductance is assumed to be
**
**		L = Xs ( ln(8 Xs/Rs) - 7/4 ) == M[i,i]
**
** Note, we must use the same conventions as Johnson, et al. (i.e. the
** PEST sign conventions.)  This is opposite of the usual conventions
** used in Jeff Freidberg's textbook.
**
** Note, also, that the shell currents have been multiplied by MU0
** as is the convention in TokaMac.
**
** The mutual inductance is given by the Green's function
**
**		M[i,j] = - G(Xi, Zi, Xj, Zj)
**
** which is symmetric.  (The negative sign relates to our sign convention.
** The energy must be positive definite.)
**
** Notice that Psi[i] is the sum of the other (and self) shell currents
** and the induced fluxes from the plasma and coils.  We write
**
**		Psi[i] = Psi[shells] + Psi[plasma] + Psi[coils]
**
** The last two items are known. The first item is expressed in terms of the
** inductance matrix, M...
**
**	 	Psi[shells, i] = SUM( M[i,j] Cur[j] )
**
**	Finally, a note about Up/Down Symmetry.
**	When td->PsiGrid->Symmetric == UpDownSymmetric, then each subshell has a
**  mirror image at -(subshell->Z).  This changes each term of the inductance
** 	matrix!
**
** File:		Find_ShellCurrent.c
** Date:		August 5, 1993
**
** Revisions:
**	Feb, 1996		Changed factor of two error in shell inductances
**				REPORTED by Ned Eisner
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdlib.h>
#include "VAX_Alloc.h"
#include <stdio.h>
#include <math.h>
#include "nrutil.h"
#include "stdio_dmatrix.h"
#include "ludcmp.h"
#include "green.h"
#include "psigrid.h"
#include "coil.h"
#include "shell.h"
#include "tokamak.h"
#include "LoadShellGreens.h"
#include "Find_ShellCurrent.h"

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif



#define cchk(c,n)	if ((c) != ((size_t) n)) nrerror("ERROR: Could not read/write Inductance Matrix.")

extern FILE  *LogFile;

/*
**	Fill_InductanceMatrix
**
**
**	This computes the inductance matrix for the
**	perfectly conducting shells.
**
*/
void          Fill_InductanceMatrix(TOKAMAK * td, double **M, int msize)
{
#pragma unused( msize )
	int           i, j, is, js, im, jm;
	int           sym = 0;
	SHELL        *shell, *shell_j;
	SUBSHELL     *subshell, *subshell_j;

	printf("INFO:	New_InductanceMatrix for %d perfectly conducting shells.\n",
		   td->NumShells);
	fprintf(LogFile, "INFO:	New_InductanceMatrix for %d perfectly conducting shells.\n",
			td->NumShells);

	sym = td->PsiGrid->Symmetric;

	/* F I L L   S H E L L   S E L F   I N D U C T A N C E S */
	im = 1;
	for (i = 0; i < td->NumShells; i++) {
		shell = td->Shells[i];
		for (is = 0; is < shell->NumSubShells; is++) {
			subshell = shell->SubShells[is];
			M[im][im] = Self_Inductance(subshell);
			if (sym == UpDownSymmetric)
				M[im][im] += -Green(subshell->X, subshell->Z, subshell->X, -(subshell->Z));
			im++;
		}
	}

	/* F I L L   S H E L L   M U T U A L   I N D U C T A N C E S */
	/*
	** Notice that M[i][j] = M[j][i].
	*/
	im = 1;
	for (i = 0; i < td->NumShells; i++) {
		shell = td->Shells[i];
		for (is = 0; is < shell->NumSubShells; is++) {
			subshell = shell->SubShells[is];
			jm = 1;
			for (j = 0; j < td->NumShells; j++) {
				shell_j = td->Shells[j];
				for (js = 0; js < shell_j->NumSubShells; js++) {
					subshell_j = shell_j->SubShells[js];
					if (im < jm) {
						M[im][jm] = -Green(subshell->X, subshell->Z, subshell_j->X, subshell_j->Z);
						if (sym == UpDownSymmetric)
							M[im][jm] += -Green(subshell->X, subshell->Z, subshell_j->X, -(subshell_j->Z));
					} else if (im > jm)
						M[im][jm] = M[jm][im];
					jm++;
				}
			}
			im++;
		}
	}

	/* F I L L    T H E    L A G R A N G I A N     M U L T I P L I E R S   */
	/*
	** im now contains the location of the 1st Lagrange mutliplier.
	** This should be NumSubShells + 1.
	*/
	for (i = 0; i < td->NumShells; i++) {
		shell = td->Shells[i];
		for (is = 0; is < shell->NumSubShells; is++)
			M[im][is + 1] = M[is + 1][im] = 1.0;
		im++;					/* im now points to the next Lagrange multiplier */
	}

}

/*
**
**	RewriteShellInductance
**
*/
void          RewriteShellInductance(TOKAMAK * td, double **M, int msize)
{
	FILE         *fi;
	size_t        c;

	fi = fopen(td->SMname, "wb");
	if (!fi)
		nrerror("ERROR:	Could not open file to write Shell Inductance.");

	Fill_InductanceMatrix(td, M, msize);

	/* W R I T E    S H E L L   I N D U C T A N C E */
	c = fwrite_dmatrix(M, 1, msize, 1, msize, fi);
	cchk(c, (msize * msize));

	fclose(fi);

	td->SInductStatus = ShellInductOK;
}

/*
**
**	LoadShellInductance
**
*/
void          LoadShellInductance(TOKAMAK * td, double **M, int msize)
{
	FILE         *fi;
	size_t        c;

	if (td->SInductStatus == ShellInductOK) {	/* read ShellInduct from file */
		fi = fopen(td->SMname, "rb");
		if (!fi)
			nrerror("ERROR:	Could not open shell inductance file for reading.");

		/* R E A D    S H E L L   I N D U C T */
		c = fread_dmatrix(M, 1, msize, 1, msize, fi);
		cchk(c, (msize * msize));

		fclose(fi);
	} else
		RewriteShellInductance(td, M, msize);
}

/*
**	Find_ShellCurrent
**
**
**	This solves the equation
**
**		M * Ishell = B
**
**	for Ishell where M is the inductance matrix and
**	B is the known current drivers.
**
*/
void          Find_ShellCurrent(TOKAMAK * td)
{
	int           i, j, is, im, jm, msize;
	int           num_subshells = 0;
	int           ix, iz, nmax;
	double      **J, *X, *Z, sum;
	int         **ip;
	double      **pg, *cg;		/* PlasmaGreen and CoilGreen */
	COIL         *coil;
	SHELL        *shell;
	SUBSHELL     *subshell;

	double      **M;			/* the inductance matrix for all subshells */

	double       *B;			/* the subshell "drivers" from the plasma and coils */

	double        LUD_D;
	int          *LUD_INDX;

	J = td->PsiGrid->Current;
	X = td->PsiGrid->X;
	Z = td->PsiGrid->Z;
	ip = td->PsiGrid->IsPlasma;
	nmax = td->PsiGrid->Nsize;

	printf("INFO:	Find_ShellCurrent for %d perfectly conducting shells.\n",
		   td->NumShells);
	fprintf(LogFile, "INFO:	Find_ShellCurrent for %d perfectly conducting shells.\n",
			td->NumShells);

	/* D E T E R M I N E   S I Z E   O F   P R O B L E M  */
	for (i = 0; i < td->NumShells; i++) {
		shell = td->Shells[i];
		num_subshells += shell->NumSubShells;
	}

	/* the matrix size is the sum of the inductances plus Lagrange multipliers */
	msize = num_subshells + td->NumShells;

	/* S E T U P   A R R A Y S */
	M = dmatrix(1, msize, 1, msize);
	for (im = 1; im <= msize; im++)
		for (jm = 1; jm <= msize; jm++)
			M[im][jm] = 0.0;

	LoadShellInductance(td, M, msize);

	B = dvector(1, msize);
	for (im = 1; im <= msize; im++)
		B[im] = 0.0;

	LUD_INDX = ivector(1, msize);

	/* F I N D   K N O W N   F L U X   D R I V E R S */
	/*
	** The flux "drivers" for a shell are the plasma currents and the
	** coil currents.
	*/

	LoadShellGreens(td);

	/* First, the plasma current */
	printf("		[Plasma - Shell Coupling]\n");
	fprintf(LogFile, "		[Plasma - Shell Coupling]\n");
	im = 1;
	for (i = 0; i < td->NumShells; i++) {
		shell = td->Shells[i];
		for (is = 0; is < shell->NumSubShells; is++) {
			subshell = shell->SubShells[is];
			pg = subshell->PlasmaGreen;
			sum = 0.0;
			for (ix = 1; ix < nmax; ix++)
				for (iz = 1; iz < nmax; iz++)
					sum += pg[ix][iz] * J[ix][iz];
			B[im] -= sum * (td->PsiGrid->dx) * (td->PsiGrid->dz);
			im++;
		}
	}

	/* Then, the coil sets */
	printf("		[Coil - Shell Coupling]\n");
	fprintf(LogFile, "		[Coil - Shell Coupling]\n");
	im = 1;
	for (i = 0; i < td->NumShells; i++) {
		shell = td->Shells[i];
		for (is = 0; is < shell->NumSubShells; is++) {
			subshell = shell->SubShells[is];
			cg = subshell->CoilGreen;
			for (j = 0; j < td->NumCoils; j++) {
				coil = td->Coils[j];
				if (coil->Enabled)
					B[im] -= cg[j] * coil->CoilCurrent;
			}
			im++;
		}
	}

	free_ShellGreens(td);

	/* S O L V E   T H E   L I N E A R   S Y S T E M   U S I N G   L U  D E C O M P  */
	ludcmp(M, msize, LUD_INDX, &LUD_D);
	lubksb(M, msize, LUD_INDX, B);	/* the answer is now contained in B */

	/* S E T   T H E   S U B C O I L   C U R R E N T S */
	im = 1;
	for (i = 0; i < td->NumShells; i++) {
		shell = td->Shells[i];
		for (is = 0; is < shell->NumSubShells; is++) {
			subshell = shell->SubShells[is];
			subshell->Current = B[im];
			im++;
		}
	}

	free_ivector(LUD_INDX, 1, msize);
	free_dvector(B, 1, msize);
	free_dmatrix(M, 1, msize, 1, msize);
}
