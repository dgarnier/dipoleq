/*
** TokaMac v2.0
**
** psigrid.h
**
** Interface file for PsiGrid.c
**
** File:		psigrid.h
** Date:		June 19, 1992
**
**	USAGE:
**		PSIGRID *pg;
**		pg = New_PsiGrid();
**		pg->Nsize = 64;
**		init_PsiGrid(pg);
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _PSIGRID_

#define _PSIGRID_ 1

#define UpDownSymmetric 	1
#define NotUpDownSymmetric 	0

#define YesPlasma			1
#define NoPlasma			0


/*
**
** PsiGrid Data Structure
**
*/

typedef	struct psigrid { 	/* the poloidal flux and the toroidal current */
	int 		Nsize;			/* size of the computational grid */
	int			Symmetric;		/* if non-zero, then up-down symmetric */
	double 		MaxRes;			/* the maximum residual */
	double 		PastMaxRes; 	/* fMaxRes from previous iteration */
	double		Xmax,Xmin;		/* the computational extent */
	double		Zmax,Zmin;
	double 		dx;				/* radial grid spacing */
	double		dz;				/* vertical grid spacing */
	double		BoundError;		/* size of error in boundary */
	double		BoundThreshold;	/* how small to make step change in boundary */
	double		ResThreshold;	/* how small to make residual */
	double		UnderRelax1;	/* when = 0.0, then update Current linearly */
	double		UnderRelax2;	/* when = 0.0, then update Unknowns linearly */

	double		PsiAxis;		/* value of Psi at Magnetic Axis or FCFS */
	double 		PsiMagAxis;		/* value of Psi at Magnetic Axis */
	double		PsiLim;			/* value of Psi at Plasma/Vacuum boundary */
	double		DelPsi;			/* PsiLim - PsiAxis */
	double		XMagAxis;		/* Copies of data within PLASMA */
	double		ZMagAxis;

	double      *X;				/* the major radial coordinate */
	double      *Z;				/* the major radial coordinate */
	int			**IsPlasma;		/* when == YesPlasma, then grid point contains plasma */
	double 		**Psi; 			/* poloidal flux */
	double 		**Current;	 	/* current density */
	double 		**Residual; 	/* residual */
	} PSIGRID;

/*
**
**	Public function prototypes
**
*/

#ifdef __cplusplus
extern "C" {
#endif


 PSIGRID *new_PsiGrid();
 void init_PsiGrid(PSIGRID *);
 void free_PsiGrid(PSIGRID *);

 void MakePsiSymmetric(PSIGRID *);
 void GetNewResidual(PSIGRID *);
 void NewSolution(PSIGRID *);
 void NewMSolution(PSIGRID *);
 void GoPDE(PSIGRID *);

 double GetPsi(PSIGRID *,double ,double );
 double GetIsPlasma(PSIGRID *,double ,double );

#ifdef __cplusplus
}
#endif

#endif
