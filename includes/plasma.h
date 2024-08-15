 /*
** 	TokaMac v2.0
**
** 	plasma.h
**
** 	Interface file for plasma.c
**
** 	File:		plasma.h
** 	Date:		June 19, 1992
**
**	Modifications:
**
**		April 8, 1993		Added flux parameter profile array pointers
**		Sept. 16, 1993		Added virial terms (Lao, Nuc Fusion, 1985)
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _PLASMA_

#define _PLASMA_ 1

struct CPlasmaModel;

#define MaxPolyTerms		7

#define Plasma_Std			0
#define	Plasma_IsoNoFlow	1
#define	Plasma_IsoFlow		2
#define	Plasma_AnisoNoFlow	3
#define	Plasma_AnisoFlow	4
#define Plasma_DipoleStd	5
#define Plasma_DipoleIntStable 6
#define Plasma_DipoleStablePsiN 7


/*
**
** Plasma
**
**	The parameters defining a tokamak
**
*/

typedef struct plasma {
	int					Nsize;				/* The arrays are Nsize x Nsize */
	int					ModelType;			/* profile type */
	struct CPlasmaModel	*Model;
	int 				G2pTerms;
	int					HTerms;
	int 				PpTerms;			/* number of terms in poly. expansion */
	int					RotTerms;
	int					SisoTerms;
	int					SparTerms;
	int					SperTerms;

	/***
	* The next quantities are the polynomial expansions of the flux functions.
	***/

	double				*G2p;		/* gradient of square of toroidal flux */
	double				*H;			/* Bernouilli function */
	double 				*Pp;		/* gradient of isotropic pressure */
	double				*Rot;		/* Toroidal flux flunction */
	double				*Siso;		/* Adiabatic flux function for isotropic pressure */
	double				*Spar;		/* Parallel flux function */
	double				*Sper;		/* Perpendicular flux function */

	/***
	* The next five quantities are the measurable quantities that we can derive from
	* the last five flux functions.
	***/

	double				**B2;				/* Square of the total field strength */
	double				**GradPsiX;			/* �Psi/�x                 */
	double				**GradPsiZ;			/* �Psi/dz                 */
	double				**GradPsi2;			/* |�Psi/�x|2 + |�Psi/dz|2 */
	double              **Bt;               /* toroidal magnetic field */
	double              **G;                /* G = SQRT(G2) is RELATED to the toroidal flux */
	double				**Rho;				/* mass density */
	double				**Piso;				/* isotropic pressure */
	double				**Ppar;				/* parallel pressure */
	double				**Pper;				/* perpendicular pressure */
	double				**Alpha;			/* Alpha = MU0*(Ppar - Pper)/B^2 */


	double 				StndP;				/* when a fixed profile is used */
	double  			StndG;

	double				Jedge;				/* an approximate edge current density */

	double  			ChiSqr;				/* the quality of the fit */

	double  			R0;		/* the initial major radius */
	double				Z0;		/* the initial up/down offset */
	double  			a0;		/* the initial minor radius */
	double  			Ip0;	/* the initial plasma current */
	double  			B0;		/* vacuum magnetic field at R0 */
	double  			B0R0;

	int					NumBndMomts;/* Number of boundary moments (from zero) */
    int 				NumPsiPts;	/* Number of Psi points for profile output */
    double				PsiXmax;	/* Outermost flux surface from 0.0 to 1.0 */
    double				*q_pr;		/* safety factor profile */
    double				*Volp_pr;	/* flux derivative of flux surface volume */
    double				*Vol_pr;	/* flux surface volume */
    double				*S_pr;		/* global shear of q */
    double				*B2_pr;		/* flux surface average of B2 */
    double				*Well_pr;	/* magnetic well */
    double				*J_pr;		/* flux-surface averaged toroidal current */
    double              *Beta_pr;   /* flux tube averaged beta profile. */
    double              *BetaMax_pr;/* max local beta on flux tube. */
    double              *XBetaMax_pr;/* X at max local beta on flux tube. */
    double              *ZBetaMax_pr;/* Z at max local beta on flux tube. */
    double              *BBetaMax_pr;/* B at max local beta on flux tube. */
    double              *BMax_pr;/* max B on flux tube. */
    double              *XBMax_pr;/* X at max B on flux tube. */
    double              *ZBMax_pr;/* Z at max B flux tube. */

	/* calculated at end.. store for python */
	double 				*Psi_pr; 
	double				*PsiX_pr;
	double				*P_pr;
	double				*G_pr;
	double				*Pp_pr;
	double				*G2p_pr;

	double  			Ip;		/* the output plasma current */
	double  			beta0;	/* the vaccuum toroidal beta at R0*/
	double  			beta;	/* the average toroidal beta*/
	double  			betap;	/* the poloidal beta*/
	double  			li;		/* the normalized internal inductance*/
	double  			Ltotal;	/* the total inductance*/
	double  			mu;		/* the normalized diamagnetism */

	double  			Volume;
	double				CrossSection;	/* area of flux surface */
	double				Perimeter;
	double              Diamag;
	double              q0;				/* central safety factor */
	double  			qCircular;		/* 5 a^2 B0 / Rp Ip     */
	double				qStar;			/* qCircular *(1+k^2)/2 */
	double  			XMagAxis;		/* the magnetic axis */
	double  			ZMagAxis;
	double				PsiMagAxis;		/* Psi at magnetic axis */
	double				PsiAxis,PsiLim;	/* copies of values named in PsiGrid */
	double  			HalfWidth;
	double  			Elongation;
	double  			RStar;			/* the major radius as 2V/(perimeter^2) */
	double  			RCentroid;		/* Current centroid */
	double  			RCenter;		/* Center of outer flux surface */
	double  			RSurfaceAvg;	/* Surface average major radius */

	double				R_vr;			/* Energy-weighted radius */
	double				Alpha_vr;		/* elongation parameter */
	double				S1_vr,S2_vr,S3_vr;	/* Virial moments of plasma boundary */

	double  			TotKinEnergy;	/* the total plasma kinetic energy */
	double  			TotMagEnergy;	/* li* (1/4 mu R Ip^2) */
	} PLASMA;

/*
**
** Public subroutines
**
*/

#ifdef __cplusplus
extern "C" {
#endif

PLASMA *new_Plasma();
void free_Plasma(PLASMA *);
void init_Plasma(PLASMA *);

double PlasmaP  (PLASMA *pl, double Psi);
double PlasmaPp (PLASMA *pl, double Psi);
double PlasmaG  (PLASMA *pl, double Psi);
double PlasmaG2p(PLASMA *pl, double Psi);

#ifdef __cplusplus
}
#endif




#endif
