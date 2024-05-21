/*
** TokaMac v2.0
**
** Measurement header file.
** This files defines the basic data structure used
** to store tokamak measurements and their estimated
** errors.
**
** A point of clarification is needed:
** meas->Fit, is the calculated measurement using the present
** value of the plasma current, etc.
** meas->Now, is the calculated measurement using the present
** unknowns.  These may be different since the present current
** density may not be equal to the present unknowns due
** to under-relaxation.
**
** File:		include:measurement.h
** Date:		January 26, 1993
**
** Revisions:
**
**		August 6, 1993		Added perfectly conducting shells
**		August 6, 1993		Added "Now"-->Fit with present unknowns
**		October 27, 1993	Added meas_ppsix
**		November 15, 1993	Added meas_bpangle
**
**
** Example:
**
**		MEAS aBpMeas,aSaddleMeas;
**
**		aBpMeas.mtype = meas_bp;
**		aBpMeas.Value = 2.03;
**		aBpMeas.X = 0.4;
**		aBpMeas.Z = 0.0;
**		aBpMeas.parm.bp.Angle = 0.0;
**
**		aSaddleMeas.mtype = meas_saddle;
**		aSaddleMeas.Value = 0.0;
**		aSaddleMeas.parm.saddle.X1 = 0.5;
**		aSaddleMeas.parm.saddle.X2 = 0.2;
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _MEASUREMENT_

#define _MEASUREMENT_ 1

/*
**
**	MEAS type identifiers
**
*/

#define meas_unk			0	/* some unknown meas type */
#define	meas_bp				1	/* local value of poloidal field */
#define meas_press			2	/* local isotropic pressure */
#define	meas_pperp			3	/* local perp pressure */
#define meas_ppar			4	/* local parallel pressure */
#define meas_flux			5	/* local poloidal flux */
#define meas_saddle			6	/* two-point saddle coil */
#define meas_circle			7	/* TFTR's poloidal field moment */
#define	meas_coilcur		8	/* measured current through a coil set */
#define meas_plasmacur		9	/* plasma current */
#define meas_bt				10	/* local toroidal field measurement */
#define meas_diam			11	/* diamagnetic flux */
#define meas_bangle			12	/* poloidal angle of mag. field (MSE) */
#define	meas_flowt			13	/* local plasma toroidal flow UNUSED */
#define meas_flowp			14	/* local plasma poloidal flow UNUSED */
#define meas_ne				15	/* local electron density */
#define meas_Te				16	/* local electron temperature */
#define meas_Zeff			17	/* local effective ion charge */
#define meas_Ti				18	/* local ion temperature */
#define meas_rho			19	/* local mass density */
#define meas_rot			20	/* local toroidal rotation */
#define meas_ppsix			21	/* local isotropic pressure at a value of Sqrt(PsiX) */
#define meas_bpangle		22	/* angle of magnetic field with respect to horizontal */
#define meas_pnorm			23	/* pressure at Sqrt(PsiX) normalized to P(0) */
#define meas_J0				24	/* Current density on axis, A/m^2 */

/*
**
** Measurement data flags
**
*/

#define CircleType_btcos	1
#define CircleType_brsin	2
#define CircleType_brcos	3


/*
**
** MEAS data structure
**
** Fundemental description of all experimental
** data.
**
*/

typedef struct meas {
    union {

	struct {
	    double Angle;			/* angle in degrees clockwise from vertical */
	}
	bp;					/* meas_bp */

	struct {
	    double X1;
	    double Z1;
	    double X2;
	    double Z2;
	}
	saddle;				/* meas_saddle */

	struct {
	    double Radius;
	    int Number;
	    int CircleType;
	}
	circle;				/* meas_circle */

	struct {
	    int CoilNum;			/* the index for the Coils[i] array */
	}
	coilcur;			/* meas_coilcur */

    }
    parm;

	void		(*FindFit)(void *, void *);			/* Function which computes Fit */
	void		(*FindNow)(void *, void *);			/* Function which computes Now */
	void		(*FindL)(void *, void *, double *);	/* Function which computes L */
	void		(*FindGreen)(void *, void *); 		/* Function which computes Greens */
	double 		*CoilGreen;							/* Only for magnetic measurements */
	double 		*ShellGreen;						/* Only for magnetic measurements */
	double 		**PlasmaGreen;						/* Only for magnetic measurements */
	double		X,Z;								/* for local measurements */
    double 		Value;								/* the measured value */
    double 		StdDev;								/* std deviation */
    double 		Fit;								/* calcuated "best fit" */
    double 		Now;								/* calcuated "best fit" using unknowns */
    char 		Name[32];							/* string identifier */
    int 		mType;								/* measurement type */

} MEAS;

/*
**
** Function Prototypes
**
*/

#ifdef __cplusplus
extern "C" {
#endif

MEAS *new_Measure(int );
void free_Measure(MEAS *,int , int, int );

#ifdef __cplusplus
}
#endif

#endif
