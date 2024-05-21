/* Dipole
**
**  CDipoleStd.h
**
**
**  defines a general class for plasma models.  Tries to get out
**  of having switches all over the Tokamak 2.0 code.
**
*/

class CDipoleStd : public CPlasmaModel {

//	int ModelType;
//	int NumUnkns;      // number of unknowns or free params for Least squares fit

protected:
	double	mRPeak;
	double  mZPeak;
	double  mPsiPeak;
	double  mPsiEdge;
	double  mPsiDipole;
	double  mPPeak;
	double  mPrExp;

public:
	CDipoleStd(PLASMA *pl);
	virtual void    UpdateModel(TOKAMAK *td);
	virtual double  P(double Psi);
	virtual double  Pp(double Psi);
//	virtual void 	FindJ(TOKAMAK *td, double **J);
//	virtual double 	FindJ_Loc(TOKAMAK *td, int ix, int iz);
	virtual double  ***dJdy(TOKAMAK *) { return NULL; };
	virtual void    ModelOutput(FILE *);
	virtual void 	ModelDescription(FILE *);
	virtual void    ModelInput(char *key, char *val, char *init);
};
