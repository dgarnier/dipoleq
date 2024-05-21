/* Dipole
**
**  CDipoleIntStable.h
**
**
**  defines a general class for plasma models.  Tries to get out
**  of having switches all over the Tokamak 2.0 code.
**
*/

class CDipoleIntStable : public CPlasmaModel {

//	int ModelType;
//	int NumUnkns;      // number of unknowns or free params for Least squares fit

protected:
	double 		mPEdge;
	double 		mRPeak,mZPeak;
	double 		mPsiFlat;
	double 		mFracCrit;
	int 		mNpts;
	double      mPsiPeak;
	double      mPsiDipole;
	double      mPsiEdge;
	double      mPPeak;
	double      mG2pow;
	double      *mP;
	double      *mPsi;
	double      *mSplPp;
	double      *mSigmas;

public:
	CDipoleIntStable(PLASMA *pl);
	virtual double  P(double Psi);
	virtual double  Pp(double Psi);
	virtual double  G2(double Psi);
	virtual double  G2p(double Psi);
	virtual double  ***dJdy(TOKAMAK *) { return NULL; };
	virtual void    UpdateModel(TOKAMAK *td);
	virtual void    ModelOutput(FILE *);
	virtual void 	ModelDescription(FILE *);
	virtual void    ModelInput(char *, char *, char *);
	virtual void    GetPParam(TOKAMAK *);
	virtual ~CDipoleIntStable();
};
