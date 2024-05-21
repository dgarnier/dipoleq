/* Dipole
**
**  CPlasmaModel.h
**
**
**  defines a general class for plasma models.  Tries to get out
**  of having switches all over the Tokamak 2.0 code.
**
*/

#ifndef _CPlasmaModel_

#define _CPlasmaModel_ 1

#include <stdio.h>
#include "tokamak.h"
#include "plasma.h"

class CPlasmaModel {
  public:
	int mModelType;
	int mNumUnkns;
	int isVacuum;
	CPlasmaModel() { isVacuum = 0; } ;
	static  CPlasmaModel * CreateModel(PLASMA *p);
	virtual void    UpdateModel(TOKAMAK *);
	virtual double  P(  double ) { return 0.0; };
	virtual double  Pp( double ) { return 0.0; };
	virtual double  G2( double ) { return 1.0; };
	virtual double  G2p(double ) { return 0.0; };
	virtual void 	FindJ(TOKAMAK *td, double **J);
	virtual double 	FindJ_Loc(TOKAMAK *td, int ix, int iz);
	virtual double  ***dJdy(TOKAMAK *td) = 0;
	virtual void    ModelOutput(FILE *) = 0;
	virtual void    ModelDescription(FILE *) = 0;
	virtual void	ModelInput(char *, char *, char *) {};
	virtual void    GetPParam(TOKAMAK *);
	virtual ~CPlasmaModel() {};
};

#endif
