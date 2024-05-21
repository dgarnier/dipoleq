/* Dipole
**
**  CDipoleIntStable.cc
**
**
**  defines a general class for plasma models.  Tries to get out
**  of having switches all over the Tokamak 2.0 code.
**
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "psigrid.h"
#include "plasma.h"
#include "tokamak.h"
#include "CPlasmaModel.h"

#include "CDipoleStd.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06

CDipoleStd::CDipoleStd(PLASMA *pl)
{
	mModelType = pl->ModelType;
	mNumUnkns = 4;

	// unknowns
	mRPeak = 1.0;
	mZPeak = 0.0;
	mPPeak = 1.0;
	mPrExp = -6;

	mPsiPeak = 0;
	mPsiEdge = 0;
	mPsiDipole = 0;
}

static const char * sInputWords[] = {
	"RPeak",
	"ZPeak",
	"PPeak",
	"PrExp"
};

enum sInputCodes {
	kRPeak = 0,
	kZPeak,
	kPress0,
	kPrExp,
	kDone
};

void CDipoleStd::ModelInput(char *key, char *val, char *init)
{
	int i;
// this recieves statements that one may or may not wish to
// deal with here.
	for (i = 0; i < kDone; i++) {
		if (!strcmp(sInputWords[i], key))
			break;
	}

	switch (i) {
		case kRPeak :
			sscanf(val,"%lf",&mRPeak);
			break;
		case kZPeak :
			sscanf(val,"%lf",&mZPeak);
			break;
		case kPress0 :
			sscanf(val,"%lf",&mPPeak);
			break;
		case kPrExp :
			sscanf(val,"%lf",&mPrExp);
		default :
			CPlasmaModel::ModelInput(key,val,init);
	}

//		if (terms != NULL) {
//		tok = strtok(init,",\n");
//		for (i=0; i<maxt; i++) {
//			if (tok == NULL) break;
//			sscanf(tok,"%lf",terms++);
//			tok = strtok(NULL,",\n");
//		}


}


void CDipoleStd::UpdateModel(TOKAMAK *td)
{
	if (td->VacuumOnly != 0) isVacuum = 1;

	mPsiPeak = GetPsi(td->PsiGrid, mRPeak, mZPeak);
	mPsiDipole = td->PsiGrid->PsiAxis;
	mPsiEdge = td->PsiGrid->PsiLim;
}

double CDipoleStd::P(double psi)
{
	double p, theta;

	if (isVacuum) return 0.0;

	if ((psi <= mPsiDipole) || (psi > mPsiEdge)) {
		p = 0;
	} else if (psi <= mPsiPeak) {
		theta = (PI/2) * (psi - mPsiDipole) / (mPsiPeak - mPsiDipole);
		p = mPPeak*sin(theta)*sin(theta);
	} else {
		p = mPPeak*pow( psi/mPsiPeak , mPrExp);
	}

	return MU0*p;
}


double CDipoleStd::Pp(double psi)
{
	double pp, theta;

	if (isVacuum) return 0.0;

	if ((psi <= mPsiDipole) || (psi > mPsiEdge)) {
		pp = 0;
	} if (psi <= mPsiPeak) {
		theta = (PI/2) * (psi - mPsiDipole) / (mPsiPeak - mPsiDipole);
		pp = mPPeak*sin(theta)*cos(theta)*PI/(mPsiPeak - mPsiDipole);
	} else {
		pp = mPPeak/mPsiPeak*mPrExp*pow( psi/mPsiPeak , mPrExp - 1);
	}

	return MU0*pp;
}

void CDipoleStd::ModelDescription(FILE *fi)
{
	fprintf(fi,"    Plasma_DipoleStd with RPeak, ZPeak, PPeak, PExp = %f %f %f %f\n",
		mRPeak,mZPeak,mPPeak,mPrExp);

}

void CDipoleStd::ModelOutput(FILE *fi)
{
	fprintf(fi,"    Plasma_DipoleStd with RPeak, ZPeak, PPeak, PExp = %f %f %f %f\n",
		mRPeak,mZPeak,mPPeak,mPrExp);

}
