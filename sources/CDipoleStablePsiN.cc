/* Dipole
**
**  CDipoleStablePsiN.cc
**
**
**  defines a general class for plasma models.  Tries to get out
**  of having switches all over the Tokamak 2.0 code.
**
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"
#include "contour.h"
#include "interpolate.h"
#include "nrSpline.h"
#include "CPlasmaModel.h"

#include "CDipoleStablePsiN.h"
#include "GetPlasmaParameters.h"
#include "rolldown.h"

#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

#if _USE_TSPACK_

#include "tspack.h"

#endif /* _USE_TSPACK_ */

CDipoleStablePsiN::CDipoleStablePsiN(PLASMA *pl)
{
	mModelType = pl->ModelType;
	mNumUnkns = 3;
	mPEdge = 1e4;
	mNpts = 20;
	mRPeak = 1.0;
	mZPeak = 0.0;
	mFracCrit = 0.8;
	mPsiFlat = .95;
	mP = mPsi = mSplPp = mSigmas = NULL;
	mG2pow = 0.0;
}

CDipoleStablePsiN::~CDipoleStablePsiN()
{
	if (mP) 	delete[] mP;
	if (mPsi) 	delete[] mPsi;
	if (mSplPp) delete[] mSplPp;
	if (mSigmas) delete[] mSigmas;

}

static const char * sInputWords[] = {
	"PsiNPeak",
	"PEdge",
	"PsiFlat",
	"NSurf",
	"fCrit",
	"G2Pow"
};

enum sInputCodes {
	kPsiNPeak = 0,
	kPEdge,
	kPsiFlat,
	kNSurf,
	kfCrit,
	kgPower,
	kDone
};

void CDipoleStablePsiN::ModelInput(char *key, char *val, char *init)
{
	int i;
// this recieves statements that one may or may not wish to
// deal with here.
	for (i = 0; i < kDone; i++) {
		if (!strcmp(sInputWords[i], key))
			break;
	}

	switch (i) {
		case kPsiNPeak :
			sscanf(val,"%lf",&mPsiNPeak);
			break;
		case kPEdge :
			sscanf(val,"%lf",&mPEdge);
			break;
		case kPsiFlat :
			sscanf(val,"%lf",&mPsiFlat);
			break;
		case kNSurf :
			sscanf(val,"%d",&mNpts);
			break;
		case kfCrit :
			sscanf(val,"%lf",&mFracCrit);
			break;
		case kgPower :
		    sscanf(val,"%lf",&mG2pow);
		    break;
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


void CDipoleStablePsiN::UpdateModel(TOKAMAK *td)
{
	int i;
	double *rr, *zz, *vv;
	double p, lastp;

	if (td->VacuumOnly != 0) {
		isVacuum = 1;
		return;
	}

	// we'll need the correct B field from Psi
	GetGradPsi(td);

	//mPsiPeak = GetPsi(td->PsiGrid, mRPeak, mZPeak);
	mPsiDipole = td->PsiGrid->PsiAxis;
	mPsiEdge = td->PsiGrid->PsiLim;
	mPsiPeak = mPsiDipole + mPsiNPeak*(mPsiEdge-mPsiDipole);
	if (mP == NULL) {
		mP       = new double[mNpts];
		mPsi     = new double[mNpts];
		mSplPp   = new double[mNpts];
		mSigmas  = new double[mNpts];
	}

	rr = new double[mNpts];
	zz = new double[mNpts];
	vv = new double[mNpts];

	// first lets compute Psi
	mPsi[0]     = mPsiPeak;
	mPsi[mNpts-1] = mPsiEdge;
	for (i=0;i<mNpts;i++) {
		mPsi[i] = mPsiPeak + mPsiFlat*(mPsiEdge-mPsiPeak)*i/(mNpts-1);
	}

	// now lets get some r's and z's to start our contours safely!
//	RollDownHill(td, mRPeak, mZPeak, mNpts, mPsi, rr, zz);

    // now lets compute v for each psi (or contours starting at r, z
//	CompVPoints(td, mNpts, rr, zz, vv);

	//GetRandV(td, mRPeak, mZPeak, mNpts, mPsi, rr, zz, vv);
	GetRandVfromPsi(td, mNpts, mPsi, rr, zz, vv);

	// finally, we define P = Pedge * (Vedge/V)^(frac*5/3)

	lastp = MU0*mPEdge;
	for (i=mNpts-1;i>=0;i--) {
		p = MU0*mPEdge*pow(vv[mNpts-1]/vv[i],mFracCrit*1.666666666666667);
		lastp = mP[i] = (lastp < p) ? p : lastp;
	}

	mP[mNpts-1] = MU0*mPEdge;
	mPPeak = mP[0];

#if _USE_TSPACK_
	{
		double *wk;

		long ncd = 2;
		long per = 0;
		long iendc = 1;
		long unifrm = 0;
		long n = mNpts;
		long ier, lwk;

		if (ncd == 1) {
			lwk = 0;
			wk = NULL;
		} else {
			if (unifrm) {
				if (per) {
					lwk = 2*n-2;
				} else {
					lwk = n-1;
				}
			} else {
				if (per) {
					lwk = 3*n-3;
				} else {
					lwk = 2*n-2;
				}
			}
			wk = new double[lwk];
		}

		mSplPp[0] = 0.0 ;
		mSplPp[mNpts-1] = 0.0;

		// now lets spline it for later and set the end derivatives to zero.
		TSPSI(&n, mPsi, mP, &ncd, &iendc, &per, &unifrm, &lwk, wk, mSplPp, mSigmas, &ier);
		if (ier < 0) {
			printf("Error in tension spline: %d",ier);
			nrerror("Error in tension spline routine.\n");
		}

		if (wk) delete[] wk;
	}
#else

	spline(mPsi-1,mP-1,mNpts,0.0,0.0,mSplPp-1);   // spline likes vectors from 1-n

#endif /* _USE_TSPACK_ */

	delete[] rr;
	delete[] zz;
	delete[] vv;

}

double CDipoleStablePsiN::P(double psi)
{
	double p, theta;

	if (isVacuum) return 0;

	if ((psi <= mPsiDipole) || (psi > mPsiEdge)) {
		p = 0.0;
	} else if (psi <= mPsiPeak) {
	// on the inside just use a simple sine wave
		theta = PI * (psi - mPsiDipole) / (mPsiPeak - mPsiDipole);
		p = mPPeak*(1-cos(theta))/2;
	} else if (psi >= mPsiEdge*mPsiFlat) {
		p = mPEdge*MU0;
	} else {
	// on the outside, use the computed value
#if _USE_TSPACK_
		long ier;
		long n = mNpts;
		p = HVAL(&psi, &n, mPsi, mP, mSplPp, mSigmas, &ier);
		if (ier < 0) nrerror("Error in CDipoleStablePsiN::P\n");
#else
	    splint(mPsi-1,mP-1,mSplPp-1,mNpts,psi,&p);
#endif /* _USE_TSPACK_ */

	}

	return p;
}

double CDipoleStablePsiN::Pp(double psi)
{
	double pp, theta;

	if (isVacuum) return 0;

	if ((psi <= mPsiDipole) || (psi >= mPsiEdge/mPsiFlat)) {
		pp = 0.0;
	} else if (psi <= mPsiPeak) {
		theta = PI * (psi - mPsiDipole) / (mPsiPeak - mPsiDipole);
		pp = mPPeak/2*sin(theta)*PI/(mPsiPeak - mPsiDipole);
	} else {
#if _USE_TSPACK_
		long ier;
		long n = mNpts;
		pp = HPVAL(&psi, &n, mPsi, mP, mSplPp, mSigmas, &ier);
		if (ier < 0) nrerror("Error in CDipoleStablePsiN::P\n");
#else
		double p;
		splint_dervs(mPsi-1,mP-1,mSplPp-1,mNpts,psi,&p,&pp);
#endif /* _USE_TSPACK_ */
	}

	return pp;
}


void CDipoleStablePsiN::ModelDescription(FILE *fi)
{
	fprintf(fi,"    Plasma_DipoleStablePsiN with RPeak, ZPeak, PEdge, fCrit, PsiFlat, G2Pow, NPts = \n");
	fprintf(fi,"    %f %f %f %f %f %f %d\n",
		mRPeak,mZPeak,mPEdge,mFracCrit,mPsiFlat,mG2pow,mNpts);

}

void CDipoleStablePsiN::ModelOutput(FILE *fi)
{
	int i;
	fprintf(fi,"Psi\tP\tPp\tG2\tG2p\t  [%d]\n",mNpts);
	for (i=0;i<mNpts;i++) fprintf(fi,"%f %g %g %g %g\n", 							mPsi[i],P(mPsi[i]),Pp(mPsi[i]),G2(mPsi[i]),G2p(mPsi[i]));
}

double CDipoleStablePsiN::G2(double psi)
{
	double lG2;

	if (isVacuum) return 1.0;

	if (mG2pow == 0.0) return CPlasmaModel::G2(psi);
	// this really is wrong... not sure what I was thinking of 25 years ago.
	if (psi <= mPsiDipole) {
		lG2 = 1.0;
	} else if (psi >= mPsiEdge ) {
		lG2 = 1.0;
	} else {
		lG2 = pow(1.0 - (psi - mPsiDipole) / (mPsiEdge - mPsiDipole), mG2pow);
	}
	return lG2;
}

double CDipoleStablePsiN::G2p(double psi)
{
	double lG2p;

	if (isVacuum) return 0.0;

	if (mG2pow == 0.0) return CPlasmaModel::G2p(psi);

	if ((psi <= mPsiDipole) || (psi >= mPsiEdge )) {
		lG2p = 0.0;
	} else {
		lG2p = -mG2pow*pow(1.0 - (psi - mPsiDipole) / (mPsiEdge - mPsiDipole), mG2pow-1) /
		   (mPsiEdge-mPsiDipole);
	}
	return lG2p;
}

void CDipoleStablePsiN::GetPParam(TOKAMAK *td)
{
	PLASMA       *pl;
	PSIGRID      *pg;
	int           nmax, ix, iz;
	double      **Psi, **Press, **gPsiX, **gPsiZ, **gPsi2, **G, **Bt, **B2;
	double        dx, dz, PsiD, PsiW;
	int         **ip;

	pg = td->PsiGrid;
	pl = td->Plasma;
	nmax = pg->Nsize;
	dx = pg->dx;
	dz = pg->dz;
	PsiW    = pg->PsiLim;
	PsiD    = pg->PsiAxis;
	Psi = pg->Psi;
	ip = pg->IsPlasma;
	Press = pl->Piso;
	gPsiX = pl->GradPsiX;
	gPsiZ = pl->GradPsiZ;
	gPsi2 = pl->GradPsi2;
	G = pl->G;
	Bt = pl->Bt;
	B2 = pl->B2;

	/*  P R E S S U R E,  G,  E T C . . . .  */
	for (ix = 0; ix <= nmax; ix++)
		for (iz = 0; iz <= nmax; iz++) {
			if ((ip[ix][iz]) && (isVacuum==0)){
				Press[ix][iz] = P(Psi[ix][iz]);
			} else {
				Press[ix][iz] = 0.0;
			}
			/* G not effected by pressure */
			G[ix][iz] = sqrt(G2(Psi[ix][iz]));
			Bt[ix][iz] = G[ix][iz] * pl->B0R0 / pg->X[ix];

			B2[ix][iz] = gPsi2[ix][iz] / TWOPI / pg->X[ix] / TWOPI / pg->X[ix];
			B2[ix][iz] = B2[ix][iz] + DSQR(Bt[ix][iz]);
		}
}
