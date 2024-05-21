//  CMeasure.cpp
#include <string>
#include "tokamak.h"
#define MACHINE TOKAMAK*

static const int kMeasBp = 1;


class CMeasure {
	public:
	char 	mName[32];
	int 	mType;	// type of measurement
	double mValue;  /* value of measurement */
	double mStdDev; /* standard deviation of measurement for chisqr fitting */
	double mNow;    /* calculated value from using unknowns */
	double mFit;    /* calculated value from current fit */
	double *L;      /* L vector... dependence of measurement from unknown... one column
						of the chisqr linear matrix */

	CMeasure(int type);		/* constructor */
	virtual ~CMeasure();		/* destructor */
	virtual void CompNow(MACHINE m);
	virtual void CompFit(MACHINE m);
	virtual void CompL(MACHINE m);

};

class CMagMeasure : public CMeasure {
	public:

	double	*CoilGreen;
	double  *ShellGreen;
	double  **PlasmaGreen;
        CMagMeasure(int type) : CMeasure(type) {};
} ;

CMeasure::CMeasure(int type) {
	mType = type;
}
