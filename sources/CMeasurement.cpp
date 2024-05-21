/*
 *  CMeasurement.cpp
 *  DipolEq
 *
 *  Created by Darren Garnier on 3/31/05.
 *  Copyright 2005 Columbia University. All rights reserved.
 *
 */

#include "CMeasurement.h"
//#include "CMeasFluxLoop.h"
//#include "CMeasBpCoil.h"
//#include "CMeasCoilFlux.h"
//#include "CMeasPressure.h"
#include "nrutil.h"

CMeasurement* CMeasurement::CreateWithType(char* mType)
{
    CMeasurement *mp = NULL;



    return mp;
}

CMeasurement* CMeasurement::CreateWithType(int mType)
{
    CMeasurement *mp = NULL;



    return mp;
}


void CMeasurement::SetParameter(char *key, char *val)
{
    const char * kInputWords[] = {
        "Name",
        "Value",
        "StdDev"
    };

    enum eInputCodes {
        kName,
        kValue,
        kStdDev,
        kDone
    };

    int i;

    for (i=0; i<kDone; i++)
        if (!strcmp(kInputWords[i], key)) break;

    switch (i) {
        case kValue:
            sscanf(val,"%lf",&mValue);
            break;
        case kStdDev:
            sscanf(val,"%lf",&mStdDev);
            break;
        case kName:
            strncpy(mName,val,sizeof(mName)-1);
            break;
        default:
        {
            char err[256];
            sprintf(err,"Measurement: Invalid key-value pair %s = %s\n",key,val);
            nrerror(err);
        }
    }


}
