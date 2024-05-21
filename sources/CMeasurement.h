/*
 *  CMeasurement.h
 *  DipolEq
 *
 *  Created by Darren Garnier on 3/31/05.
 *  Copyright 2005 Columbia University. All rights reserved.
 *
 *  Replaces the measurements.h types...
 */

#ifndef _CMeasurement_h_
#define _CMeasurement_h_ 1

#include <vector>

class CMeasurement;

typedef struct tokamak TOKAMAK;

typedef std::vector<CMeasurement *> CMeasurementList;

class CMeasurement {
public:
    double                  mValue;
    double                  mStdDev;
    double                  mFit;
    double                  mNow;
    char                    mName[32];

    // virtual functions to find fits from measurements
    virtual void            FindGreen(TOKAMAK *td) = 0;
    virtual void            FindFit(TOKAMAK *td) = 0;
    virtual void            FindNow(TOKAMAK *td) = 0;
    virtual void            FindL(TOKAMAK *td, double *L) = 0;

    // read and create measurement from file!
    static CMeasurement *   CreateWithType(char *);
    static CMeasurement *   CreateWithType(int);
    virtual void            SetParameter(char *, char *);

    static void             FindAllGreens(CMeasurementList &mlist);
    static void             FindAllFit(CMeasurementList &mlist);
    static void             FindAllNow(CMeasurementList &mlist);
    static void             FindAllL(CMeasurementList &mlist);
    static void             FreeAll(CMeasurementList &mlist);

    // destructor
    virtual ~CMeasurement(void);
    // virtual constructors!
    virtual CMeasurement * create() const = 0;
    virtual CMeasurement * clone() const = 0;
};

#endif // _CMeasurement_h_
