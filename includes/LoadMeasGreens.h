/*
** TokaMac v2.0
**
** LoadMeasGreens.h
**
** Reads the Meas Green's functions from disk.
** Also, it rewrites the file (erasing if it already exists),
** by calling RewriteMeasGreens.
**
** File:		LoadMeasGreens.h
** Date:		March 24, 1993
**
**	USAGE:
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _LOADMEASGREENS_

#define _LOADMEASGREENS_		1

#ifdef __cplusplus
extern "C" {
#endif

void          LoadMeasGreens(TOKAMAK * );
void          RewriteMeasGreens(TOKAMAK * );
void          free_MeasGreens(TOKAMAK * );

#ifdef __cplusplus
}
#endif

#endif
