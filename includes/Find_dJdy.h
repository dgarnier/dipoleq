/*
** TokaMac v2.0
**
** Find_dJdy.h
**
**
**
** File:		Find_dJdy.h
** Date:		March 22, 1993
**
** Routine list:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/


#ifndef _FIND_DJDY_

#define _FIND_DJDY_  1

#ifdef __cplusplus
extern "C" {
#endif

double ***new_dJdy(int , int );
void 	free_dJdy(double ***,int , int);

void Find_dJdy(TOKAMAK *, double ***);

#ifdef __cplusplus
}
#endif

#endif
