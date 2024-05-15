/*
** TokaMac v2.0
**
** Interface file for tokgreen.c.
**
** File:		include:tokgreen.h
** Date:		February 28, 1993
**
** Example:
**
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _TOKGREEN_

#define _TOKGREEN_ 1

/*
**
**	LHVEC
**
**	This stores the LH Green's functions for a coil set
**	or some current system with a fixed current distribution.
**
*/

  typedef struct lhvec {
      double *Top;
      double *Bot;
      double *In;
      double *Out;
  } LHVEC;

/*
**
**	LHARY
**
**	This stores an array of LHVECs for the Lackner & Von Hagenow
**	procedure.
**
*/

  typedef struct lhary {
      LHVEC **In;
      LHVEC **Out;
      LHVEC **Bot;
  } LHARY;

/*
**
**	Public subroutines
**
*/

#ifdef __cplusplus
extern "C" {
#endif

LHVEC *new_LHvec(int);
void free_LHvec(LHVEC *, int);

LHARY *new_LHary(int);
void free_LHary(LHARY *, int);

#ifdef __cplusplus
}  // extern "C"
#endif
#endif
