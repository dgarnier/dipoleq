/*
** TokaMac v2.0
**
** LoadBndryGreens.h
**
** Interface file for LoadBndryGreens.c
**
** File:		LoadBndryGreens.h
** Date:		March 21, 1993
**
**	USAGE:
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _LOADBNDRYGREENS_

#define _LOADBNDRYGREENS_ 		1

void LoadBndryGreens(TOKAMAK *);
void RewriteBndryGreens(TOKAMAK *);

void free_BndryGreens(TOKAMAK *);

#endif