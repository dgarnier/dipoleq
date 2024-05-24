/*
** TokaMac v2.0
**
** FileInput.c
**
** This file defines subroutines to read and parse
** ASCII files for program input and control.
**
** File:		FileInput.c
** Date:		March 7, 1993
**
** Revisions:
**
**		August 5, 1993		Added perfectly conducting shells
**		Oct. 3, 1993		Added RestartUnkns
**		November 15, 1993	Angle MSE measurements are input in degrees
**
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#include <stdlib.h>
#include "VAX_Alloc.h"
#include <stdio.h>
#include <string.h>
#include "nrutil.h"
#include "coil.h"
#include "shell.h"
#include "limiter.h"
#include "separatrix.h"
#include "measurement.h"
#include "plasma.h"
#include "psigrid.h"
#include "tokamak.h"
#include "FileInput.h"
#include "CPlasmaModel.h"

int           getline(char *line, int max, FILE * fi);
int           getkey(char *kw);
int           getControlVar(char *w);
int           getPsiGridVar(char *w);
int           getPlasmaVar(char *w);
int           getLimiterVar(char *w);
int           getSeparatrixVar(char *w);
int           getCoilVar(char *w);
int           getSubCoilVar(char *w);
int           getShellVar(char *w);
int           getSubShellVar(char *w);
int           getMeasureVar(char *w);

void          AssignControlVar(TOKAMAK * td, char *word, char *value);
void          AssignPsiGridVar(TOKAMAK * td, char *word, char *value);
void          AssignPlasmaVar(TOKAMAK * td, char *word, char *value);
void          AssignLimiterVar(TOKAMAK * td, int c, int isNew, char *word, char *value);
void          AssignSeparatrixVar(TOKAMAK * td, int c, int isNew, char *word, char *value);
void          AssignCoilVar(COIL * c, char *word, char *value);
void          AssignSubCoilVar(SUBCOIL * sc, char *word, char *value);
void          AssignShellVar(SHELL * sh, char *word, char *value);
void          AssignSubShellVar(SUBSHELL * sh, char *word, char *value);
void          AssignControlVar(TOKAMAK * td, char *word, char *value);
void          AssignPsiGridVar(TOKAMAK * td, char *word, char *value);
void          AssignMeasureVar(TOKAMAK * td, int c, int isNew, char *word, char *value);




#define PI          3.14159265358979323
#define MU0			1.25663706e-06
#define TWOPI		6.283185307

extern "C" FILE  *LogFile;



/*
**
**	getline
**
*/

int           getline(char *line, int max, FILE * fi)
{
	if (fgets(line, max, fi) == NULL)
		return 0;
	else if ((line[0] == '/' && line[1] == '/') || (line[0] == '\\' && line[1] == '\\'))
		return 0;				/* ignore C++ comment lines */
	else
		return (int) strlen(line);
}

/*
**
**	getkey, getControlVar, getPsiGridVar, etc.
**
*/

int           getkey(char *kw)
{
	int           i;
	for (i = K_CodeControl; i <= K_End; i++)
		if (!strcmp(KeyWords[i], kw))
			return i;
	return K_NoKey;
}

int           getControlVar(char *w)
{
	int           i;
	for (i = C_MaxIterFixed; i <= C_MGreenStatus; i++)
		if (!strcmp(ControlWords[i], w))
			return i;
	return NoWord;
}

int           getPsiGridVar(char *w)
{
	int           i;
	for (i = Ps_Nsize; i <= Ps_UnderRelax1; i++)
		if (!strcmp(PsiGridWords[i], w))
			return i;
	return NoWord;
}

int           getPlasmaVar(char *w)
{
	int           i;
	for (i = Pl_ModelType; i <= Pl_B0; i++)
		if (!strcmp(PlasmaWords[i], w))
			return i;
	return NoWord;
}

int           getLimiterVar(char *w)
{
	int           i;
	for (i = LS_X1; i <= LS_Z2; i++)
		if (!strcmp(LimiterWords[i], w))
			return i;
	return NoWord;
}

int           getSeparatrixVar(char *w)
{
	int           i;
	for (i = LS_X1; i <= LS_ZC; i++)
		if (!strcmp(SeparatrixWords[i], w))
			return i;
	return NoWord;
}

int           getCoilVar(char *w)
{
	int           i;
	for (i = Coil_InitialCurrent; i <= Coil_Enabled; i++)
		if (!strcmp(CoilWords[i], w))
			return i;
	return NoWord;
}

int           getSubCoilVar(char *w)
{
	int           i;
	for (i = SC_X; i <= SC_CurrentFraction; i++)
		if (!strcmp(SubCoilWords[i], w))
			return i;
	return NoWord;
}

int           getShellVar(char *w)
{
	int           i;
	for (i = Shell_Name; i <= Shell_Enabled; i++)
		if (!strcmp(ShellWords[i], w))
			return i;
	return NoWord;
}

int           getSubShellVar(char *w)
{
	int           i;
	for (i = SS_X; i <= SS_Radius; i++)
		if (!strcmp(SubShellWords[i], w))
			return i;
	return NoWord;
}

int           getMeasureVar(char *w)
{
	int           i;
	for (i = M_X1; i <= M_CoilNum; i++)
		if (!strcmp(MeasureWords[i], w))
			return i;
	return NoWord;
}

/*
**
**	AssignControlVar
**
*/

void          AssignControlVar(TOKAMAK * td, char *word, char *value)
{
	int           i;

	if ((i = getControlVar(word)) != 0)
		switch (i) {
		  case C_MaxIterFixed:
			  sscanf(value, "%d", &(td->MaxIterFixed));
			  break;
		  case C_MaxIterFree:
			  sscanf(value, "%d", &(td->MaxIterFree));
			  break;
		  case C_Name:
			  sscanf(value, "%s", (td->Name));
			  break;
		  case C_Info:
			  sscanf(value, "%s", (td->Info));
			  break;
		  case C_NumCoils:
			  sscanf(value, "%d", &(td->NumCoils));
			  break;
		  case C_NumShells:
			  sscanf(value, "%d", &(td->NumShells));
			  break;
		  case C_NumLimiters:
			  sscanf(value, "%d", &(td->NumLimiters));
			  break;
		  case C_NumSeps:
			  sscanf(value, "%d", &(td->NumSeps));
			  break;
		  case C_NumMeasures:
			  sscanf(value, "%d", &(td->NumMeasures));
			  break;
		  case C_Oname:
			  sscanf(value, "%s", (td->Oname));
			  break;
		  case C_LHname:
			  sscanf(value, "%s", (td->LHname));
			  break;
		  case C_MGname:
			  sscanf(value, "%s", (td->MGname));
			  break;
		  case C_SGname:
			  sscanf(value, "%s", (td->SGname));
			  break;
		  case C_SMname:
			  sscanf(value, "%s", (td->SMname));
			  break;
		  case C_RSname:
			  sscanf(value, "%s", (td->RSname));
			  break;
		  case C_RestartStatus:
			  sscanf(value, "%d", &(td->RestartStatus));
			  break;
		  case C_RestartUnkns:
			  sscanf(value, "%d", &(td->RestartUnkns));
			  break;
		  case C_LHGreenStatus:
			  sscanf(value, "%d", &(td->LHGreenStatus));
			  break;
		  case C_SGreenStatus:
			  sscanf(value, "%d", &(td->SGreenStatus));
			  break;
		  case C_SInductStatus:
			  sscanf(value, "%d", &(td->SInductStatus));
			  break;
		  case C_NumEqualEq:
			  sscanf(value, "%d", &(td->NumEqualEq));
			  break;
		  case C_Confidence:
			  sscanf(value, "%lf", &(td->Confidence));
			  break;
		  case C_VacuumOnly:
			  sscanf(value, "%d", &(td->VacuumOnly));
			  break;
		  case C_NumMCarloEq:
			  sscanf(value, "%d", &(td->NumMCarloEq));
			  break;
		  case C_NumMCarloData:
			  sscanf(value, "%d", &(td->NumMCarloData));
			  break;
		  case C_MaxIterMCarlo:
			  sscanf(value, "%d", &(td->MaxIterMCarlo));
			  break;
		  case C_MGreenStatus:
			  sscanf(value, "%d", &(td->MGreenStatus));
			  break;
		}
}

/*
**
**	AssignPsiGridVar
**
*/

void          AssignPsiGridVar(TOKAMAK * td, char *word, char *value)
{
	int           i;
	PSIGRID      *p;

	p = td->PsiGrid;

	if ((i = getPsiGridVar(word))!=0)
		switch (i) {
		  case Ps_Nsize:
			  sscanf(value, "%d", &(p->Nsize));
			  break;
		  case Ps_Symmetric:
			  sscanf(value, "%d", &(p->Symmetric));
			  break;
		  case Ps_Xmax:
			  sscanf(value, "%lf", &(p->Xmax));
			  break;
		  case Ps_Xmin:
			  sscanf(value, "%lf", &(p->Xmin));
			  break;
		  case Ps_Zmax:
			  sscanf(value, "%lf", &(p->Zmax));
			  break;
		  case Ps_Zmin:
			  sscanf(value, "%lf", &(p->Zmin));
			  break;
		  case Ps_BoundThreshold:
			  sscanf(value, "%lf", &(p->BoundThreshold));
			  break;
		  case Ps_ResThreshold:
			  sscanf(value, "%lf", &(p->ResThreshold));
			  break;
		  case Ps_UnderRelax1:
			  sscanf(value, "%lf", &(p->UnderRelax1));
			  break;
		  case Ps_UnderRelax2:
			  sscanf(value, "%lf", &(p->UnderRelax2));
			  break;
		}
}

/*
**
**	AssignPlasmaVar
**
*/

void          AssignPlasmaVar(TOKAMAK * td, char *word, char *value)
{
	int           i;
	PSIGRID      *pg;
	PLASMA       *p;

	p = td->Plasma;
	pg = td->PsiGrid;

	p->Nsize = pg->Nsize;		/* copy grid size from PsiGrid */

	if ((i = getPlasmaVar(word))!=0)
		switch (i) {
		  case Pl_ModelType:
			  sscanf(value, "%d", &(p->ModelType));
			  break;
		  case Pl_PpTerms:
			  sscanf(value, "%d", &(p->PpTerms));
			  if (p->PpTerms > MaxPolyTerms)
				  p->PpTerms = MaxPolyTerms;
			  if (p->PpTerms < 1)
				  p->PpTerms = 1;
			  break;
		  case Pl_G2pTerms:
			  sscanf(value, "%d", &(p->G2pTerms));
			  break;
		  case Pl_SisoTerms:
			  sscanf(value, "%d", &(p->SisoTerms));
			  break;
		  case Pl_SperTerms:
			  sscanf(value, "%d", &(p->SperTerms));
			  break;
		  case Pl_SparTerms:
			  sscanf(value, "%d", &(p->SparTerms));
			  break;
		  case Pl_RotTerms:
			  sscanf(value, "%d", &(p->RotTerms));
			  break;
		  case Pl_HTerms:
			  sscanf(value, "%d", &(p->HTerms));
			  break;
		  case Pl_StndP:
			  sscanf(value, "%lf", &(p->StndP));
			  if (p->StndP == 1.0000)
				  p->StndP = 1.00001;
			  break;
		  case Pl_StndG:
			  sscanf(value, "%lf", &(p->StndG));
			  if (p->StndG == 1.0000)
				  p->StndG = 1.00001;
			  break;
		  case Pl_NumBndMomts:
			  sscanf(value, "%d", &(p->NumBndMomts));
			  if (p->NumBndMomts < 2)
				  p->NumBndMomts = 2;
			  break;
		  case Pl_NumPsiPts:
			  sscanf(value, "%d", &(p->NumPsiPts));
			  if (p->NumPsiPts < 4)
				  p->NumPsiPts = 4;
			  break;
		  case Pl_PsiXmax:
			  sscanf(value, "%lf", &(p->PsiXmax));
			  if (p->PsiXmax > 0.9999999)
				  p->PsiXmax = 0.9999999;
			  break;
		  case Pl_R0:
			  sscanf(value, "%lf", &(p->R0));
			  break;
		  case Pl_a0:
			  sscanf(value, "%lf", &(p->a0));
			  break;
		  case Pl_Z0:
			  sscanf(value, "%lf", &(p->Z0));
			  break;
		  case Pl_Ip0:
			  sscanf(value, "%lf", &(p->Ip0));
			  break;
		  case Pl_Jedge:
			  sscanf(value, "%lf", &(p->Jedge));
			  p->Jedge *= MU0;		/* stored as MU0*J */
			  break;
		  case Pl_B0:
			  sscanf(value, "%lf", &(p->B0));
			  break;
		}
	p->B0R0 = p->R0 * p->B0;
}


/*
**
**	AssignLimiterVar
**
*/

void          AssignLimiterVar(TOKAMAK * td, int c, int isNew, char *word, char *value)
{
	int           i;
	LIMITER      *lm;

	if (isNew)
		td->Limiters[c] = new_Limiter();
	lm = td->Limiters[c];

	if ((i = getLimiterVar(word))!=0)
		switch (i) {
		  case LS_X1:
			  sscanf(value, "%lf", &(lm->X1));
			  break;
		  case LS_Z1:
			  sscanf(value, "%lf", &(lm->Z1));
			  break;
		  case LS_X2:
			  sscanf(value, "%lf", &(lm->X2));
			  break;
		  case LS_Z2:
			  sscanf(value, "%lf", &(lm->Z2));
			  break;
		  case LS_Enabled:
			  sscanf(value, "%d", &(lm->Enabled));
			  break;
		  case LS_Name:
			  sscanf(value, "%s", (lm->Name));
			  break;
		}
}

/*
**
**	AssignSeparatrixVar
**
*/

void          AssignSeparatrixVar(TOKAMAK * td, int c, int isNew, char *word, char *value)
{
	int           i;
	SEPARATRIX   *sp;

	if (isNew)
		td->Seps[c] = new_Separatrix();
	sp = td->Seps[c];

	if ((i = getSeparatrixVar(word))!=0)
		switch (i) {
		  case LS_X1:
			  sscanf(value, "%lf", &(sp->X1));
			  break;
		  case LS_Z1:
			  sscanf(value, "%lf", &(sp->Z1));
			  break;
		  case LS_X2:
			  sscanf(value, "%lf", &(sp->X2));
			  break;
		  case LS_Z2:
			  sscanf(value, "%lf", &(sp->Z2));
			  break;
		  case LS_XC:
			  sscanf(value, "%lf", &(sp->XC));
			  break;
		  case LS_ZC:
			  sscanf(value, "%lf", &(sp->ZC));
			  break;
		  case LS_Enabled:
			  sscanf(value, "%d", &(sp->Enabled));
			  break;
		  case LS_Name:
			  sscanf(value, "%s", (sp->Name));
			  break;
		}
}

/*
**
**	AssignMeasureVar
**
*/

void          AssignMeasureVar(TOKAMAK * td, int c, int isNew, char *word, char *value)
{
	int           i, mt;
	MEAS         *m;

	/* The first value of a Measure input record must be the mType */

	if ((i = getMeasureVar(word))!=0) {
		if (isNew) {
			if (i != M_mType) nrerror("Input file error!  mType not first in Measure input record.\n");
			sscanf(value, "%d", &mt);	/* read mType */
			m = new_Measure(mt);
			m->mType = mt;
			td->Measures[c] = m;
		} else {
			m = td->Measures[c];
			switch (i) {
			  case M_Value:
				  sscanf(value, "%lf", &(m->Value));
				  break;
			  case M_StdDev:
				  sscanf(value, "%lf", &(m->StdDev));
				  break;
			  case M_Name:
				  sscanf(value, "%s", (m->Name));
				  break;
			  case M_mX:
				  sscanf(value, "%lf", &(m->X));
				  break;
			  case M_mZ:
				  sscanf(value, "%lf", &(m->Z));
				  break;
			}
			switch (m->mType) {
			  case meas_bp:
				  switch (i) {
					case M_Angle:
						sscanf(value, "%lf", &(m->parm.bp.Angle));
						break;
				  }
				  break;
			  case meas_saddle:
				  switch (i) {
					case M_X1:
						sscanf(value, "%lf", &(m->parm.saddle.X1));
						break;
					case M_Z1:
						sscanf(value, "%lf", &(m->parm.saddle.Z1));
						break;
					case M_X2:
						sscanf(value, "%lf", &(m->parm.saddle.X2));
						break;
					case M_Z2:
						sscanf(value, "%lf", &(m->parm.saddle.Z2));
						break;
				  }
				  break;
			  case meas_circle:
				  switch (i) {
					case M_Radius:
						sscanf(value, "%lf", &(m->parm.circle.Radius));
						break;
					case M_Number:
						sscanf(value, "%d", &(m->parm.circle.Number));
						break;
					case M_CircleType:
						sscanf(value, "%d", &(m->parm.circle.CircleType));
						break;
				  }
				  break;
			  case meas_coilcur:
				  switch (i) {
					case M_CoilNum:
						sscanf(value, "%d", &(m->parm.coilcur.CoilNum));
						m->parm.coilcur.CoilNum -= 1;
						break;
				  }
				  break;
			  case meas_bpangle:
				  switch (i) {
					case M_Value:
						m->Value = m->Value*(PI/180.0); /* convert to radians */
						break;
					case M_StdDev:
						m->StdDev = m->StdDev*(PI/180.0); /* convert to radians */
						break;
				  }
				  break;
			}
		}
	}
}

/*
**
**	AssignCoilVar
**
*/

void          AssignCoilVar(COIL * c, char *word, char *value)
{
	int           i;

	if ((i = getCoilVar(word))!=0)
		switch (i) {
		  case Coil_InitialCurrent:
			  sscanf(value, "%lf", &(c->CoilCurrent));
			  c->CoilCurrent *= MU0;
			  break;
		  case Coil_Enabled:
			  sscanf(value, "%d", &(c->Enabled));
			  break;
		  case Coil_Name:
			  sscanf(value, "%s", (c->Name));
			  break;
		  case Coil_X:
			  sscanf(value, "%lf", &(c->X));
			  break;
		  case Coil_DX:
			  sscanf(value, "%lf", &(c->dX));
			  break;
		  case Coil_Z:
			  sscanf(value, "%lf", &(c->Z));
			  break;
		  case Coil_DZ:
			  sscanf(value, "%lf", &(c->dZ));
			  break;

		}
}

/*
**
**	AssignSubCoilVar
**
*/

void          AssignSubCoilVar(SUBCOIL * sc, char *word, char *value)
{
	int           i;

	if ((i = getSubCoilVar(word))!=0)
		switch (i) {
		  case SC_X:
			  sscanf(value, "%lf", &(sc->X));
			  break;
		  case SC_Z:
			  sscanf(value, "%lf", &(sc->Z));
			  break;
		  case SC_CurrentFraction:
			  sscanf(value, "%lf", &(sc->CurrentFraction));
			  break;
		  case SC_Name:
			  sscanf(value, "%s", (sc->Name));
			  break;
		}
}

/*
**
**	AssignShellVar
**
*/

void          AssignShellVar(SHELL * c, char *word, char *value)
{
	int           i;

	if ((i = getShellVar(word))!=0)
		switch (i) {
		  case Shell_Enabled:
			  sscanf(value, "%d", &(c->Enabled));
			  break;
		  case Shell_Name:
			  sscanf(value, "%s", (c->Name));
			  break;
		}
}

/*
**
**	AssignSubShellVar
**
*/

void          AssignSubShellVar(SUBSHELL * sc, char *word, char *value)
{
	int           i;

	if ((i = getSubShellVar(word))!=0)
		switch (i) {
		  case SS_X:
			  sscanf(value, "%lf", &(sc->X));
			  break;
		  case SS_Z:
			  sscanf(value, "%lf", &(sc->Z));
			  break;
		  case SS_Radius:
			  sscanf(value, "%lf", &(sc->Radius));
			  break;
		  case SS_Name:
			  sscanf(value, "%s", (sc->Name));
			  break;
		}
}

/*
**
**	FileInput
**
*/

TOKAMAK      *FileInput(const char *fname)
{
	FILE         *fi;
	TOKAMAK      *td;
	PSIGRID      *pg;
	PLASMA       *pl;
	SUBCOIL      *sc;
	SUBSHELL     *ss;
	int           i = 0;		/* line number */
	int           linelen, count, isNew;
	int           isOK;
	int           theKey;		/* the active input group */
	char          linebuf[128], keystr[24], valstr[24], *p;

	td = new_Tokamak();
	pl = td->Plasma;
	pg = td->PsiGrid;

	strncpy(td->Iname, fname, sizeof(td->Iname) - 1);
	fi = fopen(fname, "r");		/* open textfile for read only */
	if (!fi)
		nrerror("ERROR:	Could not open file in FileInput.");

	/* C O D E C O N T R O L */
	i++;
	linelen = getline(linebuf, 128, fi);
        isOK = FALSE;
	while (!feof(fi)) {
		if (linelen && sscanf(linebuf, "K_%s", keystr))
			isOK = (getkey(keystr) == K_CodeControl);
		else if (linelen && isOK && (sscanf(linebuf, "%s = %s", keystr, valstr) == 2))
			AssignControlVar(td, keystr, valstr);
		i++;
		linelen = getline(linebuf, 128, fi);
	}




	/* P S I G R I D */
	rewind(fi);
	i = 1;
	linelen = getline(linebuf, 128, fi);
        isOK = FALSE;
	while (!feof(fi)) {
		if (linelen && sscanf(linebuf, "K_%s", keystr))
			isOK = (getkey(keystr) == K_PsiGrid);
		else if (linelen && isOK && (sscanf(linebuf, "%s = %s", keystr, valstr) == 2))
			AssignPsiGridVar(td, keystr, valstr);
		i++;
		linelen = getline(linebuf, 128, fi);
	}

	/* P L A S M A */
	rewind(fi);
	i = 1;
	linelen = getline(linebuf, 128, fi);
        isOK = FALSE;
	while (!feof(fi)) {
		if (linelen && sscanf(linebuf, "K_%s", keystr))
			isOK = (getkey(keystr) == K_Plasma);
		else if (linelen && isOK && (sscanf(linebuf, "%s = %s", keystr, valstr) == 2))
			AssignPlasmaVar(td, keystr, valstr);
		i++;
		linelen = getline(linebuf, 128, fi);
	}

	/* Allocate arrays for PsiGrid and Plasma */
	init_Tokamak(td);

	/* P L A S M A model initialization */
	/* this loads individual model definitions */
	if (pl->Model) {
		rewind(fi);
		i = 1;
		linelen = getline(linebuf, 128, fi);
                isOK = FALSE;
		while (!feof(fi)) {
			if (linelen && sscanf(linebuf, "K_%s", keystr)) {
				isOK = (getkey(keystr) == K_Plasma);
				linelen = getline(linebuf, 128, fi);
			} else if (linelen && isOK && (sscanf(linebuf, "%s = %s", keystr, valstr) == 2)) {
				linelen = getline(linebuf, 128, fi);
				if (linelen && (sscanf(linebuf, "K_%s", keystr) == 0) &&
						(strstr(linebuf,"Init ->"))){
					p = strstr(linebuf,"->") + 2;
					if (!feof(fi)) linelen = getline(linebuf, 128, fi);
				} else p = NULL;
				pl->Model->ModelInput(keystr,valstr,p);
			} else
				linelen = getline(linebuf, 128, fi);
		}
	}

	/* L I M I T E R S */
	count = -1;
	rewind(fi);
	i = 1;
	linelen = getline(linebuf, 128, fi);
        isOK = FALSE;
	while (!feof(fi)) {
		if (linelen && sscanf(linebuf, "K_%s", keystr)) {
			isOK = (getkey(keystr) == K_Limiter);
			if (isOK) {
				count++;
				isNew = TRUE;
			}
		} else if (linelen && isOK && (sscanf(linebuf, "%s = %s", keystr, valstr) == 2)) {
			if (count < td->NumLimiters) {
				AssignLimiterVar(td, count, isNew, keystr, valstr);
				isNew = FALSE;
			}
		}
		i++;
		linelen = getline(linebuf, 128, fi);
	}
	if (count < td->NumLimiters - 1) {	/* not enough Limiters */
		printf("\nWARNING from FileInput.  Not enough Limiters in file.\n");
		fprintf(LogFile, "\nWARNING from FileInput.  Not enough Limiters in file.\n");
	}
	/* S E P A R A T R I C I E S */
	count = -1;
	rewind(fi);
	i = 1;
	linelen = getline(linebuf, 128, fi);
        isOK = FALSE;
	while (!feof(fi)) {
		if (linelen && sscanf(linebuf, "K_%s", keystr)) {
			isOK = (getkey(keystr) == K_Separatrix);
			if (isOK) {
				count++;
				isNew = TRUE;
			}
		} else if (linelen && isOK && (sscanf(linebuf, "%s = %s", keystr, valstr) == 2)) {
			if (count < td->NumSeps) {
				AssignSeparatrixVar(td, count, isNew, keystr, valstr);
				isNew = FALSE;
			}
		}
		i++;
		linelen = getline(linebuf, 128, fi);
	}
	if (count < td->NumSeps - 1) {	/* not enough Limiters */
		printf("\nWARNING from FileInput.  Not enough Separarticies in file.\n");
		fprintf(LogFile, "\nWARNING from FileInput.  Not enough Separarticies in file.\n");
	}
	/* M E A S U R E M E N T S */
	count = -1;
        isOK = FALSE;
	rewind(fi);
	i = 1;
	linelen = getline(linebuf, 128, fi);
	while (!feof(fi)) {
		if (linelen && sscanf(linebuf, "K_%s", keystr)) {
			isOK = (getkey(keystr) == K_Measure);
			if (isOK) {
				count++;
				isNew = TRUE;
			}
		} else if (linelen && isOK && (sscanf(linebuf, "%s = %s", keystr, valstr) == 2)) {
			if (count < td->NumMeasures) {
				AssignMeasureVar(td, count, isNew, keystr, valstr);
				isNew = FALSE;
			}
		}
		i++;
		linelen = getline(linebuf, 128, fi);
	}
	if (count < td->NumMeasures - 1) {	/* not enough Limiters */
		printf("\nWARNING from FileInput.  Not enough Measurements in file.\n");
		fprintf(LogFile, "\nWARNING from FileInput.  Not enough Measurements in file.\n");
	}
	/* C O I L S and S U B C O I L S */
	count = -1;
        isOK = FALSE;
	rewind(fi);
	i = 1;
	linelen = getline(linebuf, 128, fi);
	while (!feof(fi)) {
		if (linelen && sscanf(linebuf, "K_%s", keystr)) {
			theKey = getkey(keystr);
			isOK = (theKey == K_Coil) || (theKey == K_SubCoil);
			if (theKey == K_Coil)
				count++;
			if (theKey == K_SubCoil) {
				sc = new_SubCoil();
				add_SubCoil(td->Coils[count], sc);
			}
		} else if (linelen && isOK && (sscanf(linebuf, "%s = %s", keystr, valstr) == 2)) {
			if ((theKey == K_Coil) && (count < td->NumCoils))
				AssignCoilVar(td->Coils[count], keystr, valstr);
			if (theKey == K_SubCoil)
				AssignSubCoilVar(sc, keystr, valstr);
		}
		i++;
		linelen = getline(linebuf, 128, fi);
	}
	if (count < td->NumCoils - 1) {
		printf("\nWARNING from FileInput.  Not enough Coils in file.\n");
		fprintf(LogFile, "\nWARNING from FileInput.  Not enough Coils in file.\n");
	}
	for (i=0;i<=count;i++)
		if (td->Coils[i]->dX >= 0) compute_SubCoils(td->Coils[i], td->PsiGrid);
	/* S H E L L S and S U B S H E L L S */
	count = -1;
        isOK = FALSE;
	rewind(fi);
	i = 1;
	linelen = getline(linebuf, 128, fi);
	while (!feof(fi)) {
		if (linelen && sscanf(linebuf, "K_%s", keystr)) {
			theKey = getkey(keystr);
			isOK = (theKey == K_Shell) || (theKey == K_SubShell);
			if (theKey == K_Shell)
				count++;
			if (theKey == K_SubShell) {
				ss = new_SubShell();
				add_SubShell(td->Shells[count], ss);
			}
		} else if (linelen && isOK && (sscanf(linebuf, "%s = %s", keystr, valstr) == 2)) {
			if ((theKey == K_Shell) && (count < td->NumShells))
				AssignShellVar(td->Shells[count], keystr, valstr);
			if (theKey == K_SubShell)
				AssignSubShellVar(ss, keystr, valstr);
		}
		i++;
		linelen = getline(linebuf, 128, fi);
	}
	if (count < td->NumShells - 1) {
		printf("\nWARNING from FileInput.  Not enough Shells in file.\n");
		fprintf(LogFile, "\nWARNING from FileInput.  Not enough Shells in file.\n");
	}
	fclose(fi);
	return td;
}
