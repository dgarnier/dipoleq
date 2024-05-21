/*
** TokaMac v2.0
**
** FileInput.h
**
** This file defines subroutines to read and parse
** ASCII files for program input and control.
**
** File:		FileInput.h
** Date:		March 7, 1993
**
** Revisions:
**
**		August 5, 1993		Added perfectly conducting shells
**		Oct. 3, 1993		Added RestartUnkns
**
** (c) L. Bai and M. Mauel -- Columbia University
*/

#ifndef _FILEINPUT_

#define _FILEINPUT_ 1

#ifdef __cplusplus
extern "C" {
#endif

#define TRUE	1
#define FALSE	0

/*
**
** 	Major Input KeyWord Headings
**
*/

static const char *KeyWords[] =
{
    "NoKey",
    "CodeControl",
    "PsiGrid",
    "Plasma",
    "Measure",
    "Limiter",
    "Separatrix",
    "Shell",
    "SubShell",
    "Coil",
    "SubCoil",
    "End"
};

#define K_NoKey			0
#define K_CodeControl	1
#define K_PsiGrid		2
#define K_Plasma		3
#define K_Measure		4
#define K_Limiter		5
#define K_Separatrix	6
#define K_Shell			7
#define K_SubShell		8
#define K_Coil			9
#define K_SubCoil		10
#define K_End			11


/*
**
**	Control input variables
**
*/

static const char *ControlWords[] =
{
    "NoWord",
    "MaxIterFixed",
    "MaxIterFree",
    "Name",
    "Info",
    "NumCoils",
    "NumShells",
    "NumLimiters",
    "NumSeps",
    "NumMeasures",
    "Oname",
    "LHname",
    "MGname",
    "SGname",
    "SMname",
    "RSname",
    "RestartStatus",
    "RestartUnkns",
    "LHGreenStatus",
    "SGreenStatus",
    "SInductStatus",
    "NumEqualEq",
    "Confidence",
    "VacuumOnly",
    "NumMCarloEq",
    "NumMCarloData",
    "MaxIterMCarlo",
    "MGreenStatus"

};

#define NoWord				0
#define C_MaxIterFixed		1
#define C_MaxIterFree		2
#define C_Name				3
#define C_Info				4
#define C_NumCoils			5
#define C_NumShells			6
#define C_NumLimiters		7
#define C_NumSeps			8
#define C_NumMeasures		9
#define C_Oname				10
#define C_LHname			11
#define C_MGname			12
#define C_SGname			13
#define C_SMname			14
#define C_RSname			15
#define C_RestartStatus		16
#define C_RestartUnkns		17
#define C_LHGreenStatus		18
#define C_SGreenStatus		19
#define C_SInductStatus		20
#define C_NumEqualEq		21
#define C_Confidence		22
#define C_VacuumOnly        23
#define C_NumMCarloEq		24
#define C_NumMCarloData		25
#define C_MaxIterMCarlo		26
#define C_MGreenStatus		27

/*
**
**	PsiGrid input variables
**
*/

static const char *PsiGridWords[] =
{
    "NoWord",
    "Nsize",
    "Symmetric",
    "Xmax",
    "Xmin",
    "Zmax",
    "Zmin",
    "BoundThreshold",
    "ResThreshold",
    "UnderRelax2",
    "UnderRelax1"
};

#define Ps_Nsize 	 		1
#define Ps_Symmetric		2
#define Ps_Xmax				3
#define Ps_Xmin				4
#define Ps_Zmax				5
#define Ps_Zmin				6
#define Ps_BoundThreshold	7
#define Ps_ResThreshold		8
#define Ps_UnderRelax2		9
#define Ps_UnderRelax1		10


/*
**
**	Plasma input variables
**
*/

static const char *PlasmaWords[] =
{
    "NoWord",
    "ModelType",
    "PpTerms",
    "G2pTerms",
    "SisoTerms",
    "SperTerms",
    "SparTerms",
    "RotTerms",
    "HTerms",
    "StndP",
    "StndG",
    "NumBndMomts",
    "NumPsiPts",
    "PsiXmax",
    "R0",
    "Z0",
    "a0",
    "Ip0",
    "Jedge",
    "B0",
};

#define Pl_ModelType 	 	1
#define Pl_PpTerms			2
#define Pl_G2pTerms			3
#define Pl_SisoTerms		4
#define Pl_SperTerms		5
#define Pl_SparTerms		6
#define Pl_RotTerms			7
#define Pl_HTerms			8
#define Pl_StndP			9
#define Pl_StndG			10
#define Pl_NumBndMomts		11
#define Pl_NumPsiPts		12
#define Pl_PsiXmax			13
#define Pl_R0				14
#define Pl_Z0				15
#define Pl_a0				16
#define Pl_Ip0				17
#define Pl_Jedge			18
#define Pl_B0				19


/*
**
**	Limiter input variables
**
*/

static const char *LimiterWords[] =
{
    "NoWord",
    "X1",
    "Z1",
    "Name",
    "Enabled",
    "X2",
    "Z2"
};

#define LS_X1 	 	1
#define LS_Z1 	 	2
#define LS_Name 	3
#define LS_Enabled 	4
#define LS_X2 	 	5
#define LS_Z2 	 	6
#define LS_XC 	 	7
#define LS_ZC 	 	8


/*
**
**	Separatrix input variables
**
*/

static const char *SeparatrixWords[] =
{
    "NoWord",
    "X1",
    "Z1",
    "Name",
    "Enabled",
    "X2",
    "Z2",
    "XC",
    "ZC"
};


/*
**
**	Coil input variables
**
*/

static const char *CoilWords[] =
{
    "NoWord",
    "NoWord",
    "InitialCurrent",
    "Name",
    "X","dX","Z","dZ",
    "Enabled"
};

#define Coil_InitialCurrent 	2
#define Coil_Name			 	3
#define Coil_X					4
#define Coil_DX					5
#define Coil_Z					6
#define Coil_DZ					7
#define Coil_Enabled		 	8



/*
**
**	SubCoil input variables
**
*/

static const char *SubCoilWords[] =
{
    "NoWord",
    "X",
    "Z",
    "Name",
    "CurrentFraction"
};

#define SC_X 				1
#define SC_Z 				2
#define SC_Name 			3
#define SC_CurrentFraction	4


/*
**
**	Shell input variables
**
*/

static const char *ShellWords[] =
{
    "NoWord",
    "NoWord",
    "NoWord",
    "Name",
    "Enabled"
};

#define Shell_Name			 	3
#define Shell_Enabled		 	4


/*
**
**	SubShell input variables
**
*/

static const char *SubShellWords[] =
{
    "NoWord",
    "X",
    "Z",
    "Name",
    "Radius"
};

#define SS_X 				1
#define SS_Z 				2
#define SS_Name 			3
#define SS_Radius			4


/*
**
**	Measurement input variables
**
*/

static const char *MeasureWords[] =
{
    "NoWord",
    "X1",
    "Z1",
    "Name",
    "mType",
    "X2",
    "Z2",
    "X",
    "Z",
    "Radius",
    "Number",
    "CircleType",
    "Angle",
    "Value",
    "StdDev",
    "CoilNum"
};

#define M_X1 	 		1
#define M_Z1 	 		2
#define M_Name 			3
#define M_mType			4
#define M_X2			5
#define M_Z2			6
#define M_mX			7
#define M_mZ			8
#define M_Radius		9
#define M_Number		10
#define M_CircleType	11
#define M_Angle			12
#define M_Value			13
#define M_StdDev		14
#define M_CoilNum		15


/*
**
**	Function Prototypes
**
*/

TOKAMAK *FileInput(const char *);

#ifdef __cplusplus
}
#endif

#endif
