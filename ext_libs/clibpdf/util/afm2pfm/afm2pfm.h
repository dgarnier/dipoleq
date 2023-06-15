/*  afm2pfm.h 
 *  afm2pfm - convert Adobe .afm file to Windows binary .pfm file
 *  Copyright (C) 1994  Russell Lang <rjl@eng.monash.edu.au>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/


#define FALSE 0
#define TRUE 1
typedef int BOOL;
typedef unsigned char  BYTE;
typedef unsigned short WORD;
typedef unsigned long  DWORD;

#ifdef NOTUSED
    PFMHEADER       Header;          /* font header */
    WORD            WidthTable[];    /* character widths */
    PFMEXTENSION    Extensions;      /* extensions */
    char            DeviceName[];    /* printer-device name */
    char            FaceName[];      /* font name */
    EXTTEXTMETRIC   ExtTextMetrics;  /* extended-text metrics */
    WORD            ExtentTable[];   /* unscaled-character widths */
    DRIVERINFO      DriverInfo;      /* driver-specific information */
    PAIRKERN        KerningPairs[];  /* pair-kerning table */
    KERNTRACK       KerningTracks[]; /* track-kerning table */
#endif

#pragma pack(1)

typedef struct PFMHEADER {
        WORD  dfVersion;
/*        WORD  dfSize; */
/* c700/drivrkit/document/pfk31.txt says dfSize is of type WORD, 
   but PFM files use DWORD */
	DWORD  dfSize;
        char  dfCopyright[60];
        WORD  dfType;
        WORD  dfPoints;
        WORD  dfVertRes;
        WORD  dfHorizRes;
        WORD  dfAscent;
        WORD  dfInternalLeading;
        WORD  dfExternalLeading;
        BYTE  dfItalic;
        BYTE  dfUnderline;
        BYTE  dfStrikeOut;
        WORD  dfWeight;
        BYTE  dfCharSet;
        WORD  dfPixWidth;
        WORD  dfPixHeight;
        BYTE  dfPitchAndFamily;
        WORD  dfAvgWidth;
        WORD  dfMaxWidth;
        BYTE  dfFirstChar;
        BYTE  dfLastChar;
        BYTE  dfDefaultChar;
        BYTE  dfBreakChar;
        WORD  dfWidthBytes;
        DWORD dfDevice;
        DWORD dfFace;
        DWORD dfBitsPointer;
        DWORD dfBitsOffset;
} PFMHEADER;

typedef struct tagPFMEXTENSION {
        WORD  dfSizeFields;
        DWORD dfExtMetricsOffset;
        DWORD dfExtentTable;
        DWORD dfOriginTable;
        DWORD dfPairKernTable;
        DWORD dfTrackKernTable;
/* c700/drivrkit/document/pfk31.txt says dfDriverInfo should point to 
 * the PostScript name for the font, NOT a DRIVERINFO structure */
        DWORD dfDriverInfo;
        DWORD dfReserved;
} PFMEXTENSION;

typedef struct tagKERNPAIR {
        union {
            BYTE each [2];
            WORD both;
        } kpPair;
        short  kpKernAmount;
} KERNPAIR;

typedef struct tagKERNTRACK {
        short ktDegree;
        short ktMinSize;
        short ktMinAmount;
        short ktMaxSize;
        short ktMaxAmount;
} KERNTRACK;


#ifdef NOTUSED
typedef struct tagDRIVERINFO {
        WORD epSize;      /* size of this data structure */
        WORD epVersion;   /* number indicating structure version */
        DWORD epMemUsage; /* amount of memory font takes up in printer */
        DWORD epEscape;   /* pointer to escape that selects the font */
        TRANSTABLE xtbl;  /* character-set translation information */
} DRIVERINFO;
#endif

typedef struct tagEXTTEXTMETRIC {
        short etmSize;
        short etmPointSize;
        short etmOrientation;
        short etmMasterHeight;
        short etmMinScale;
        short etmMaxScale;
        short etmMasterUnits;
        short etmCapHeight;
        short etmXHeight;
        short etmLowerCaseAscent;
        short etmUpperCaseDescent;
        short etmSlant;
        short etmSuperScript;
        short etmSubScript;
        short etmSuperScriptSize;
        short etmSubScriptSize;
        short etmUnderlineOffset;
        short etmUnderlineWidth;
        short etmDoubleUpperUnderlineOffset;
        short etmDoubleLowerUnderlineOffset;
        short etmDoubleUpperUnderlineWidth;
        short etmDoubleLowerUnderlineWidth;
        short etmStrikeOutOffset;
        short etmStrikeOutWidth;
        WORD  etmKernPairs;
        WORD  etmKernTracks;
} EXTTEXTMETRIC;
#pragma pack()

typedef struct tagBBOX {
	int	llx;
	int	lly;
	int	urx;
	int	ury;
} BBOX;

