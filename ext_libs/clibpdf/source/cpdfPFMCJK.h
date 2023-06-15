/* cpdfPFM.h -- PFM (printer font metrics) header
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.txt or LICENSE.pdf.

1999-03-12 [IO]
----
For PFM structure specification, see:
Adobe Technical Note #5178 "Building PFM Files for PostScript-Language CJK Fonts"

*/

/* PFM Header: 117 bytes
#### WARNING ###################################################################
	This struct may come out occupying 120 bytes unless #pragma or
	compiler option to use BYTE alignment is enabled.
	Who created this incompetent struct?
	Check sizeof(CPDF_PFM_Header) by testing with the main() at the end
	which is enabled by defining MAINDEF.

For Windows/Microsoft env use:
#pragma pack(1)
   packing enabled here
#gragma pack()

For GCC (GNU C compiler), use:
compiler option "-fpack-struct".

For other compilers, you must find a way to achieve packed structures for it.
################################################################################

*/
#pragma pack(1)

typedef struct {
    short dfVersion;
    long  dfSize;
    char  dfCopyright[60];
    short dfType;
    short dfPoint;
    short dfVertRes;
    short dfHorizRes;
    short dfAscent;
    short dfInternalLeading;
    short dfExternalLeading;
    unsigned char dfItalic;
    unsigned char dfUnderline;
    unsigned char dfStrikeOut;
    short dfWeight;
    unsigned char dfCharSet;
    short dfPixWidth;
    short dfPixHeight;
    unsigned char dfPitchAndFamily;
    short dfAvgWidth;
    short dfMaxWidth;
    unsigned char dfFirstChar;
    unsigned char dfLastChar;
    unsigned char dfDefaultChar;
    unsigned char dfBreakChar;
    short dfWidthBytes;
    long dfDevice;
    long dfFace;
    long dfBitsPointer;
    long dfBitsOffset;
} CPDF_PFM_Header;
#pragma pack()

/* PFM Extension: 30 bytes */
typedef struct {
    short dfSizeFields;
    long  dfExtMetricsOffset;
    long  dfExtentTable;
    long  dfOriginTable;
    long  dfPairKernTable;
    long  dfTrackKernTable;
    long  dfDriverInfo;
    long  dfReserved;
} CPDF_PFM_Extension;

/* PFM Extended Text Metrics: 52 bytes */
typedef struct {

} CPDF_PFM_ExtTextMetrics;

/* PFM PostScript Information */
typedef struct {
    /* actually 3 pointers to char are probably not what I need */
    char *device_type;		/* static string "PostScript" */
    char *windows_name;		/* Windows' menu name as it should appear in app's font menus */
    char *postscript_name;	/* fully qualified PS font name */
} CPDF_PFM_PSInfo;

/* PFM Extent Table: 448 byes */
typedef struct {
    short width[224];
} CPDF_PFM_ExtTab;

