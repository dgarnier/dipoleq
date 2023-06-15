/* cpdfPFM.h -- PFM (printer font metrics) header
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

1999-03-12 [IO]
----
For PFM structure specification, see:
Windows Device Driver Delopver Kit:
Also there is:
Adobe Technical Note #5178 "Building PFM Files for PostScript-Language CJK Fonts"

Aargh. Gave up on trying to make compiler pack structures.  This may be
done on a given platform, but is going to be a nightmare for cross-platform
portability.  It may be slower, but we will just read invdividual items
within a struct one-by-one.
*/

/*
#### WARNING ###################################################################
1999-12-04 [io]  Struct packing pragma's disabled for Sun CC

	This struct may come out occupying wrong number of bytes
	unless #pragma or compiler option to use BYTE alignment is enabled.
	Check sizeof(CPDF_PFM_Header) by testing with the main() at the end
	which is enabled by defining MAINDEF.

For Windows/Microsoft env use:
#pragma pack(1)
   packing enabled here
#gragma pack()

For GCC (GNU C compiler):
Use compiler option "-fpack-struct" <-- FORGET IT.

For other compilers, you must find a way to achieve packed structures for it.
################################################################################
Parts of this file have been derived from the public domain pfm2afm.c,
the original version of which is available from:
ftp://tug.ctan.org/pub/tex-archive/fonts/utilities/
 ********************************************************************
 *                                                                  *
 *  Title:  pfm2afm - Convert Windows .pfm files to .afm files      *
 *                                                                  *
 *  Author: Ken Borgendale   10/9/91  Version 1.0                   *
 *                                                                  *
 *  Function:                                                       *
 *      Convert a Windows .pfm (Printer Font Metrics) file to a     *
 *      .afm (Adobe Font Metrics) file.  The purpose of this is     *
 *      to allow fonts put out for Windows to be used with OS/2.    *
 *                                                                  *
 *  Syntax:                                                         *
 *      pfm2afm  infile  [outfile] -a                               *
 *                                                                  *
 *  Copyright:                                                      *
 *      pfm2afm - Copyright (C) IBM Corp., 1991                     *
 *                                                                  *
 *      This code is released for public use as long as the         *
 *      copyright remains intact.  This code is provided asis       *
 *      without any warrenties, express or implied.                 *
 *                                                                  *
 *  Notes:                                                          *
 *      1. Much of the information in the original .afm file is     *
 *         lost when the .pfm file is created, and thus cannot be   *
 *         reconstructed by this utility.  This is especially true  *
 *         of data for characters not in the Windows character set. *
 *                                                                  *
 *      2. This module is coded to be compiled by the MSC 6.0.      *
 *         For other compilers, be careful of the packing of the    *
 *         PFM structure.                                           *
 *                                                                  *
 ********************************************************************
 ********************************************************************
 *  Modified:  Russell Lang <rjl@eng.monash.edu.au>                 *
 *             1994-01-06  Version 1.1                              *
 *  Compiles with EMX/GCC                                           *
 *  Changed to AFM 3.0                                              *
 *  Put PFM Copyright in Notice instead of Comment                  *
 *  Added ItalicAngle                                               *
 *  Added UnderlinePosition                                         *
 *  Added UnderlineThickness                                        *
 *                                                                  *
 *  Modified 1995-03-10  rjl (fixes from Norman Walsh)              *
 *  Dodge compiler bug when creating descender                      *
 ********************************************************************

*/

#ifndef SolarisSunCC
#pragma pack(1)
#endif

/* PFM Header: 143 bytes */
#define PFM_SIZE 143

typedef struct {
    unsigned short  vers;
    unsigned long   len;             /* Total length of .pfm file */
    unsigned char   copyright[60];   /* Copyright string */
    unsigned short  type;
    unsigned short  points;
    unsigned short  verres;
    unsigned short  horres;
    unsigned short  ascent;
    unsigned short  intleading;
    unsigned short  extleading;
    unsigned char   italic;
    unsigned char   uline;
    unsigned char   overs;
    unsigned short  weight;
    unsigned char   charset;         /* 0=windows, otherwise nomap */
    unsigned short  pixwidth;        /* Width for mono fonts */
    unsigned short  pixheight;
    unsigned char   kind;            /* Lower bit off in mono */
    unsigned short  avgwidth;        /* Mono if avg=max width */
    unsigned short  maxwidth;        /* Use to compute bounding box */
    unsigned char   firstchar;       /* First char in table */
    unsigned char   lastchar;        /* Last char in table */
    unsigned char   defchar;
    unsigned char   brkchar;
    unsigned short  widthby;
    unsigned long   device;
    unsigned long   face;            /* Face name */
    unsigned long   bits;
    unsigned long   bitoff;
    unsigned short  extlen;
    unsigned long   psext;           /* PostScript extension */
    unsigned long   chartab;         /* Character width tables */
    unsigned long   res1;
    unsigned long   kernpairs;       /* Kerning pairs */
    unsigned long   kerntrack;	     /* Track Kern table */ /* rjl */
    unsigned long   fontname;        /* Font name */
} CPDF_PFM_Header;

#ifndef SolarisSunCC
#pragma pack()
#endif

/* PFM Extension: 36 bytes */
#define PSX_SIZE 36

typedef struct {
    unsigned short  len;
    unsigned char   res1[12];
/*  unsigned char   res1[14];           replaced by above two lines by rjl */
    unsigned short  capheight;       /* Cap height */
    unsigned short  xheight;         /* X height */
    unsigned short  ascender;        /* Ascender */
    unsigned short  descender;       /* Descender (positive) */
    /* extra entries added by rjl */
    short  slant;	     /* CW italic angle */
    short  superscript;
    short  subscript;
    short  superscriptsize;
    short  subscriptsize;
    short  underlineoffset;  /* +ve down */
    short  underlinewidth;   /* width of underline */
} CPDF_PFM_Extension;


