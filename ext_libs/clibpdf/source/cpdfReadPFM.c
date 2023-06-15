/* cpdfReadPFM.c -- PFM (printer font metrics) file read module
 * Copyright (C) 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

## MULTI_THREADING_OK  ##

For testing this module stand-alone (with main() for testing):
cc -Wall -g -DNEXTSTEP -DMAINDEF -o cpdfReadPFM cpdfReadPFM.c

NOTE: "-fpack-struct" is no longer needed. Just don't save the PFM structs in one write.
cc -Wall -g -fpack-struct -DNEXTSTEP -DMAINDEF -o cpdfReadPFM cpdfReadPFM.c

cpdfReadPFM hv______.pfm
cpdfReadPFM a010013l.pfm


Sheesh. Gave up on trying to make compiler pack structures.  This may be
done on a given platform, but it is going to be a nightmare for cross-platform
portability.  It may be slower, but we will just read invdividual items
within a struct one-by-one.

#### WARNING ###################################################################

	PFM structure design is screwed.  Many compilers will have problems
	with byte alignment requirement unless #pragma or
	compiler option to use BYTE alignment is enabled.  You must find out
	how to do this for each compiler.  Please see cpdfPFM.h for details.
	Check sizeof(CPDF_PFM_Header) by testing with the main() at the end
	which is enabled by defining MAINDEF.

For Windows/Microsoft env use:
#pragma pack(1)
   packing enabled here
#gragma pack()

For GCC (GNU C compiler: gcc version 2.7.2.1 or later):
Use compiler option "-fpack-struct" <-- FORGET IT.

For other compilers, you must find a way to achieve packed structures for it.
--
1999-08-15 [io] Should be multi-threading ready.  No globals used.
		Must use cpdf_GolbalError().

1999-03-12 [IO]

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

/* #include "endian.h" */
#include "config.h"
#include "version.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "cpdflib.h"		/* This must be included before all other local include files */
#include "cpdf_private.h"
#include "cpdfPFM.h"

#define BUFSIZE 8192

int _cpdf_readPFM(const char *pfmfile, CPDFfontInfo *fI, CPDFfontDescriptor *fD);
int _cpdf_readPFM_Header(FILE *fp, CPDF_PFM_Header *pfm);
int _cpdf_readPFM_Extension(FILE *fp, CPDF_PFM_Extension *psx);
void _cpdf_byteRevPFM_Header(CPDF_PFM_Header *pfm);
void _cpdf_byteRevPFM_Extension(CPDF_PFM_Extension *pfmx);
unsigned long LbyteRev(unsigned long Vin);
unsigned short SbyteRev(unsigned short Vin);

int _cpdf_readPFM(const char *pfmfile, CPDFfontInfo *fI, CPDFfontDescriptor *fD)
{
int len=0, isMono, i;
int byteRev = _Big_Endian();	/* PFM files are stored in the little-endian format. */
FILE *fp;
CPDF_PFM_Header *pfm, pfmbuf;
CPDF_PFM_Extension *psx, psxbuf;
char lbuf[256], fntemp[64];

    pfm = &pfmbuf;
    psx = &psxbuf;

    /* Open the file */
    if((fp = fopen(pfmfile, BINARY_READ)) == NULL) {
#ifdef MAINDEF
        fprintf(stderr, "Unable to open input file: %s\n", pfmfile);
	exit(1);
#else
	/* cpdf_GlobalError(  1, "ClibPDF", "ReadPFM - Unable to open input file: %s", pfmfile); */
	return(1);
#endif
    }
    /* File opened OK */

    /* Read the file */
    _cpdf_readPFM_Header(fp, pfm);
    if(byteRev)
	_cpdf_byteRevPFM_Header(pfm);		/* endian adjustment */

    /* Do consistency check */
    if(len != (int)pfm->len &&  			/* Check length field matches file length */
        pfm->extlen != 30 &&     			/* Check length of PostScript extension   */
        pfm->fontname > 75 && pfm->fontname < 512) {	/* Font name specified */
#ifdef MAINDEF
        fprintf(stderr, "Not a valid Windows Type 1 .pfm file: %s\n", pfmfile);
        exit(2);
#else
	cpdf_GlobalError(  1, "ClibPDF", "ReadPFM - Not a valid Type1 PFM file:  %s", pfmfile);
	fclose(fp);
	return(2);
#endif
    }  /* end of consistency check */

    /* GOOD - We have a OK PFM file, copy values to structs passed from caller */
    /* Now, seek to offset pfm->psext */

    fseek(fp, (long int)pfm->psext, SEEK_SET);
    _cpdf_readPFM_Extension(fp, psx);
    if(byteRev)
	_cpdf_byteRevPFM_Extension(psx);		/* endianness */

#if defined(MAINDEF) || defined(CPDF_DEBUG)
    printf("Copyright: <%s>\n", pfm->copyright);
    /* printf("Font name: <%s>\n", (char *)(pfm+pfm->fontname)); */
    printf("Char tab offset: %ld\n", pfm->chartab);
#endif

    /* seek to font name */
    fseek(fp, (long int)pfm->fontname, SEEK_SET);
    fread(lbuf, 1, sizeof(lbuf), fp);
    sscanf(lbuf, "%s", fntemp);
    fD->fontName = (char *)malloc(sizeof(fntemp)+1);
    _cpdf_malloc_check((void *)fD->fontName);
    strcpy(fD->fontName, fntemp);
#ifdef	CPDF_DEBUG
    fprintf(stderr, "Font found: %s\n", fD->fontName);
#endif

    /*
     *  The mono flag in the pfm actually indicates whether there is a
     *  table of font widths, not if they are all the same.
     */
    if (!(pfm->kind&1) ||			/* Flag for mono */
        pfm->avgwidth == pfm->maxwidth ) {	/* Avg width = max width */
        isMono = 1;
    } else {
        isMono = 0;
    }

    /* Assign values to Font Descriptor */
    fD->ascent = pfm->ascent;
    fD->capHeight = psx->capheight;
    fD->descent = psx->descender;
    fD->flags = 6;	/* can't get this info from PFM. Doesn't matter if font embedded */
    /*
     * The font bounding box is lost, but try to reconstruct it.
     * Much of this is just guess work.  The bounding box is required in
     * the .afm, but is not used by the PM font installer.
     */
    if(isMono)
        fD->fontBBox[0] = -20;
    else
	fD->fontBBox[0] = -100;
    fD->fontBBox[1] = -(psx->descender+5);  /* Descender is given as positive value */
    fD->fontBBox[2] = pfm->maxwidth+10;
    fD->fontBBox[3] = pfm->ascent+5;

    /* fD->fontname = ; Caller can stuff this. */
    fD->italicAngle = (float)(psx->slant)/10.0;
    fD->stemV = 80;		/* not in PFM */
    fD->xHeight = psx->xheight;

    /* Assign values to entries in Font Info struct */
    fI->firstchar = pfm->firstchar;
    fI->lastchar = pfm->lastchar;
    fI->descLevel = 2;

    /* If the charset is not the Windows standard, just put out
     * unnamed entries.
     * Now, seek to offset pfm->chartab.
     */
    fI->width = (short *)malloc(sizeof(short)*268);	/* make sure to deallocate (256 or 268 ?) */
    _cpdf_malloc_check((void *)fI->width);
    for(i=0; i< fI->firstchar; i++)
	fI->width[i] = 0;
    for(i= fI->lastchar; i< 268; i++)
	fI->width[i] = 0;

    if(fseek(fp, (long int)pfm->chartab, SEEK_SET)) {
	fclose(fp);
	return(3);
    }
    len = fI->lastchar - fI->firstchar +1;		/* number of widths to read */
    if(fread(fI->width + fI->firstchar, sizeof(unsigned short), len, fp) != (unsigned long)len) {
	fclose(fp);
	return(4);
    }
    if(byteRev) {
	for(i= fI->firstchar; i<= fI->lastchar; i++)
	    fI->width[i] = SbyteRev(fI->width[i]);
    }
    fclose(fp);
    return(0);
}


/* We give up on packing structure to overcome Microsoft design.
 * It is just too messy for making this work for all platforms.
 * Just read each item separately.
*/
int _cpdf_readPFM_Header(FILE *fp, CPDF_PFM_Header *pfm)
{
int bytecount = 0;
    if(sizeof(CPDF_PFM_Header) == PFM_SIZE) {
	/* the struct is packed */
	bytecount = fread(pfm, 1, sizeof(CPDF_PFM_Header), fp);
    }
    else {
	/* struct did not get packed by compiler */
	bytecount += fread(&(pfm->vers), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(pfm->len), 1, sizeof(unsigned long), fp);
 	bytecount += fread(&(pfm->copyright), sizeof(char), 60, fp);
	bytecount += fread(&(pfm->type), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(pfm->points), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(pfm->verres), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(pfm->horres), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(pfm->ascent), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(pfm->intleading), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(pfm->extleading), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(pfm->italic), 1, sizeof(unsigned char), fp);
	bytecount += fread(&(pfm->uline), 1, sizeof(unsigned char), fp);
	bytecount += fread(&(pfm->overs), 1, sizeof(unsigned char), fp);
	bytecount += fread(&(pfm->weight), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(pfm->charset), 1, sizeof(unsigned char), fp);
	bytecount += fread(&(pfm->pixwidth), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(pfm->pixheight), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(pfm->kind), 1, sizeof(unsigned char), fp);
	bytecount += fread(&(pfm->avgwidth), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(pfm->maxwidth), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(pfm->firstchar), 1, sizeof(unsigned char), fp);
	bytecount += fread(&(pfm->lastchar), 1, sizeof(unsigned char), fp);
	bytecount += fread(&(pfm->defchar), 1, sizeof(unsigned char), fp);
	bytecount += fread(&(pfm->brkchar), 1, sizeof(unsigned char), fp);
	bytecount += fread(&(pfm->widthby), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(pfm->device), 1, sizeof(unsigned long), fp);
	bytecount += fread(&(pfm->face), 1, sizeof(unsigned long), fp);
	bytecount += fread(&(pfm->bits), 1, sizeof(unsigned long), fp);
	bytecount += fread(&(pfm->bitoff), 1, sizeof(unsigned long), fp);
	bytecount += fread(&(pfm->extlen), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(pfm->psext), 1, sizeof(unsigned long), fp);
	bytecount += fread(&(pfm->chartab), 1, sizeof(unsigned long), fp);
	bytecount += fread(&(pfm->res1), 1, sizeof(unsigned long), fp);
	bytecount += fread(&(pfm->kernpairs), 1, sizeof(unsigned long), fp);
	bytecount += fread(&(pfm->kerntrack), 1, sizeof(unsigned long), fp);
	bytecount += fread(&(pfm->fontname), 1, sizeof(unsigned long), fp);
    }
#if defined(MAINDEF) || defined(CPDF_DEBUG)
    fprintf(stderr,"Header section read %d bytes\n", bytecount);
#endif
    return(bytecount);
}

int _cpdf_readPFM_Extension(FILE *fp, CPDF_PFM_Extension *psx)
{
int bytecount = 0;
    if(sizeof(CPDF_PFM_Extension) == PFM_SIZE) {
	/* the struct is packed */
	bytecount = fread(psx, 1, sizeof(CPDF_PFM_Extension), fp);
    }
    else {
	/* the struct did not get packed by compiler */
	bytecount += fread(&(psx->len), 1, sizeof(unsigned short), fp);
 	bytecount += fread(&(psx->res1), sizeof(char), 12, fp);
	bytecount += fread(&(psx->capheight), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(psx->xheight), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(psx->ascender), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(psx->descender), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(psx->slant), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(psx->superscript), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(psx->subscript), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(psx->superscriptsize), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(psx->subscriptsize), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(psx->underlineoffset), 1, sizeof(unsigned short), fp);
	bytecount += fread(&(psx->underlinewidth), 1, sizeof(unsigned short), fp);
   }
#if defined(MAINDEF) || defined(CPDF_DEBUG)
    fprintf(stderr,"Extension section read %d bytes\n", bytecount);
#endif
    return(bytecount);
}


/* ===== Utility Functions to deal with endian conversion =============================== */
/* These are used to reverse the byte order if the files byte-order does not
 * match the CPU's byte order.
 */

void _cpdf_byteRevPFM_Header(CPDF_PFM_Header *pfm)
{
    pfm->vers = SbyteRev(pfm->vers);
    pfm->len  = LbyteRev(pfm->len);
    pfm->type = SbyteRev(pfm->type);
    pfm->points = SbyteRev(pfm->points);
    pfm->verres = SbyteRev(pfm->verres);
    pfm->horres = SbyteRev(pfm->horres);
    pfm->ascent = SbyteRev(pfm->ascent);
    pfm->intleading = SbyteRev(pfm->intleading);
    pfm->extleading = SbyteRev(pfm->extleading);
    pfm->weight = SbyteRev(pfm->weight);
    pfm->pixwidth = SbyteRev(pfm->pixwidth);
    pfm->pixheight = SbyteRev(pfm->pixheight);
    pfm->avgwidth = SbyteRev(pfm->avgwidth);
    pfm->maxwidth = SbyteRev(pfm->maxwidth);
    pfm->widthby = SbyteRev(pfm->widthby);
    pfm->extlen = SbyteRev(pfm->extlen);
    pfm->device  = LbyteRev(pfm->device);
    pfm->face  = LbyteRev(pfm->face);
    pfm->bits  = LbyteRev(pfm->bits);
    pfm->bitoff  = LbyteRev(pfm->bitoff);
    pfm->psext  = LbyteRev(pfm->psext);
    pfm->chartab  = LbyteRev(pfm->chartab);
    pfm->res1  = LbyteRev(pfm->res1);
    pfm->kernpairs  = LbyteRev(pfm->kernpairs);
    pfm->kerntrack  = LbyteRev(pfm->kerntrack);
    pfm->fontname  = LbyteRev(pfm->fontname);
}

void _cpdf_byteRevPFM_Extension(CPDF_PFM_Extension *pfmx)
{
    pfmx->len = SbyteRev(pfmx->len);
    pfmx->capheight = SbyteRev(pfmx->capheight);
    pfmx->xheight = SbyteRev(pfmx->xheight);
    pfmx->ascender = SbyteRev(pfmx->ascender);
    pfmx->descender = SbyteRev(pfmx->descender);
    pfmx->slant = SbyteRev(pfmx->slant);
    pfmx->superscript = SbyteRev(pfmx->superscript);
    pfmx->subscript = SbyteRev(pfmx->subscript);
    pfmx->superscriptsize = SbyteRev(pfmx->superscriptsize);
    pfmx->subscriptsize = SbyteRev(pfmx->subscriptsize);
    pfmx->underlineoffset = SbyteRev(pfmx->underlineoffset);
    pfmx->underlinewidth = SbyteRev(pfmx->underlinewidth);
}


unsigned long LbyteRev(unsigned long Vin)
{
unsigned long temp;
    temp=(Vin<<24)|((Vin&0xff00)<<8)|((Vin&0xff0000)>>8)|((Vin&0xff000000)>>24);
    return(temp);
}

unsigned short SbyteRev(unsigned short Vin)
{
unsigned short temp;
    temp = (Vin>>8)|((Vin&0x00ff)<<8);
    return(temp);
}

/* ====================================================================================== */
/*   ONLY FOR TESTING */

#ifdef MAINDEF
void Usage(char *pname)
{
    fprintf(stderr, "Usage: %s pfm_filename\n", pname);
    exit(1);
}

void main(int argc, char *argv[])
{
int i;
CPDFfontInfo fntI;
CPDFfontDescriptor fntD;

    printf("Size of PFM Header:    %ld  (should be 143)\n", (long)sizeof(CPDF_PFM_Header));
    printf("Size of PFM Extension: %ld  (should be 36)\n\n", (long)sizeof(CPDF_PFM_Extension));

    if(argc != 2)
	Usage(argv[0]);

    // fntI.width = (short *)malloc(sizeof(short)*256);

    _cpdf_readPFM(argv[1], &fntI, &fntD);

    printf("\nFontBBox [%d %d %d %d]\n", fntD.fontBBox[0], fntD.fontBBox[1], fntD.fontBBox[2], fntD.fontBBox[3]);
    printf("Font width array = [\n");
    for(i=0; i<= 255; i++) {
	printf("%d ", fntI.width[i]);
	if((i+1)%10 == 0) putchar('\n');
    }
    printf("]\n");
}

#endif		/* MAINDEF */
