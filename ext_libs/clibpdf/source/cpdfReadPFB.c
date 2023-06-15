/* cpdfReadPFB.c - Read PFB (binary font file)
 * Copyright (C) 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

## MULTI_THREADING_OK  ##

Format of PFB file:
---------------------------------------------------------------------------------
    uchar  magic	(must be 128 or 0x80)
    uchar  type		(1=ASCII, 2=BINARY): this block should have type = 1
    ulong  length	(little endian length of ASCII portion)
     ... ASCII portion of PFB file ...
    uchar  magic	(must be 128 or 0x80)
    uchar  type		(1=ASCII, 2=BINARY): this block should have type = 2
    ulong  length	(little endian length of BINARY portion)
     ... BINARY portion of PFB file ...
    uchar  magic	(must be 128 or 0x80)
    uchar  type		(1=ASCII, 2=BINARY): this block should have type = 1
    ulong  length	(little endian length of tail in ASCII portion)
     ... 512 '0's (line breaks not counted) followed by
	"cleartomark" and possibly some more data (?)
---------------------------------------------------------------------------------

1999-08-15 [io] Should be multi-threading ready.  No globals used.
		Must use cpdf_GolbalError().

1999-06-12 [io]	Rewritten with no third-party code.
*/

/* #include "endian.h" */
#include "config.h"
#include "version.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "cpdflib.h"		/* This must be included before all other local include files */
#include "cpdf_private.h"

#define PFB_ASCII	1
#define PFB_BINARY	2
#define PFB_END		3
#define PFB_MAGIC	128

unsigned long LbyteRev(unsigned long Vin);
static int _check_PFBmagic_and_type(FILE *fp, const char *pfile, int pfb_type, int section);
static unsigned char *_cr_to_lf(unsigned char *buf, unsigned long length);



int _cpdf_readPFB(const char *pfbfile, CPDFextFont *eFI)
{
int  retval = 0;
unsigned long length;
unsigned char *buf;
FILE *fpb;
char *memBuffer;
int memLen, bufSize;
int byteRev = _Big_Endian();	/* PFM files are stored in the little-endian format. */

	eFI->fontMemStream = cpdf_openMemoryStream();	/* memory stream to hold font data */
	/* Open the file */
	if((fpb = fopen(pfbfile, BINARY_READ)) == NULL) {
	    cpdf_GlobalError(  1, "ClibPDF", "ReadPFB - Unable to open PFB file: %s", pfbfile);
	    retval =1;
	    goto err_exit;
	}

	/* === Intial ASCII section of PFB font */
	if((retval = _check_PFBmagic_and_type(fpb, pfbfile, PFB_ASCII, 1)))	/* check magic and type */
	    goto err_exit;
	fread(&length, sizeof(unsigned long), 1, fpb);
	if(byteRev)
	    length = LbyteRev(length);
	buf = (unsigned char *) malloc(length);
	_cpdf_malloc_check((void *)buf);
	if(fread(buf, 1, (size_t)length, fpb) != length) {
	    cpdf_GlobalError(  1, "ClibPDF", "ReadPFB - ASCII section-1 read error: %s", pfbfile);
	    retval =4;
	    free(buf);
	    goto err_exit;
	}
	_cr_to_lf(buf, length);
	eFI->length1 = length;
	cpdf_writeMemoryStream(eFI->fontMemStream, (char *)buf, length);
	free(buf);
    
	/* === BINARY section of PFB font -- actual font data */
	if((retval = _check_PFBmagic_and_type(fpb, pfbfile, PFB_BINARY, 2)))	/* check magic and type */
	    goto err_exit;
	fread(&length, sizeof(unsigned int), 1, fpb);
	if(byteRev)
	    length = LbyteRev(length);
	buf = (unsigned char *) malloc(length);
	_cpdf_malloc_check((void *)buf);
	if(fread (buf, 1, (size_t)length, fpb) != length) {
	    cpdf_GlobalError(  1, "ClibPDF", "ReadPFB - BINARY section-2 read error: %s", pfbfile);
	    retval =5;
	    free(buf);
	    goto err_exit;
	}
	eFI->length2 = length;
	cpdf_writeMemoryStream(eFI->fontMemStream, (char *)buf, length);
	free(buf);
    
	/* === Last ASCII section of PFB font (512 zeros and "cleartomark"). */
	if((retval = _check_PFBmagic_and_type(fpb, pfbfile, PFB_ASCII, 3)))	/* check magic and type */
	    goto err_exit;
	fread(&length, sizeof(unsigned long), 1, fpb);
	if(byteRev)
	    length = LbyteRev(length);
	buf = (unsigned char *) malloc(length);
	_cpdf_malloc_check((void *)buf);
	if(fread(buf, 1, (size_t)length, fpb) != length) {
	    cpdf_GlobalError(  1, "ClibPDF", "ReadPFB - ASCII section-3 read error: %s, %d", pfbfile);
	    retval =6;
	    free(buf);
	    goto err_exit;
	}
	_cr_to_lf(buf, length);
	eFI->length3 = length;
	cpdf_writeMemoryStream(eFI->fontMemStream, (char *)buf, length);
	free(buf);
    
	/* All sections successfully read. */
	fclose(fpb);	/* Close PFB file */
	cpdf_getMemoryBuffer(eFI->fontMemStream, &memBuffer, &memLen, &bufSize);
	eFI->length = memLen;	/* total length of font data */
	return(0);

err_exit:
	fclose(fpb);
	if(eFI->fontMemStream)
	    cpdf_closeMemoryStream(eFI->fontMemStream);
	return(retval);
}

static int _check_PFBmagic_and_type(FILE *fp, const char *pfile, int pfb_type, int section)
{
int abtype = 0;
	if((abtype = getc(fp)) != PFB_MAGIC) {
	    cpdf_GlobalError(  1, "ClibPDF",
		"ReadPFB - Bad magic number: %d (128 expected) in section %d of file %s",
		abtype, pfile, section);
	    return(2);
	}
	/* Check ASCII or BINARY type */
	if((abtype = getc(fp)) != pfb_type) {
	    cpdf_GlobalError(  1, "ClibPDF",
		"ReadPFB - Unexpected type: %d (%d expected) in section %d of file %s",
		abtype, pfb_type, section, pfile);
	    return(3);
	}
	return(0);	/* Magic OK, and expected type */
}

static unsigned char *_cr_to_lf(unsigned char *buf, unsigned long length)
{
unsigned long i;
char *p;
	p = (char *)buf;
	for(i = 0; i < length; i++)
	    if( *p == '\r') *p++ = '\n';
	return(buf);
}

