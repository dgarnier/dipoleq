/* cpdfMemBuf.c -- Memory Stream --- for ClibPDF v2.00-r1 or later
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

	Write to memory stream as if writing to file via fprintf(), fwrite() etc.
	without worrying about how long the memory buffer is (no buffer over-runs).
	The buffer will be expanded as needed to accomodate more data.
	At the end, the entire buffer content will be available with the address
	of the buffer, and the length of the content.
	Modelled after NeXT's streams, but I implemented write related functions only,
	i.e., equivalents of fscanf(), fread() are not available.
	See near the end of file for usage example.

1999-09-18 [io] Windows PDF output to stdout fixed (using setmode()).

1999-08-17 [IO] 
	## MULTI_THREADING_OK ##

1998-08-27 [IO]

For testing this module stand-alone (with main() for testing):
cc -Wall -g -DMAINDEF -o cpdfMemBuf cpdfMemBuf.c

---
*/

#include "config.h"
#include "version.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "cpdflib.h"		/* This must be included before all other local include files */
#include "cpdf_private.h"

#if defined(_WIN32) || defined(WIN32)
#include <fcntl.h>
#endif

#ifdef MacOS8
#include <stat.mac.h>
#include <unix.mac.h>
#include <StandardFile.h>
#include <Memory.h>
#include <Files.h>
#include <Strings.h>	/* c2pstr() */
/* extern size_t strlen(const char *s);
extern int strcmp(const char *s1, const char *s2);
extern void *memcpy(void *s1, const void *s2, size_t n);
*/
#endif


#define MEMSTREAM_MAGIC		0xa5b5cafe
#define MINBUFSIZE 		256
/* don't increase buffer size by more than this much at a time */
#define MAXBUFBUMP 		65536

static char *memErrorFmt = "Corrupted memory stream: 0x%lx (%lu)";

/* For write operation, address must be NULL, and size = 0 */
CPDFmemStream *cpdf_openMemoryStream(void)
{
CPDFmemStream *aMstrm;
	aMstrm = (CPDFmemStream *) malloc((size_t)sizeof(CPDFmemStream));
	_cpdf_malloc_check((void *)aMstrm);
	aMstrm->magic_number = (unsigned long)MEMSTREAM_MAGIC;
	aMstrm->buffer = (char *)malloc((size_t)MINBUFSIZE);
	_cpdf_malloc_check((void *)aMstrm->buffer);
 	aMstrm->bufSize = MINBUFSIZE;
	aMstrm->count = 0;
	return aMstrm;
}

/* clears content of memory stream while keeping the buffer as is */
void cpdf_clearMemoryStream(CPDFmemStream *aMstrm)
{
	if(aMstrm) {
	    if(aMstrm->magic_number != MEMSTREAM_MAGIC)
	        cpdf_GlobalError(3, "ClibPDF", memErrorFmt, (long unsigned)aMstrm, (long unsigned)aMstrm);
	    aMstrm->count = 0;
	    *aMstrm->buffer = 0;	/* null terminate */
	}
}

void cpdf_closeMemoryStream(CPDFmemStream *aMstrm)
{
	/* alex_v: should check on NULL pointer */
	if(aMstrm) {
	    if(aMstrm->magic_number != MEMSTREAM_MAGIC)
		cpdf_GlobalError(3, "ClibPDF", memErrorFmt, (long unsigned)aMstrm, (long unsigned)aMstrm);
	    if(aMstrm->buffer) free(aMstrm->buffer);
	    if(aMstrm) free(aMstrm);
	    aMstrm = NULL;
	}
}

/* Write (append) more data of length 'len' to the memory stream */

int cpdf_writeMemoryStream(CPDFmemStream *memStream, const char *data, int len)
{
int errcode = 0;
    /* expand buffer as needed */
    if(memStream->magic_number != MEMSTREAM_MAGIC)
	cpdf_GlobalError(3, "ClibPDF", memErrorFmt, (long unsigned)memStream, (long unsigned)memStream);

    if((len + memStream->count + 1) > memStream->bufSize) {    /* Must expand buffer */
        if (memStream->buffer == NULL) {
            memStream->buffer = (char *) malloc(MINBUFSIZE);
	    _cpdf_malloc_check((void *)memStream->buffer);
            memStream->bufSize = MINBUFSIZE;
        }
        else {
	    while((len + memStream->count + 1) > memStream->bufSize) {
		if(memStream->bufSize < MAXBUFBUMP)
                    memStream->bufSize *= 2;		/* geometrical growth upto a limit */
		else
               	    memStream->bufSize += MAXBUFBUMP;	/* switch to linear increment */
	    }
            memStream->buffer = (char *) realloc(memStream->buffer, memStream->bufSize);
	    _cpdf_malloc_check((void *)memStream->buffer);
        }
    }
    if(errcode) return(errcode);

    /* Put the new data in the buffer and update the count */
    memcpy(memStream->buffer + memStream->count, data, len);
    memStream->count += len;

    /* null-terminate buffer (for string applications) */
    *(memStream->buffer + memStream->count) = 0;
     return(0);
}


int cpdf_memPutc(int ch, CPDFmemStream *memStream)
{
char c1[2];
	c1[0] = ch; c1[1] = 0;
	cpdf_writeMemoryStream(memStream, c1, 1);
	return (1);
}

int cpdf_memPuts(char *str, CPDFmemStream *memStream)
{
int numchar = strlen(str);
	cpdf_writeMemoryStream(memStream, str, numchar);
	return numchar;
}

void cpdf_getMemoryBuffer(CPDFmemStream *memStream, char **streambuf, int *len, int *maxlen)
{
	if(memStream->magic_number != MEMSTREAM_MAGIC)
	    cpdf_GlobalError(3, "ClibPDF", memErrorFmt, (long unsigned)memStream, (long unsigned)memStream);
	*streambuf = memStream->buffer;
	*len = memStream->count;
	*maxlen = memStream->bufSize;
}

/* Using a special filename "-" will send the buffer to stdout */

int cpdf_saveMemoryStreamToFile(CPDFmemStream *memStream, const char *name)
{
FILE *fpout;
char *memBuffer;
int memLen, bufSize;
#if defined(_WIN32) || defined(WIN32)
extern int setmode(int filenum, int mode);
#endif

    if( strcmp(name, "-") == 0) {
#if defined(_WIN32) || defined(WIN32)
	fpout = stdout;
	setmode(fileno(fpout), O_BINARY);
#else
	fpout = fdopen(1, BINARY_WRITE);	/* fd of 1 == stdout */
#endif
    }
    else {
	if((fpout = fopen(name, BINARY_WRITE)) == NULL) {
	    cpdf_GlobalError(1, "ClibPDF", "Cannot open output file: %s", name);
	    return(-1);
	}
    }
    cpdf_getMemoryBuffer(memStream, &memBuffer, &memLen, &bufSize);
    /* fprintf(stderr, "file=%s\nlength=%d, bufferSize=%d\n", name, memLen, bufSize); */
    fwrite(memBuffer, 1, (size_t)memLen, fpout);
    fclose(fpout);

    return(0);
}


/* memory stream debug dump function */

void _checkMemMagic(char *idstr, CPDFmemStream *memStream)
{
    if(memStream->magic_number != MEMSTREAM_MAGIC) {
	cpdf_GlobalError(3, idstr, "stream=%lu, magic= %lu, buffer= %lu, count= %d, bufSize=%d",
	    (long unsigned)memStream, (long unsigned)memStream->magic_number,
	    (long unsigned)memStream->buffer, memStream->count, memStream->bufSize);
    }
}


#ifdef MAINDEF
static char *string1 = "0123456789 This is a test of memory stream 0123456789";

void main(void)
{
char sbuf[256];
int i;
CPDFmemStream *memStream = { NULL };
char *memBuffer;
int memLen, bufSize;
	memStream = cpdf_openMemoryStream();
	for(i=0; i<500; i++) {
	    sprintf(sbuf, "%05d - %s\n", i, string1);
	    cpdf_writeMemoryStream(memStream, sbuf, strlen(sbuf));
	}
	cpdf_getMemoryBuffer(memStream, &memBuffer, &memLen, &bufSize);
	cpdf_saveMemoryStreamToFile(memStream, "junkmem.out");
	fprintf(stderr,"memLen=%d,  bufSize=%d\n", memLen, bufSize);
	cpdf_closeMemoryStream(memStream);
}


#endif

