/* zlibcomp.c  - Zlib/Flate compression program.
See http://www.cdrom.com/pub/infozip/zlib/
for the library and header files.
cc -o zlibcomp zlibcomp.c -lz
------------
This is a highly stripped and edited version of:
   example.c -- usage example of the zlib compression library
 * Copyright (C) 1995-1998 Jean-loup Gailly.
 * For conditions of distribution and use, see copyright notice in zlib.h */

#include <stdio.h>
#include "zlib.h"
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

#if defined(_WIN32) || defined(__MWERKS__)
#define BINARY_WRITE	"wb"
#define BINARY_READ	"rb"
#else
#define UNIX 		1
#define BINARY_WRITE	"w"
#define BINARY_READ	"r"
#endif

#define CHECK_ERR(err, msg) { \
    if (err != Z_OK) { \
        fprintf(stderr, "%s error: %d\n", msg, err); \
        exit(1); \
    } \
}

void main(int argc, char *argv[]) {
unsigned long	filesize;
FILE *fpi, *fpo;
#ifndef _WIN32
struct stat filestat;
#else
struct _stat filestat;		/* For Windows VC++ */
#endif
Byte *inputBuf, *outputBuf;
uLong comprLen = 0;
static const char* myVersion = ZLIB_VERSION;
char ofile[1026];
int err;
	if(argc != 2) {
	    fprintf(stderr, "%s afile\n", argv[0]);
	    fprintf(stderr, "  Note: afile will be deleted and replaced by the compressed version afile.zlib\n");
	    exit(1);
	}
	if (zlibVersion()[0] != myVersion[0]) {
	    fprintf(stderr, "incompatible zlib version\n");
	    exit(1);
	} else if (strcmp(zlibVersion(), ZLIB_VERSION) != 0) {
	    fprintf(stderr, "warning: different zlib version\n");
	}
#ifndef _WIN32
	if(stat(argv[1], &filestat))
#else
	if(_stat(argv[1], &filestat))	/* Windows */
#endif
	    filesize = 0;		/* stat failed -- no such file */
	else {
	    if( !(filestat.st_mode & S_IFREG) )
		filesize = 0;	/* not a regular file */
#ifdef UNIX
	    else if( !(filestat.st_mode & S_IREAD) )
		filesize = 0;	/* not readable */
#endif
	    else
		filesize = filestat.st_size;
	}
	inputBuf = (Byte *)malloc((size_t)(filesize+16));
	comprLen = filesize+1024;
	outputBuf = (Byte *)malloc((size_t)comprLen);
	if((fpi = fopen(argv[1], BINARY_READ)) == NULL) {
	    fprintf(stderr, "Cant open input file: %s\n", argv[1]);
	    exit(1);
	}
	strncpy(ofile, argv[1], 1000);
	strcat(ofile, ".zlib");
	if((fpo = fopen(ofile, BINARY_WRITE)) == NULL) {
	    fprintf(stderr, "Cant open output file: %s\n", ofile);
	    fclose(fpi);
	    exit(1);
	}
	fread(inputBuf, 1, (size_t)filesize, fpi);
	fclose(fpi);
	err = compress(outputBuf, &comprLen, inputBuf, filesize);
    	CHECK_ERR(err, "zlib compress");
	fwrite(outputBuf, 1, (size_t)comprLen, fpo);
	fclose(fpo);
	remove(argv[1]);
        exit(0);
}
