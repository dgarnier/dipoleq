/* jpegsize.c -- modified from below to obtain JPEG image size and other
	information from header.  Comments will be ignored.
	For use with inclusion of JPEG images into PDF file.

## MULTI_THREADING_OK  ##

The original README file from Independent JPEG Group's distribution is
included in this directory as README-jpeg.

1999-08-12 [IO, FastIO Systems]
	Eliminated static FILE *infile; for multi-threaded use.

cc -Wall -DJPEG_DEBUG  -o jpegsize jpegsize.c

--------------------------------------------------------------------------
 Original header below
--------------------------------------------------------------------------
 * 
 * rdjpgcom.c
 *
 * Copyright (C) 1994-1997, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains a very simple stand-alone application that displays
 * the text in COM (comment) markers in a JFIF file.
 * This may be useful as an example of the minimum logic needed to parse
 * JPEG markers.
 */

#define JPEG_CJPEG_DJPEG	/* to get the command-line config symbols */
#include "jinclude.h"		/* get auto-config symbols, <stdio.h> */

#include "cpdflib.h"		/* for CPDFimageInfo typedef -- FastIO Systems IO/1998-08-29 */
#include "cpdf_private.h"

#include <ctype.h>		/* to declare isupper(), tolower() */
#ifdef USE_SETMODE
#include <fcntl.h>		/* to declare setmode()'s parameter macros */
/* If you have setmode() but not <io.h>, just delete this line: */
#include <io.h>			/* to declare setmode() */
#endif

#ifdef USE_CCOMMAND		/* command-line reader for Macintosh */
#ifdef __MWERKS__
#include <SIOUX.h>              /* Metrowerks needs this */
#include <console.h>		/* ... and this */
#endif
#ifdef THINK_C
#include <console.h>		/* Think declares it here */
#endif
#endif

#ifdef DONT_USE_B_MODE		/* define mode parameters for fopen() */
#define READ_BINARY	"r"
#else
#ifdef VMS			/* VMS is very nonstandard */
#define READ_BINARY	"rb", "ctx=stm"
#else				/* standard ANSI-compliant case */
#define READ_BINARY	"rb"
#endif
#endif

#ifndef EXIT_FAILURE		/* define exit() codes if not provided */
#define EXIT_FAILURE  1
#endif
#ifndef EXIT_SUCCESS
#ifdef VMS
#define EXIT_SUCCESS  1		/* VMS is very nonstandard */
#else
#define EXIT_SUCCESS  0
#endif
#endif

extern int read_JPEG_header(const char *filename, CPDFimageInfo *jInfo);

/*
 * These macros are used to read the input file.
 * To reuse this code in another application, you might need to change these.
 */


/* Error exit handler */
#define ERREXIT(msg)  (cpdf_GlobalError(1, "ClibPDF jpegsize", msg))


/* Read one byte, testing for EOF */
static int
read_1_byte (FILE *infile)
{
  int c;

  c = getc(infile);
  if (c == EOF)
    ERREXIT("Premature EOF in JPEG file");
  return c;
}

/* Read 2 bytes, convert to unsigned int */
/* All 2-byte quantities in JPEG markers are MSB first */
static unsigned int
read_2_bytes (FILE *infile)
{
  int c1, c2;

  c1 = getc(infile);
  if (c1 == EOF)
    ERREXIT("Premature EOF in JPEG file");
  c2 = getc(infile);
  if (c2 == EOF)
    ERREXIT("Premature EOF in JPEG file");
  return (((unsigned int) c1) << 8) + ((unsigned int) c2);
}


/*
 * JPEG markers consist of one or more 0xFF bytes, followed by a marker
 * code byte (which is not an FF).  Here are the marker codes of interest
 * in this program.  (See jdmarker.c for a more complete list.)
 */

#define M_SOF0  0xC0		/* Start Of Frame N */
#define M_SOF1  0xC1		/* N indicates which compression process */
#define M_SOF2  0xC2		/* Only SOF0-SOF2 are now in common use */
#define M_SOF3  0xC3
#define M_SOF5  0xC5		/* NB: codes C4 and CC are NOT SOF markers */
#define M_SOF6  0xC6
#define M_SOF7  0xC7
#define M_SOF9  0xC9
#define M_SOF10 0xCA
#define M_SOF11 0xCB
#define M_SOF13 0xCD
#define M_SOF14 0xCE
#define M_SOF15 0xCF
#define M_SOI   0xD8		/* Start Of Image (beginning of datastream) */
#define M_EOI   0xD9		/* End Of Image (end of datastream) */
#define M_SOS   0xDA		/* Start Of Scan (begins compressed data) */
#define M_APP0	0xE0		/* Application-specific marker, type N */
#define M_APP12	0xEC		/* (we don't bother to list all 16 APPn's) */
#define M_COM   0xFE		/* COMment */


/*
 * Find the next JPEG marker and return its marker code.
 * We expect at least one FF byte, possibly more if the compressor used FFs
 * to pad the file.
 * There could also be non-FF garbage between markers.  The treatment of such
 * garbage is unspecified; we choose to skip over it but emit a warning msg.
 * NB: this routine must not be used after seeing SOS marker, since it will
 * not deal correctly with FF/00 sequences in the compressed image data...
 */

static int
next_marker (FILE *infile)
{
  int c;
  int discarded_bytes = 0;

  /* Find 0xFF byte; count and skip any non-FFs. */
  c = read_1_byte(infile);
  while (c != 0xFF) {
    discarded_bytes++;
    c = read_1_byte(infile);
  }
  /* Get marker code byte, swallowing any duplicate FF bytes.  Extra FFs
   * are legal as pad bytes, so don't count them in discarded_bytes.
   */
  do {
    c = read_1_byte(infile);
  } while (c == 0xFF);

  if (discarded_bytes != 0) {
    cpdf_GlobalError(1, "ClibPDF jpegsize", "Warning: garbage data found in JPEG file");
  }

  return c;
}


/*
 * Read the initial marker, which should be SOI.
 * For a JFIF file, the first two bytes of the file should be literally
 * 0xFF M_SOI.  To be more general, we could use next_marker, but if the
 * input file weren't actually JPEG at all, next_marker might read the whole
 * file and then return a misleading error message...
 */

static int
first_marker (FILE *infile)
{
  int c1, c2;

  c1 = getc(infile);
  c2 = getc(infile);
  if (c1 != 0xFF || c2 != M_SOI)
    ERREXIT("Not a JPEG file");
  return c2;
}


/*
 * Most types of marker are followed by a variable-length parameter segment.
 * This routine skips over the parameters for any marker we don't otherwise
 * want to process.
 * Note that we MUST skip the parameter segment explicitly in order not to
 * be fooled by 0xFF bytes that might appear within the parameter segment;
 * such bytes do NOT introduce new markers.
 */

static void
skip_variable (FILE *infile)
/* Skip over an unknown or uninteresting variable-length marker */
{
  unsigned int length;

  /* Get the marker parameter length count */
  length = read_2_bytes(infile);
  /* Length includes itself, so must be at least 2 */
  if (length < 2)
    ERREXIT("Erroneous JPEG marker length");
  length -= 2;
  /* Skip over the remaining bytes */
  while (length > 0) {
    (void) read_1_byte(infile);
    length--;
  }
}


#ifdef DEFINE_ME_TO_ENABLE

/*
 * Process a COM marker.
 * We want to print out the marker contents as legible text;
 * we must guard against non-text junk and varying newline representations.
 */

static void
process_COM (FILE *infile)
{
  unsigned int length;
  int ch;
  int lastch = 0;

  /* Get the marker parameter length count */
  length = read_2_bytes(infile);
  /* Length includes itself, so must be at least 2 */
  if (length < 2)
    ERREXIT("Erroneous JPEG marker length");
  length -= 2;

  while (length > 0) {
    ch = read_1_byte(infile);
    /* Emit the character in a readable form.
     * Nonprintables are converted to \nnn form,
     * while \ is converted to \\.
     * Newlines in CR, CR/LF, or LF form will be printed as one newline.
     */
    if (ch == '\r') {
      printf("\n");
    } else if (ch == '\n') {
      if (lastch != '\r')
	printf("\n");
    } else if (ch == '\\') {
      printf("\\\\");
    } else if (isprint(ch)) {
      putc(ch, stdout);
    } else {
      printf("\\%03o", ch);
    }
    lastch = ch;
    length--;
  }
  printf("\n");
}

#endif		/* DEFINE_ME_TO_ENABLE */

/*
 * Process a SOFn marker.
 * This code is only needed if you want to know the image dimensions...
 */

static void
process_SOFn (FILE *infile, int marker, CPDFimageInfo *jInfo)
{
  unsigned int length;
  unsigned int image_height, image_width;
  int data_precision, num_components;
  int ci;

  length = read_2_bytes(infile);	/* usual parameter length count */

  data_precision = read_1_byte(infile);
  image_height = read_2_bytes(infile);
  image_width = read_2_bytes(infile);
  num_components = read_1_byte(infile);

  jInfo->process = marker;
  jInfo->width = image_width;
  jInfo->height = image_height;
  jInfo->ncomponents = num_components;
  jInfo->bitspersample = data_precision;

/*
  printf("JPEG image is %uw * %uh, %d color components, %d bits per sample\n",
	 image_width, image_height, num_components, data_precision);
  printf("JPEG process: %s\n", process);
*/

  if (length != (unsigned int) (8 + num_components * 3))
    ERREXIT("Bogus SOF marker length");

  for (ci = 0; ci < num_components; ci++) {
    (void) read_1_byte(infile);	/* Component ID code */
    (void) read_1_byte(infile);	/* H, V sampling factors */
    (void) read_1_byte(infile);	/* Quantization table number */
  }
}


/*
 * Parse the marker stream until SOS or EOI is seen;
 * display any COM markers.
 * While the companion program wrjpgcom will always insert COM markers before
 * SOFn, other implementations might not, so we scan to SOS before stopping.
 * If we were only interested in the image dimensions, we would stop at SOFn.
 * (Conversely, if we only cared about COM markers, there would be no need
 * for special code to handle SOFn; we could treat it like other markers.)
 */

static int
scan_JPEG_header (FILE *infile, CPDFimageInfo *jInfo)
{
  int marker;

  /* Expect SOI at start of file */
  if (first_marker(infile) != M_SOI)
	return(-2);

/*    ERREXIT("Expected SOI marker first"); */

  /* Scan miscellaneous markers until we reach SOS. */
  for (;;) {
    marker = next_marker(infile);
    switch (marker) {
      /* Note that marker codes 0xC4, 0xC8, 0xCC are not, and must not be,
       * treated as SOFn.  C4 in particular is actually DHT.
       */
    case M_SOF0:		/* Baseline */
    case M_SOF1:		/* Extended sequential, Huffman */
    case M_SOF2:		/* Progressive, Huffman */
    case M_SOF3:		/* Lossless, Huffman */
    case M_SOF5:		/* Differential sequential, Huffman */
    case M_SOF6:		/* Differential progressive, Huffman */
    case M_SOF7:		/* Differential lossless, Huffman */
    case M_SOF9:		/* Extended sequential, arithmetic */
    case M_SOF10:		/* Progressive, arithmetic */
    case M_SOF11:		/* Lossless, arithmetic */
    case M_SOF13:		/* Differential sequential, arithmetic */
    case M_SOF14:		/* Differential progressive, arithmetic */
    case M_SOF15:		/* Differential lossless, arithmetic */
      process_SOFn(infile, marker, jInfo);
      break;

    case M_SOS:			/* stop before hitting compressed data */
      return marker;

    case M_EOI:			/* in case it's a tables-only JPEG stream */
      return marker;

    case M_COM:
      skip_variable(infile);
      /* process_COM(infile); */	/* disabled */
      break;

    case M_APP12:
      /* Some digital camera makers put useful textual information into
       * APP12 markers, so we print those out too when in -verbose mode.
       */
/*
      if (verbose) {
	printf("APP12 contains:\n");
	process_COM(infile);
      } else
*/
      skip_variable(infile);
      break;

    default:			/* Anything else just gets skipped */
      skip_variable(infile);		/* we assume it has a parameter count... */
      break;
    } /* end switch () */
  } /* end loop */
  return(0);
}

/* Stuff the CPDFimageInfo struct with the JPEG file attributes
   such as image size in pixels, # of color components, bits/component
*/

int read_JPEG_header(const char *filename, CPDFimageInfo *jInfo)
{
int errcode = 0;
FILE *infile;		/* input JPEG file */

    /* Open the input file. */
    if ((infile = fopen(filename, READ_BINARY)) == NULL) {
      	cpdf_GlobalError(1, "ClibPDF read_JPEG_header", "Can't open image file: %s", filename);
      	return(-1);
    }
    /* Scan the JPEG headers. */
    errcode =  scan_JPEG_header(infile, jInfo);
    fclose(infile);
    return(errcode);
}

/* ================= Stuff below this line is for stand-alone testing if this file ============= */
#ifdef JPEG_DEBUG

static const char * progname;	/* program name for error messages */

/*
 * The main program.
 */

int
main (int argc, char **argv)
{
CPDFimageInfo jInfo;
char *procDesc = "?";

#ifdef USE_CCOMMAND
    /* On Mac, fetch a command line. */
    argc = ccommand(&argv);
#endif

    progname = argv[0];
    if (progname == NULL || progname[0] == 0)
    	progname = "jpegsize";	/* in case C library doesn't provide it */

    if(argc != 2) {
	printf("Usage: %s jpegfilename\n", progname);
	exit(0);
    }

    read_JPEG_header( argv[1], &jInfo);

    switch (jInfo.process) {
	case M_SOF0:	procDesc = "Baseline";  break;
	case M_SOF1:	procDesc = "Extended sequential";  break;
	case M_SOF2:	procDesc = "Progressive";  break;
	case M_SOF3:	procDesc = "Lossless";  break;
	case M_SOF5:	procDesc = "Differential sequential";  break;
	case M_SOF6:	procDesc = "Differential progressive";  break;
	case M_SOF7:	procDesc = "Differential lossless";  break;
	case M_SOF9:	procDesc = "Extended sequential, arithmetic coding";  break;
	case M_SOF10:	procDesc = "Progressive, arithmetic coding";  break;
	case M_SOF11:	procDesc = "Lossless, arithmetic coding";  break;
	case M_SOF13:	procDesc = "Differential sequential, arithmetic coding";  break;
	case M_SOF14:	procDesc = "Differential progressive, arithmetic coding"; break;
	case M_SOF15:	procDesc = "Differential lossless, arithmetic coding";  break;
	default:	procDesc = "Unknown";  break;
    }

    /* print out header information */
    printf("Width * Height= %d * %d,  number of components= %d,  bits/sample= %d\n",
	jInfo.width, jInfo.height, jInfo.ncomponents, jInfo.bitspersample);
    printf("JPEG process= %s\n", procDesc);
    /* PDF can take only Baseline JPEG encoding */

    /* All done. */
    exit(EXIT_SUCCESS);
    return 0;			/* suppress no-return-value warnings */
}

void cpdf_GlobalError(int level, const char* module, const char* fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);
        vfprintf(stderr, fmt, ap);
	va_end(ap);
}

#endif		/* JPEG_DEBUG */

