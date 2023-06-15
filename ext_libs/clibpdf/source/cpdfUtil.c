/* cpdfUtil.c -- tiny misc functions.
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

1999-12-04 [io]
     Use below where available:
     struct tm *localtime_r(const time_t * const timep, struct tm *p_tm)

1999-08-17 [IO] 
	## MULTI_THREADING_OK ##

*/

/* Enable below on systems with POSIX thread-safe time functions */
/* #define  THREAD_SAFE_TIME_FUNCS		1 */

/* struct tm *localtime_r(const time_t * const timep, struct tm *p_tm); */

#include "config.h"
#include "version.h"

#include <string.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "cpdflib.h"
#include "cpdf_private.h"
#include "cglobals.h"

#ifdef MacOS8
#include <stat.mac.h>
#include <unix.mac.h>
#include <Memory.h>
#include <Strings.h>	/* c2pstr() */
/* extern size_t strlen(const char *s);
extern int strcmp(const char *s1, const char *s2);
extern void *memcpy(void *s1, const void *s2, size_t n);
*/
#endif

/* ---------------------------------------------------------------------------------------- */

/* The following 3 functions may need to be customized for thread safety
   if the time functions of your C run time library is not multi-threaded.
*/

/* This function should be called after cpdf_open() and before cpdf_init().
   It increments a unique document ID (within a given process with multiple threads of
   PDF document creation).
*/

int _cpdf_inc_docID(CPDFdoc *pdf)
{
int rval;

    /* If pdf->docID is already a positive integer, something else, like:
	cpdf_setDocumentID(int dID); must have already set it, so don't
	mess with it here.
     */
    if(pdf->docID > 0)
	return(pdf->docID);

    /* ##### Mutex Lock before here */
    if(global_docID < 0 || global_docID > 2147483645)		/* 32-bit 'int' max */
	global_docID = 0;
    rval = ++global_docID;
    /* ##### Mutex Unlock here */

    pdf->docID = rval;
    return(rval);
}

/* This is a wrapper around the standard localtime(), which may not be thread safe. */

struct tm *cpdf_localtime(const time_t *timer, struct tm *p_tm)
{
struct tm *rval;

#ifdef THREAD_SAFE_TIME_FUNCS
    rval = localtime_r(timer, struct tm *p_tm);
#else
    rval = localtime(timer);
#endif
    return(rval);
}

/* This is a wrapper around the standard mktime(). */

time_t cpdf_mktime(struct tm *timeptr)
{
time_t rval;

    rval = mktime(timeptr);
    return(rval);
}







/* ---------------------------------------------------------------------------------------- */

/* This allows you set document ID (must be unique within a process) of a given thread
   of PDF creation.  If called at all, this must be called between cpdf_open() and
   cpdf_init().  There is no need to worry about document ID unless you set up
   to use temporary page content files.  The default is to use memory streams.
   If this is not called, cpdf_init() will try to serialize this ID using a
   mutex-locked global variable (see above for _cpdf_inc_docID()).
*/

int cpdf_setDocumentID(CPDFdoc *pdf, int documentID)
{
    if(documentID < 1)
	return(0);
    pdf->docID = documentID;
    return(pdf->docID);
}

/* 3 x 3 matrix for CTM : z0 = z1 = 0, z2 = 1
typedef struct {
    float a; float b;
    float c; float d;
    float x; float y;
} CPDFctm;
*/

/* Will do matrix multiply: T = S x T 
   We use it to give Text CTM concatenation capability, as the "Tm" operator
   does not concatenate CTM, but sets the new one.
*/
void multiplyCTM(CPDFctm *T, const CPDFctm *S)
{
CPDFctm tempM, *T1;
    T1 = &tempM;	/* to avoid confusion, uniformly use -> */
    memcpy(T1, T, (size_t)sizeof(CPDFctm));
    /* using temp copy T1, we will do: T = S x T1
    */
    T->a = S->a * T1->a + S->b * T1->c;
    T->b = S->a * T1->b + S->b * T1->d;
    T->c = S->c * T1->a + S->d * T1->c;
    T->d = S->c * T1->b + S->d * T1->d;
    T->x = S->x * T1->a + S->y * T1->c + T1->x;
    T->y = S->x * T1->b + S->y * T1->d + T1->y;
}


/* Convert a float into M and E, where M * 10 ^ E, and 1 <= M < 10 */
float getMantissaExp(float v, int *iexp)
{
int ie = 0;
float vv = fabs(v);
    if(v == 0.0) {
	*iexp = 0;
	return(0.0);
    }

    while( vv >= 10.0) {
	vv /= 10.0;
	ie++;
    }
    while( vv < 1.0) {
	vv *= 10.0;
	ie--;
    }
    *iexp = ie;
    if(v>=0.0) return(vv);
    else       return(-vv);
}

/* Used in hex to binary conversion function below */
unsigned char _cpdf_nibbleValue(unsigned char hexchar)
{
unsigned char val;
    if(hexchar >= '0' && hexchar <= '9')
	val = hexchar - '0';
    else if(hexchar >= 'A' && hexchar <= 'F')
	val = hexchar - 'A' + 10;
    else if(hexchar >= 'a' && hexchar <= 'f')
	val = hexchar - 'a' + 10;
    else
	val = 0xff;
    return(val);
}

/* binary to HEX conversion function -- for Unicode support */
static char *_cpdf_HexTab= "0123456789ABCDEF00";

/*  datain - binary data of size "length"
    hexout - Hex string representing binary data pointed to by datain.
    length - byte length of binary data pointed to by datain.
    Buffer for "hexout" must be allocated by the caller with a storage for a string
	of size at least 2*length+5 (considering possible addition of "FEFF" and
	string terminator).
*/
char *cpdf_convertBinaryToHex(const unsigned char *datain, char *hexout, long length, int addFEFF)
{
long i;
unsigned char *pin = (unsigned char *)datain;
unsigned short *pcheck = (unsigned short *)datain;
char *pout = hexout;
unsigned char ch;
    if(addFEFF && (*pcheck != 0xfeff && *pcheck != 0xfffe)) {
	*pout++ = 'F';
	*pout++ = 'E';
	*pout++ = 'F';
	*pout++ = 'F';
    }
    for(i=0; i < length; i++) {
	ch = *pin >> 4;
	*pout++ = _cpdf_HexTab[ch];
	ch = *pin++ & 0x0f;
	*pout++ = _cpdf_HexTab[ch];
    }
    *pout = '\0';
    return(hexout);
}

/* After the call, length will contain the byte length of binary data. */
unsigned char *cpdf_convertHexToBinary(const char *hexin, unsigned char *binout, long *length)
{
long bcount=0;
char *pin = (char *)hexin;
unsigned char *pout = binout;
unsigned char ch, chin;
int HighNibble = 1;
    while( (chin = *pin++) != '\0') {
	ch = _cpdf_nibbleValue(chin);
	if(ch > 15)
	    continue;		/* Skip all non-HEX character */

	if(HighNibble) {
	    *pout = ch << 4;
    	    HighNibble = 0;
	}
	else {
	    *pout++ |= ch;
    	    HighNibble = 1;
	    bcount++;
	}
    }
    *length = bcount;	/* byte count */
    return(binout);
}

void cpdf_hexStringMode(CPDFdoc *pdf, int flag)
{
    pdf->hex_string = flag;
}

void cpdf_useContentMemStream(CPDFdoc *pdf, int flag)
{
#ifdef EXT_ZLIBCOMP
	pdf->useContentMemStream = 0;	/* can't use external command if memory stream is used */
#else
	pdf->useContentMemStream = flag;
#endif

}


void cpdf_useStdout(CPDFdoc *pdf, int flag)
{
    pdf->useStandardOutput = flag;	/* send output to stdout if non-zero */
}

/* This function will set /Creator entry in the Info dictionary of the PDF file. */

void cpdf_setCreator(CPDFdoc *pdf, char *pname)
{
    strncpy(pdf->creator_name, pname, 62);
}


void cpdf_setTitle(CPDFdoc *pdf, char *pname)
{
    strncpy(pdf->file_title, pname, 62);
}



void cpdf_setSubject(CPDFdoc *pdf, char *pname)
{
    strncpy(pdf->file_subject, pname, 62);
}


void cpdf_setKeywords(CPDFdoc *pdf, char *pname)
{
    strncpy(pdf->file_keywords, pname, 126);
}



void cpdf_setPDFLevel(CPDFdoc *pdf, int major, int minor)
{
    if(pdf->pdfLevelMin < minor || pdf->pdfLevelMaj < major)
	pdf->pdfLevelMin = minor;
    if(pdf->pdfLevelMaj < major)	/* PDF level only gets bumped up as needed */
	pdf->pdfLevelMaj = major;	/* Do not allow going back to lower level. */
}

/* Other comments that can go anywhere in PS stream */

int cpdf_comments(CPDFdoc *pdf, char *comments) 
{
	if(comments != NULL) {
	    if(pdf->useContentMemStream)
		cpdf_writeMemoryStream(pdf->currentMemStream, comments, strlen(comments));
	    else
	        fprintf(pdf->fpcontent,"%s", comments);
	}
	return(0);
}

/* User must copy the string to retain it */

char *cpdf_getOutputFilename(CPDFdoc *pdf)
{
   return pdf->filenamepath;
}


/* this must be called before each call to cpdf_init() */

void cpdf_setOutputFilename(CPDFdoc *pdf, const char *file)
{
    strncpy(pdf->filenamepath, file, 1022);
    pdf->filename_set  = 1;
    pdf->usePDFMemStream = 0;	/* do not use memory stream for PDF file itself */
}


void cpdf_setDefaultDomainUnit(CPDFdoc *pdf, float defunit)
{
    pdf->defdomain_unit = defunit;
}


char *cpdf_version(void)
{
    return(CPDF_VERSION);
}

char *cpdf_platform(void)
{
    return(PLATFORM_NAME);
}

/* Convert Unix pathname to pathname for other platforms with dirseparator specified. */
char *cpdf_convertUpathToOS(char *pathbuf, char *upath)
{
char *ip, *op;
    ip = upath;
    op = pathbuf;
    while( *ip != '\0') {
	if(*ip == '/')
	    *op++ = DIR_SEPARATOR;
	else
	    *op++ = *ip;		/* copy */
	ip++;
    }
    *op = '\0';
    return(pathbuf);
}


CPDFerrorHandler cpdf_setErrorHandler(CPDFdoc *pdf, CPDFerrorHandler handler)
{
    CPDFerrorHandler old = pdf->_cpdfErrorHandler;
    pdf->_cpdfErrorHandler = handler;
    return (old);
}



void cpdf_Error(CPDFdoc *pdf, int level, const char* module, const char* fmt, ...)
{
    if(pdf->_cpdfErrorHandler) {
	va_list ap;
	va_start(ap, fmt);
	(*(pdf->_cpdfErrorHandler))(pdf, level, module, fmt, ap);
	va_end(ap);
    }
}


void _cpdf_defaultErrorHandler(CPDFdoc *pdf, int level, const char* module, const char* fmt, va_list ap)
{
    if(module != NULL)
	fprintf(stderr, "%s: ", module);
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, ".\n");
    if(level < 0)
	exit(level);
}


/* ===== Global Error Handler ======================================== */
/* This error handler is shared by all threads and meant for critical errors
   that probably require termination of the entire process.
*/

CPDFglobalErrorHandler cpdf_setGlobalErrorHandler(CPDFglobalErrorHandler handler)
{
    CPDFglobalErrorHandler old = _cpdfGlobalErrorHandler;
    _cpdfGlobalErrorHandler = handler;
    return (old);
}

void cpdf_GlobalError(int level, const char* module, const char* fmt, ...)
{
    if(_cpdfGlobalErrorHandler) {
	va_list ap;
	va_start(ap, fmt);
	(*_cpdfGlobalErrorHandler)(level, module, fmt, ap);
	va_end(ap);
    }
}

void _cpdf_defaultGlobalErrorHandler(int level, const char* module, const char* fmt, va_list ap)
{
    if(module != NULL)
	fprintf(stderr, "%s: ", module);
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, ".\n");
    if(level < 0)
	exit(level);
}


/* Big thanks for Rodney Sitton for the code of this function. */
int _Big_Endian(void)
{
    union be{
	    char    byte[2];
	    short   word;
    } be;

	be.word = 0x00FF;
	return((int)(be.byte[0] == 0x00));
}


/* private functions --------------------------------------------------------------------- */
void str_append_int(char *buf, int num)
{
char tnumbuf[256];
	sprintf(tnumbuf, "%d", num);
	strncat(buf, tnumbuf, 250);
}

void _cpdf_malloc_check(void *buf)
{
    if(!buf) {
	cpdf_GlobalError(-2, "ClibPDF malloc", "Memory allocation failure");
    }
}

