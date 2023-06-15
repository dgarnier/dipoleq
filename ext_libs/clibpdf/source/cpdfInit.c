/* cpdfInit.c  -- Initialization and Primary Processing Module ----------------------------
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

1999-08-21 [IO] 
	## MULTI_THREADING_OK ##
1999-08-12 [IO]
	-- Major rework for multi-threaded version started.
1998-06-30 [IO]
	-- First version, NEXTSTEP3.x
*/

#include "config.h"
#include "version.h"

#define	    COMPRESSION_CCITTFAX4	4	/* CCITT Group 4 fax encoding in TIFF */

#include <string.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef UNIX
#include <unistd.h>		/* for access() etc */
#include <sys/stat.h>
#include <sys/types.h>
#include <pwd.h>
#endif

#if defined(NEXTSTEP)
	typedef int	pid_t;
#endif

#if defined(_WIN32) || defined(WIN32)
/* #include <sys/file.h> */
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif


#ifdef MacOS8
#include <stat.mac.h>
#include <unix.mac.h>
#include <StandardFile.h>
#include <Memory.h>
#include <Files.h>
#include <Strings.h>	/* c2pstr() */
/* extern void *memcpy(void *s1, const void *s2, size_t n); */
#endif

#ifndef EXT_ZLIBCOMP
#include "zlib.h"
#endif

#define MAINDEF 1

#include "cpdflib.h"		/* This must be included before all other local include files */
#include "cpdf_private.h"
#include "cglobals.h"

extern short cpdf_encode2WinAnsi[][256];	/* in cpdfAFM.h */

/* VC++ malloc debug: run the executable from VC++ workspace: Build->Debug->Go */

#if defined(_WIN32) || defined(WIN32)
#ifdef _DEBUG
#ifndef __BORLANDC__
#include <crtdbg.h>
#endif

static init_debug_Done = 0;		/* We should keep this, since this is GLOBAL to all threads */
#endif

void my_init_debug() {
#ifdef _DEBUG
    if(init_debug_Done == 0) {
	/* do this only once */
        int tmpFlag = _CrtSetDbgFlag( _CRTDBG_REPORT_FLAG);
        tmpFlag |= _CRTDBG_CHECK_ALWAYS_DF; /* check memory at every alloc & dealloc */
        tmpFlag |= _CRTDBG_LEAK_CHECK_DF; /* checks for leaks at end of program */
        _CrtSetDbgFlag(tmpFlag);
	init_debug_Done = 1;
    }
#endif
}
#endif

/* ---------------------------------------------------------------------------------------- */
/* Key puclic API functions for ClibPDF library */
/* ---------------------------------------------------------------------------------------- */

/* This sets and works GLOBALLY.  Changes default for all threads. */
/* If this function is used, it must be called before cpdf_open(). */
void cpdf_setGlobalDocumentLimits(int maxPages, int maxFonts, int maxImages, int maxAnnots, int maxObjects)
{
    if(maxPages > 0)
	g_nMaxPages = maxPages;
    if(maxFonts > 0)
	g_nMaxFonts = maxFonts;	/* FIXMED -- malloc size bug elsewhere */
    if(maxImages > 0)
	g_nMaxImages = maxImages;
    if(maxAnnots > 0)
	g_nMaxAnnots = maxAnnots;
    if(maxObjects > 0)
	g_nMaxObjects = maxObjects;
}


/* #define  CPDF_DEBUG	1 */

/* New multi-threaded open takes document limits as a struct, and returns a document context. */
CPDFdoc *cpdf_open(int pspdf, CPDFdocLimits *docLimits)
{
int i;
CPDFdoc *nPDF;

#if defined(_WIN32) || defined(WIN32)
	my_init_debug();
#endif

	nPDF = (CPDFdoc *)malloc((size_t)sizeof(CPDFdoc));	/* allocate memory for PDF doc struct */
	_cpdf_malloc_check((void *)nPDF);

	if(docLimits == NULL) {
	    nPDF->nMaxPages = g_nMaxPages;		/* if NULL, take from default values */
	    nPDF->nMaxFonts = g_nMaxFonts;
	    nPDF->nMaxImages = g_nMaxImages;
	    nPDF->nMaxAnnots = g_nMaxAnnots;
	    nPDF->nMaxObjects = g_nMaxObjects;
	}
	else {
	    if(docLimits->nMaxPages > 0)
		nPDF->nMaxPages = docLimits->nMaxPages;  /* copy from struct */
	    else
		nPDF->nMaxPages = g_nMaxPages;
	    if(docLimits->nMaxFonts > 0)
		nPDF->nMaxFonts = docLimits->nMaxFonts;	/* FIXED */
	    else
	        nPDF->nMaxFonts = g_nMaxFonts;
	    if(docLimits->nMaxImages > 0)
		nPDF->nMaxImages = docLimits->nMaxImages;
	    else
		nPDF->nMaxImages = g_nMaxImages;
	    if(docLimits->nMaxAnnots > 0)
		nPDF->nMaxAnnots = docLimits->nMaxAnnots;
	    else
		nPDF->nMaxAnnots = g_nMaxAnnots;
	    if(docLimits->nMaxObjects)
		nPDF->nMaxObjects = docLimits->nMaxObjects;
	    else
		nPDF->nMaxObjects = g_nMaxObjects;
	}
	/* Allocate arrays for various PDF objects */
	nPDF->pageInfos = (CPDFpageInfo *)malloc((size_t)((nPDF->nMaxPages+1)*sizeof(CPDFpageInfo)));
	_cpdf_malloc_check((void *)nPDF->pageInfos);
	nPDF->kidsIndex = (int *)malloc((size_t)((nPDF->nMaxPages+1)*sizeof(int)));
	_cpdf_malloc_check((void *)nPDF->kidsIndex);
	nPDF->fontInfos = (CPDFfontInfo *)malloc((size_t)(nPDF->nMaxFonts*sizeof(CPDFfontInfo)));
	_cpdf_malloc_check((void *)nPDF->fontInfos);
	nPDF->imageInfos = (CPDFimageInfo *)malloc((size_t)(nPDF->nMaxImages*sizeof(CPDFimageInfo)));
	_cpdf_malloc_check((void *)nPDF->imageInfos);
	nPDF->annotInfos = (CPDFannotInfo *)malloc((size_t)(nPDF->nMaxAnnots*sizeof(CPDFannotInfo)));
	_cpdf_malloc_check((void *)nPDF->annotInfos);
	nPDF->objByteOffset = (long *)malloc((size_t)(nPDF->nMaxObjects*sizeof(long)));
	_cpdf_malloc_check((void *)nPDF->objByteOffset);
	nPDF->objIndex = (int *)malloc((size_t)(nPDF->nMaxObjects*sizeof(int)));
	_cpdf_malloc_check((void *)nPDF->objIndex);

	nPDF->extFontList = (CPDFextFont *)malloc((size_t)(nPDF->nMaxFonts*sizeof(CPDFextFont)));
	_cpdf_malloc_check((void *)nPDF->extFontList);
#ifdef CPDF_DEBUG
	fprintf(stderr, "sizeof(CPDFpageInfo)= %ld\n", sizeof(CPDFpageInfo));
	fprintf(stderr, "sizeof(CPDFfontInfo)= %ld\n", sizeof(CPDFfontInfo));
	fprintf(stderr, "sizeof(CPDFimageInfo)= %ld\n", sizeof(CPDFimageInfo));
	fprintf(stderr, "sizeof(CPDFannotInfo)= %ld\n", sizeof(CPDFannotInfo));
	fprintf(stderr, "sizeof(CPDFextFont)= %ld\n", sizeof(CPDFextFont));
	fprintf(stderr, "sizeof(fontInfos)= %ld\n", sizeof(fontInfos));

#endif
#undef CPDF_DEBUG
	/* Initialize global data for document */
	_cpdf_initDocumentStruct(nPDF);
	nPDF->ps_pdf_mode = pspdf;
	cpdf_setCreator(nPDF, "A ClibPDF program");
	cpdf_setTitle(nPDF, "No Title");
	cpdf_setSubject(nPDF, "None");
	cpdf_setKeywords(nPDF, "ClibPDF");

	cpdf_setFontMapFile(nPDF, DEF_FONTMAPFILE);
	cpdf_setFontDirectories(nPDF, DEF_PFM_DIR, DEF_PFB_DIR);

	/* set default compression command: zlibcomp */
	cpdf_setCompressionFilter(nPDF, ZLIB_COMPRESS_PATH, "/FlateDecode");
	nPDF->monthName = (char **)malloc((size_t)(12*sizeof(char *)));
	for(i=0; i<12; i++)
	    nPDF->monthName[i] = NULL;
	cpdf_setMonthNames(nPDF, monthNameEnglish);
	cpdf_setErrorHandler(nPDF, _cpdf_defaultErrorHandler);	/* setup default error handler */

	/* provide non-overflowing memory stream for scratch pad use */
	nPDF->scratchMem = cpdf_openMemoryStream();
#if defined(_WIN32) || defined(WIN32)
    #ifdef __BORLANDC__
        _fmode = O_BINARY;	/* FIXME: -- what to do for multi-threading, though this is probably OK. */
    #else
    	_fmode = _O_BINARY;	/* Why is this so hard in _WIN32 ?? */
    #endif
#endif
#ifdef EXT_ZLIBCOMP
	nPDF->useContentMemStream = 0;    /* can't use external command if memory stream is used */
#endif
	nPDF->inlineImages = 0;
	nPDF->init_check = 1;
	return(nPDF);
}


void cpdf_enableCompression(CPDFdoc *pdf, int cmpON)
{
    if(cmpON == 0) {
	pdf->compressionON = 0;
	return;
    }
/* Check only if we are using external compression command */
#ifdef EXT_ZLIBCOMP
#ifndef MacOS8			/* why Metrowerks CW IDE 1.7.4 doesn't have this */
    else if( access(pdf->compress_command, F_OK) != 0) {
	/* Check if the compression program specifiled in compress_command exists.
	   If the compress_command does not exist, we will fall back on NO compression.
	   Things should still work.  This check provides some safety.
	*/
	pdf->compressionON = 0;	/* compression program does not exist, so don't use it */
	return;
    }
#endif
#endif
    /* OK, we want compression, and the compression program exists */
    /* Let's see what kind of decode filter we specify */
    if(pdf->streamFilterList != NULL) {
	pdf->compressionON = cmpON;		/* Only if filter command is set */
	if(cmpON && (strstr(pdf->streamFilterList, "FlateDecode") != NULL))
	    cpdf_setPDFLevel(pdf, 1, 2);	/* FlateDecode requires PDF-1.2 */
    }
    else {
	pdf->compressionON = 0;		/* compression filter has not been set */
    }
}


/* "Command" must be a FULL pathname of the executable that implements Acrobat Reader-
   compatible compression, of the type decodeded by "filterlist" decode-filter chain.
   Existence of the compression command file is checked.  If it does not exist,
   no compression will be used, but things should still work.
   Default is: cpdf_setCompressionFilter(ZLIB_COMPRESS_PATH, "/FlateDecode");
   where ZLIB_COMPRESS_PATH is specified in config.h.
*/

void cpdf_setCompressionFilter(CPDFdoc *pdf, const char *command, const char *filterlist)
{
	if(pdf->compress_command) free(pdf->compress_command);
	pdf->compress_command = (char *)malloc((size_t)(strlen(command) + 1));
	_cpdf_malloc_check((void *)pdf->compress_command);
	strcpy(pdf->compress_command, command);
	/* Decode filter list as specified in PDF /Filter [] list */
	if(pdf->streamFilterList) free(pdf->streamFilterList);
	pdf->streamFilterList = (char *)malloc((size_t)(strlen(filterlist) +1));
	_cpdf_malloc_check((void *)pdf->streamFilterList);
	strcpy(pdf->streamFilterList, filterlist);
}



void cpdf_init(CPDFdoc *pdf)
{
int i;
CPDFpageInfo *pT;
#ifdef UNIX
uid_t myuid;
struct passwd *pw;
extern uid_t getuid(void);
	myuid = getuid();
	pw = getpwuid(myuid);
	strncpy(pdf->username, pw->pw_gecos, 32);
	strcat(pdf->username, " [");
	strncat(pdf->username, pw->pw_name, 16);
	strcat(pdf->username, "]");

	/* set docID (sort of thread ID) in case temp files are used.
	   This should not matter if PDF is generated entirely in memory.
	 */
	_cpdf_inc_docID(pdf);

#endif
#ifdef MacOS8
	strncpy(pdf->username, "MacOS 8.* or less User", 60);
#endif
#if defined(_WIN32) || defined(WIN32)
	strncpy(pdf->username, "Windows NT/95/98 User", 60);
#endif

	if(pdf->init_check != 1) {
	    cpdf_Error(pdf, -1, "ClibPDF", "cpdf_open() has not been called.");
	}
	pdf->init_check = 2;

	/* For later sanity check to detect undefined pages, e.g.,
	   cpdf_pageInit() may be called with page number skipped.  If that
	   happens, we skip undefined pages when writing the final file.
	   pageInfos[0] is unused.  Initialize it here, but don't free.
	*/
	for(i=0; i<= pdf->nMaxPages; i++) {
	    pT = &pdf->pageInfos[i];
	    pT->pagenum = -1;
	    pT->npFont = 0;
	    pT->npImage = 0;
	    pT->npAnnot = 0;
	    pT->status = 0;			/* bit flags: bit-1=compressed, bit-0=closed */

	    pT->length_compressed = 0;
	    pT->compressedStream = NULL;
	    pT->pageMemStream = NULL;	/* initialize all malloc items to NULL */
	    pT->defDomain = NULL;
	    pT->fontIdx = NULL;
	    pT->imageIdx = NULL;
	    pT->annotIdx = NULL;
	    pT->mediaBox = NULL;
	    pT->cropBox = NULL;
	    pT->fppage = NULL;
	    pT->contentfile = NULL;
	    pT->duration = -1.0;
	    pT->transition = NULL;
	}

	/* If we use memory stream for the entire PDF (file), then creat that too. */
	if(pdf->usePDFMemStream)
	    pdf->pdfMemStream = cpdf_openMemoryStream();	/* write PDF to memory stream */
	else
	    _cpdf_file_open(pdf);				/* open a stream to file for PDF output */
}


/* cpdf_pageInit() must be called for each page.
   This creats a memory stream and sets it as currentStream.
*/

int cpdf_pageInit(CPDFdoc *pdf, int pagenum, int rot, const char *mediaboxstr, const char *cropboxstr)
{	
float llxt, llyt, urxt, uryt;
CPDFpageInfo  *nP;		/* pointer to page info struct */
/* CPDFmemStream *contentStream; */	/* memory stream */
#if defined(UNIX) && !defined(SunOS5x)
extern pid_t getpid(void);
int mypid = getpid();
#else
int mypid = 123;
#endif
	if(pagenum < 1) {
	    cpdf_Error(pdf, 1, "ClibPDF", "Page number must be 1 or greater");
	    return(-1);
	}
	else if(pagenum >= pdf->nMaxPages) {
	    cpdf_Error(pdf, 1, "ClibPDF", "Too many pages: %d. Increase limit by cpdf_open()", pagenum);
	    return(-1);
	}
	nP = &pdf->pageInfos[pagenum];
	if(nP->pagenum != -1) {
	    cpdf_Error(pdf, 1, "ClibPDF", "Page %d is already initialized", pagenum);
	    cpdf_setCurrentPage(pdf, pagenum);
	    return(0);
	}
	nP->pagenum = pagenum;
	pdf->currentPage = pagenum;
	if(pdf->numPages < pagenum) pdf->numPages = pagenum;

	if(rot == ALTLANDSCAPE) 
	    pdf->display_rotation = 90;		/* rot=2: landscape means rotation for display */
	else if(rot == LANDSCAPE)
	    pdf->display_rotation = 270;	/* rot=1: landscape means rotation for display */
	else
	    pdf->display_rotation = 0;		/* rot=0: portrait means no rotation for display */
	nP->orientation = pdf->display_rotation;

	/* nP->mediaBox, and nP->cropBox are given malloced storage for strings in cpdf_setPageSize() */
	cpdf_setPageSize(pdf, mediaboxstr, cropboxstr);
	/* create and set default domain in inch coordinate system */
	/* get default domain size from cropBox string (defined in %%BoundingBox format) in cglobals.h */
	sscanf(pdf->cropBox, "%f %f %f %f", &llxt, &llyt, &urxt, &uryt);

	/* Linked list of font, image, annotation for this page.  Stores index into xxxxInfos[] */
	nP->fontIdx = (CPDFintList *)malloc((size_t)sizeof(CPDFintList));
	_cpdf_malloc_check((void *)nP->fontIdx);
	nP->fontIdx->value = -1;
	nP->fontIdx->next = NULL;
	
	nP->imageIdx = (CPDFintList *)malloc((size_t)sizeof(CPDFintList));
	_cpdf_malloc_check((void *)nP->imageIdx);
	nP->imageIdx->value = -1;
	nP->imageIdx->next = NULL;

	nP->annotIdx = (CPDFintList *)malloc((size_t)sizeof(CPDFintList));
	_cpdf_malloc_check((void *)nP->annotIdx);
	nP->annotIdx->value = -1;
	nP->annotIdx->next = NULL;

	/* Each page has its own domains */
	/* (CPDFplotDomain *)cpdf_createPlotDomain(float x, float y, float w, float h,
				float xL, float xH, float yL, float yH,
				int xtype, int ytype, int rsvd) */
	pdf->defaultDomain = cpdf_createPlotDomain(llxt, llyt, urxt-llxt, uryt-llyt,
				llxt/pdf->defdomain_unit, (urxt-llxt)/pdf->defdomain_unit,
				llyt/pdf->defdomain_unit, (uryt-llyt)/pdf->defdomain_unit, LINEAR, LINEAR, 0);
	nP->defDomain = pdf->defaultDomain;

	/* These default settings probably should be moved into createDomain itself */
	cpdf_setLinearMeshParams(pdf->defaultDomain, Y_MESH, 0.0, 1.0, 0.0, 0.5);
	cpdf_setLinearMeshParams(pdf->defaultDomain, X_MESH, 0.0, 1.0, 0.0, 0.5);
	cpdf_setPlotDomain(pdf, pdf->defaultDomain);
        _cpdf_resetTextCTM(pdf);

	/* If we use memory stream for Content, then create stream and set is as current stream. */
	if(pdf->useContentMemStream) {
	    nP->pageMemStream = cpdf_openMemoryStream();
	    cpdf_setCurrentMemoryStream(pdf, nP->pageMemStream);
	}
	/* === Now, set filename for temporary file or memory stream to store Content === */
	else {	/* use temporary file for Content */
	    strcpy(pdf->contentfile, TEMP_DIRECTORY);
	    strcat(pdf->contentfile, "_cpdf");
	    str_append_int(pdf->contentfile, mypid);
	    strcat(pdf->contentfile, "-");
	    str_append_int(pdf->contentfile, pdf->docID);	/* give different thread a diffent number */
	    strcat(pdf->contentfile, "-");
	    str_append_int(pdf->contentfile, pagenum);
	    pdf->fpcontent = fopen(pdf->contentfile, BINARY_WRITE);	/* open for actual drawing code */
	    if(pdf->fpcontent == NULL) {
		    cpdf_Error(pdf, 2, "ClibPDF", "Cannot open %s as temporary content file", pdf->contentfile);
		    return(-2);
	    }

	    nP->fppage = pdf->fpcontent;
	    nP->contentfile = (char *)malloc((size_t)(strlen(pdf->contentfile) + 8));  /* ".zlib" */
	    _cpdf_malloc_check((void *)nP->contentfile);
	    strcpy(nP->contentfile, pdf->contentfile);
	}
	return(0);
}



void cpdf_setPageDuration(CPDFdoc *pdf, float seconds)
{
	pdf->pageInfos[pdf->currentPage].duration = seconds;
}


int cpdf_setPageTransition(CPDFdoc *pdf, int type, float duration, float direction, int HV, int IO)
{
char trbuf[128];
char *hv, *io;
    if(HV) hv = "/H";		/* /Dm value */
    else   hv = "/V";

    if(IO) io = "/I";		/* /M value */
    else   io = "/O";

    switch(type) {
	default:
	case TRANS_NONE:
			sprintf(trbuf, "/S /R");
			break;
	case TRANS_SPLIT:
			sprintf(trbuf, "/S /Split /D %.3f /Dm %s /M %s", duration, hv, io);
			break;
	case TRANS_BLINDS:
			sprintf(trbuf, "/S /Blinds /D %.3f /Dm %s", duration, hv);
			break;
	case TRANS_BOX:
			sprintf(trbuf, "/S /Box /D %.3f /M %s", duration, io);
			break;
	case TRANS_WIPE:
			sprintf(trbuf, "/S /Wipe /D %.3f /Di %.2f", duration, direction);
			break;
	case TRANS_DISSOLVE:
			sprintf(trbuf, "/S /Dissolve /D %.3f", duration);
			break;
	case TRANS_GLITTER:
			sprintf(trbuf, "/S /Glitter /D %.3f /Di %.2f", duration, direction);
			break;
    }
    /* now store transition spec string in pageInfos[] */
    if(pdf->pageInfos[pdf->currentPage].transition)
	free(pdf->pageInfos[pdf->currentPage].transition);
    pdf->pageInfos[pdf->currentPage].transition = (char *)malloc((size_t)(strlen(trbuf) +1));
    _cpdf_malloc_check((void *)(pdf->pageInfos[pdf->currentPage].transition));
    strcpy(pdf->pageInfos[pdf->currentPage].transition, trbuf);
    return(0);
}


/* Switch current page to another previously initialized page.
   With this function, you can draw into multiple pages in an interleaved manner.
   Switching pages will reset the current domain to the default domain for that page.
*/

int cpdf_setCurrentPage(CPDFdoc *pdf, int page)
{
CPDFpageInfo *cP;

	if(page == pdf->currentPage) return(0);	/* no change */
	cP = &pdf->pageInfos[page];
	if( cP->pagenum == -1) {
	    cpdf_Error(pdf, 2, "ClibPDF", "Cannot switch to page %d because it has not been initialized", page);
	    return(-1);
	}
	else if( cP->status & 1) {
	    cpdf_Error(pdf, 2, "ClibPDF", "Cannot switch to page %d because it has already been finalized", page);
	    return(-2);
	}
	/* Now, we should be able to switch. */
	pdf->currentPage = page;
	pdf->defaultDomain = cP->defDomain;			/* switch default domain */
	cpdf_setPlotDomain(pdf, pdf->defaultDomain);

	if(pdf->useContentMemStream) {
	    cpdf_setCurrentMemoryStream(pdf, cP->pageMemStream);
	}
	else {
	    pdf->fpcontent = cP->fppage;			/* switch content file descriptor */
	    strcpy(pdf->contentfile, cP->contentfile);		/* sync current content filename */
	}

	/* Also copy mediaBox, cropBox */
	strncpy(pdf->mediaBox, cP->mediaBox, 62);
	strncpy(pdf->cropBox, cP->cropBox, 62);

	return(0);
}

/* This function closes the content stream file or memory stream for a given page.
   The purpose of this function is to close and finalize completed pages so
   that memory streams for completed pages to be compressed.  This will
   conserve memory usage when generating a PDF file with a large number of pages.
   If you are generating PDF files with a few pages, this function does not make
   much difference.
   Compression will take place within this call only if content memory streams are
   used, but not when temp files are used (assumes there is no problem with disk space).
   If temp files are used for page content, this function simply closes the FILE pointer
   associated with the page.
   (This has an effect even when using content memory stream starting with
	version 1.01c.)
*/

void cpdf_finalizePage(CPDFdoc *pdf, int page)
{
    /* only if page is not empty and it has not yet been closed */
    if((pdf->pageInfos[page].pagenum != -1) && ((pdf->pageInfos[page].status & 3) == 0)) {
	if(pdf->useContentMemStream)
	    _cpdf_closeContentMemStreamForPage(pdf, page);
	else
	    fclose(pdf->pageInfos[page].fppage);
	pdf->pageInfos[page].status |= 1;		/* set status bit-0: page closed flag */
    }
}


int cpdf_savePDFmemoryStreamToFile(CPDFdoc *pdf, const char *file)
{
int retcode= -1;
#ifdef MacOS8
    OSType filetype='PDF ', creator='CARO';	/* Acrobat/PDF file */
    /* The convention of storing four characters in a long, and of specifying
    them with single quotes is not part of standard C, but is a standard
    Macintosh extension to C. */
#endif

    strncpy(pdf->filenamepath, file, 1022);
    pdf->filename_set  = 1;
    if(pdf->usePDFMemStream && pdf->pdfMemStream)
	retcode = cpdf_saveMemoryStreamToFile(pdf->pdfMemStream, file);

#ifdef MacOS8
    SetFileInfo(file, filetype, creator);	/* set PDF file type */
#endif
    return(retcode);
}


char *cpdf_getBufferForPDF(CPDFdoc *pdf, int *length)
{
int bufsize;
char *mbuff;
    if(pdf->usePDFMemStream && pdf->pdfMemStream) {
        cpdf_getMemoryBuffer(pdf->pdfMemStream, &mbuff, length, &bufsize);
    }
    else {
	mbuff = NULL;
	*length = 0;
    }
    return mbuff;
}

void _cpdf_pdfWrite(CPDFdoc *pdf, char *s)
{
int bcount = strlen(s);
    if(pdf->usePDFMemStream)
	cpdf_writeMemoryStream(pdf->pdfMemStream, s, bcount);
    else
	fputs(s, pdf->fpcg);
    pdf->currentByteCount += bcount;
}

/* Mainly for writing Unicode string data */
void _cpdf_pdfBinaryWrite(CPDFdoc *pdf, char *s, long length)
{
    if(pdf->usePDFMemStream)
	cpdf_writeMemoryStream(pdf->pdfMemStream, s, length);
    else
	fwrite((void *)s, 1, length, pdf->fpcg);
    pdf->currentByteCount += length;
}


/* Close this PDF document context */

void cpdf_close(CPDFdoc *pdf)
{
    /* The following 2 lines is duplicated in cpdf_finalizeAll() which should be the
	place where the file is closed in normal processing.  It is here also just in case
	the program flow jumps here due to exceptions, so there is no leak in file
	descriptors.
    */
    if(!(pdf->usePDFMemStream))
	_cpdf_file_close(pdf);		/* close the PDF file (duplicated here in case of error) */
				        /* but keep PDF memory stream open if that is used */

    /* The stuff below used to be in cpdf_finalizeAll(), but moved here to make error handling easier. */
    _cpdf_freeAllFontInfos(pdf);	/* release memory for malloc'ed CPDFfontInfo's */
    _cpdf_freeAllImageInfos(pdf);	/* release memory for malloc'ed CPDFimageInfo's */
    _cpdf_freeAllAnnotInfos(pdf);
    if(pdf->numOutlineEntries)
	_cpdf_freeAllOutlineEntries(pdf->firstOLentry, pdf->lastOLentry);
    _cpdf_freeAllPageInfos(pdf);

    if(pdf->compress_command) free(pdf->compress_command);
    pdf->compress_command = NULL;
    if(pdf->streamFilterList) free(pdf->streamFilterList);
    pdf->streamFilterList = NULL;
    pdf->numAnnots = 0;		/* reset annotation for next time */
    /* cpdf_freePlotDomain(pdf, pdf->defaultDomain); */	/* now done in _cpdf_freeAllPageInfos() */
    pdf->defaultDomain = NULL;
    pdf->currentDomain = NULL;
    pdf->init_check = 0;
    /* Previously (until 2.00-beta1-0), the section above was in cpdf_finalizeAll() */

	cpdf_closeMemoryStream(pdf->scratchMem);	/* release scratch pad buffer */
	if(pdf->usePDFMemStream && pdf->pdfMemStream)
	    cpdf_closeMemoryStream(pdf->pdfMemStream);	/* release memory stream for PDF output */
	_cpdf_freeMonthNames(pdf);
	free(pdf->monthName);
	/* free arrays for document objects */
	free(pdf->pageInfos);
	free(pdf->kidsIndex);
	free(pdf->fontInfos);
	free(pdf->imageInfos);
	free(pdf->annotInfos);
	free(pdf->objByteOffset);
	free(pdf->objIndex);
	free(pdf->extFontList);
	free(pdf);		/* free the PDf document struct */
}

/* cpdf_finalizeAll() is where the PDF file is actually written.
   The drawing code (Contents) is first written to a temp file
   or memory stream before this function is called.
   The size of the temp file is determined here and copied to the PDF file.
   If the PDF stream is to a file, the file is closed in this function.
   Put some binary data near the beginning so binary/text autodetection
   in some software will treat PDF files as binary.
   Do not modify the bmagic string below.  */
char *_cpdf_bmagic ="zG_\325\371\337J\244\030\267\260#s6\037\246dR L\204s\037";

void cpdf_finalizeAll(CPDFdoc *pdf)
{
int i, objectCount;
#ifdef MacOS8
    OSType filetype='PDF ', creator='CARO';	/* Acrobat/PDF file */
    /* The convention of storing four characters in a long, and of specifying
    them with single quotes is not part of standard C, but is a standard
    Macintosh extension to C. */
#endif
    if(pdf->linearizeON) {
	/* === Generate linearized PDF file === */
#ifdef USE_LINEARIZER
	_cpdf_finalizeLinearAll(pdf);
#endif
    }
    else {
	/* === Do not linearize PDF file === */
	pdf->currentByteCount = 0;
	sprintf(pdf->spbuf, "%%PDF-%d.%d\n%s\n", pdf->pdfLevelMaj, pdf->pdfLevelMin, _cpdf_bmagic);
	_cpdf_pdfWrite(pdf, pdf->spbuf);
	pdf->objByteOffset[0] = 0;
	pdf->objByteOffset[1] = pdf->currentByteCount;

	/* --- SERIALIZE OBJECTS ----------------------------------------------
	   Scan through all objects (except for xref and trailer) in order, and
	   assign serialized object number for xref entries.
	*/
	objectCount = 1;
	pdf->objIndex[CPDF_Catalog] = objectCount++;
	pdf->objIndex[CPDF_Outlines] = objectCount++;
	pdf->objIndex[CPDF_Pages] = objectCount++;
	pdf->objIndex[CPDF_ProcSet] = objectCount++;
/*
	pdf->objIndex[CPDF_Page] = objectCount++;
	pdf->objIndex[CPDF_Contents] = objectCount++;
*/
	for(i=1; i<= pdf->numPages; i++) {
	    /* only those pages with pagenum != -1 are valid */
	    if(pdf->pageInfos[i].pagenum != -1) {
		pdf->pageInfos[i].objIndex = objectCount++;	/* Page object for this page */
		pdf->kidsIndex[pdf->numKids++] = pdf->pageInfos[i].objIndex;
		pdf->pageInfos[i].parent = pdf->objIndex[CPDF_Pages];	/* Only one Pages object (see above) */
		/* Contents object for this page is (pageInfos[i].objIndex + 1) */
		pdf->pageInfos[i].contents = objectCount++;	/* advance count for Contents for this page */
		cpdf_finalizePage(pdf, i);				/* close page if not already */
	    }
	}
	for(i=0; i< pdf->numFonts; i++) {
	    pdf->fontInfos[i].objIndex = objectCount++;		/* for Page object Resources entry */
	    if(pdf->fontInfos[i].descLevel > 0 && pdf->fontInfos[i].fontDesc->doneAlready == 0) {
		pdf->fontInfos[i].fontDesc->objIndex = objectCount++;	/* for non-base-14 fonts */
		pdf->fontInfos[i].fontDesc->doneAlready = 1;
		/* descendant font for CJK CIDFontType2 font */
		if(pdf->fontInfos[i].descLevel == 3 && pdf->fontInfos[i].descendFont->doneAlready == 0) {
		    pdf->fontInfos[i].descendFont->objIndex = objectCount++;	/* for CJK font descendant font */
		    pdf->fontInfos[i].descendFont->doneAlready = 1;
		}
		/* object for external font data */
		if(pdf->fontInfos[i].descLevel == 2 && pdf->fontInfos[i].fontDesc->extFont
					&& pdf->fontInfos[i].fontDesc->extFont->doneAlready == 0) {
		    pdf->fontInfos[i].fontDesc->extFont->objIndex = objectCount++;
		    pdf->fontInfos[i].fontDesc->extFont->doneAlready = 1;
		}
	    }
	}
	for(i=0; i< pdf->numImages; i++)
	    pdf->imageInfos[i].objIndex = objectCount++;
	for(i=0; i< pdf->numAnnots; i++)
	    pdf->annotInfos[i].objIndex = objectCount++;

	/* now do outline entries */
	if(pdf->numOutlineEntries)
	    _cpdf_serializeOutlineEntries(&objectCount, pdf->firstOLentry, pdf->lastOLentry);

	pdf->objIndex[CPDF_Info] = objectCount++;

	/* =========== Now, object indexed, so write objects to output file ======================= */
	_cpdf_WriteCatalogObject(pdf, pdf->objIndex[CPDF_Catalog]);
	_cpdf_WriteOutlinesObject(pdf, pdf->objIndex[CPDF_Outlines]);
	_cpdf_WritePagesObject(pdf, pdf->objIndex[CPDF_Pages]);
	_cpdf_WriteProcSetArray(pdf, pdf->objIndex[CPDF_ProcSet]);

	/* Write out all pages */
	/* Now, multiple pages -- pageInfos[0] is unused. */
	for(i=1; i<= pdf->numPages; i++) {
	    /* only those pages with pagenum != -1 are valid */
	    if(pdf->pageInfos[i].pagenum != -1) {
		_cpdf_WritePageObject(pdf, &pdf->pageInfos[i]);
		if(pdf->useContentMemStream)
		    _cpdf_WriteContentsFromMemory(pdf, &pdf->pageInfos[i]);
		else
		    _cpdf_WriteContentsFromFile(pdf, &pdf->pageInfos[i]);
	    }
	}

	for(i=0; i< pdf->numFonts; i++) {
#ifdef CPDF_DEBUG
	    fprintf(stderr, "Font num: %d\n", i);
#endif
	    _cpdf_WriteFont(pdf, &pdf->fontInfos[i]);
	    if(pdf->fontInfos[i].descLevel > 0 && pdf->fontInfos[i].fontDesc->doneAlready < 2) {
		_cpdf_WriteFontDescriptor(pdf, &pdf->fontInfos[i]);	/* for non-base-14 fonts */
		pdf->fontInfos[i].fontDesc->doneAlready = 2;
		/* descendant font for CJK CIDFontType2 font */
		if(pdf->fontInfos[i].descLevel == 3 && pdf->fontInfos[i].descendFont->doneAlready < 2) {
		    _cpdf_WriteDescendantFont(pdf, &pdf->fontInfos[i]); /* for CJK font descendant font */
		    pdf->fontInfos[i].descendFont->doneAlready = 2;
		}
		/* object for external font data */
		if(pdf->fontInfos[i].descLevel == 2 && pdf->fontInfos[i].fontDesc->extFont
					&& pdf->fontInfos[i].fontDesc->extFont->doneAlready < 2) {
		    _cpdf_WriteFontData(pdf, pdf->fontInfos[i].fontDesc->extFont);
		    pdf->fontInfos[i].fontDesc->extFont->doneAlready = 2;
		}
	    }
	}
	for(i=0; i< pdf->numImages; i++)
	    _cpdf_WriteImage(pdf, &pdf->imageInfos[i]);

	for(i=0; i< pdf->numAnnots; i++)
	    _cpdf_WriteAnnotation(pdf, &pdf->annotInfos[i]);

	if(pdf->numOutlineEntries)
	    _cpdf_WriteOutlineEntries(pdf, pdf->firstOLentry, pdf->lastOLentry);

	_cpdf_WriteProducerDate(pdf, pdf->objIndex[CPDF_Info]);
	_cpdf_WriteXrefTrailer(pdf, objectCount);
    } /* end if(linearizeON) {} else {} */
    /* ==== PDF FILE/STREM WRITING DONE --- All objects written to PDF memory stream or file ==== */
    if(!(pdf->usePDFMemStream))
	_cpdf_file_close(pdf);		/* close the PDF file */
				        /* but keep PDF memory stream open if that is used */
#ifdef MacOS8
    if( pdf->filename_set )
        SetFileInfo(pdf->filenamepath, filetype, creator);	/* set PDF file type */
#endif
}


/* set page size as a string in %%BoundingBox: format in PS */

void cpdf_setPageSize(CPDFdoc *pdf, const char *mboxstr, const char *cboxstr)
{
CPDFpageInfo *pI;
	strncpy(pdf->mediaBox, mboxstr, 62);
	strncpy(pdf->cropBox, cboxstr, 62);
	/* Iffy, but allow changing page size after cpdf_pageInit() */
        pI = &pdf->pageInfos[pdf->currentPage];
	if(pI->mediaBox)
	    free(pI->mediaBox);
	pI->mediaBox = (char *)malloc((size_t)(strlen(pdf->mediaBox) + 1));
	_cpdf_malloc_check((void *)pI->mediaBox);
	strcpy(pI->mediaBox, pdf->mediaBox);
	if(pI->cropBox)
	    free(pI->cropBox);
	pI->cropBox = (char *)malloc((size_t)(strlen(pdf->cropBox) + 1));
	_cpdf_malloc_check((void *)pI->cropBox);
	strcpy(pI->cropBox, pdf->cropBox);
	return;
}


void cpdf_setBoundingBox(CPDFdoc *pdf, int LLx, int LLy, int URx, int URy)
{
	cpdf_setMediaBox(pdf, LLx, LLy, URx, URy);
	cpdf_setCropBox(pdf, LLx, LLy, URx, URy);
}


void cpdf_setMediaBox(CPDFdoc *pdf, int LLx, int LLy, int URx, int URy)
{
char tstr[64];
	sprintf(tstr, "%d %d %d %d", LLx, LLy, URx, URy);
	strncpy(pdf->mediaBox, tstr, 62);
}


void cpdf_setCropBox(CPDFdoc *pdf, int LLx, int LLy, int URx, int URy)
{
char tstr[64];
	sprintf(tstr, "%d %d %d %d", LLx, LLy, URx, URy);
	strncpy(pdf->cropBox, tstr, 62);
}


CPDFmemStream *cpdf_setCurrentMemoryStream(CPDFdoc *pdf, CPDFmemStream *memStream)
{
CPDFmemStream *oldMS;
	oldMS = pdf->currentMemStream;
	pdf->currentMemStream = memStream;
	return oldMS;
}

/* Returns the size of file in bytes */

long getFileSize(const char *file)
{
long filesize = 0;
#if defined(_WIN32) || defined(WIN32)
    #ifndef __BORLANDC__
        struct _stat filestat;		/* For Windows VC++ */
    #else
        struct stat filestat;
    #endif
#else
struct stat filestat;
#endif

#if defined(_WIN32) || defined(WIN32)
	if(_stat(file, &filestat))	/* Windows */
#else
	if(stat(file, &filestat))
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
	return(filesize);
}

/* This functions sets preferences that control viewer display options when the document
   is opened, such as whether outline (book marks) display or thumnails should be open
   upon doc opening.
*/

/* API changed more than the addition of PDF context argument.  Now, preferences are
   passed as a struct for possible future expansion.  Backward compatibility retained as above.
*/
void cpdf_setViewerPreferences(CPDFdoc *pdf, CPDFviewerPrefs *vP)
{
    if(vP == NULL) return;
    pdf->viewerPrefs.pageMode = vP->pageMode;	/* This belongs directly to Catalog obj, but here for convenience */
    pdf->viewerPrefs.hideToolbar = vP->hideToolbar;	/* when YES, tool bar in viewer will be hidden */
    pdf->viewerPrefs.hideMenubar = vP->hideMenubar;	/* when YES, menu bar in viewer will be hidden */
    pdf->viewerPrefs.hideWindowUI= vP->hideWindowUI;	/* when YES, GUI elements in doc window will be hidden */
    pdf->viewerPrefs.fitWindow= vP->fitWindow;		/* when YES, viewer resizes the doc window to the page size */
    pdf->viewerPrefs.centerWindow= vP->centerWindow;	/* when YES, viewer moves doc window to screen's center */
    pdf->viewerPrefs.pageLayout= vP->pageLayout;	/* 1-page display or two-page columnar displays */
    pdf->viewerPrefs.nonFSPageMode= vP->nonFSPageMode;  /* pageMode coming out of the full-screen mode */
}


/* ---------------------------------------------------------------------------------------- */
/*   Private functions -- Do not call functions below directly.
     These private functions don't need to be backward compatible.

   This function is called from cpdf_finalizePage(page), and compresses the
   (uncompressed) page content memory stream, alloc memory for compressed buffer,
   deallocates memory stream, and set various flags that indicates this fact.
   This fuction is impletemted to reduce memory usage for an application that
   generates large number of pages, so that uncompressed memory streams are
   not kept hanging around until the end.
*/
int  _cpdf_closeContentMemStreamForPage(CPDFdoc *pdf, int page)
{
Byte *inputBuf, *outputBuf;
CPDFpageInfo *pInf;
uLong comprLen = 0;
int filesize;
float headroom;
int iheadroom;
int bufSize;
int err;

	pInf = &pdf->pageInfos[page];
	cpdf_getMemoryBuffer(pInf->pageMemStream, (char **)&inputBuf, &filesize, &bufSize);

	if(pdf->compressionON) {
	    /* zlib compress needs a little head room for output buffer size */
	    headroom = (float)filesize * 0.001 + 16.0;	/* 0.1% plus 12 bytes */
	    iheadroom = (int)headroom;
	    /* fprintf(stderr, "filesize=%d, bufSize=%d, headroom=%d\n", filesize, bufSize, iheadroom); */
	    comprLen = filesize+iheadroom;		/* available buf size - must set this for input */
	    outputBuf = (Byte *)malloc((size_t)comprLen);	/* working buffer */
	    _cpdf_malloc_check((void *)outputBuf);
	    err = compress(outputBuf, &comprLen, inputBuf, filesize);	/* #### COMPRESSION #### */
	    if(err != Z_OK)
		cpdf_Error(pdf, -1, "ZLIB compress", "code=%d", err);
	    /* Compression done, malloc new memory just the size of compressed stream and copy content */
	    pInf->compressedStream = (char *)malloc((size_t)(comprLen+1));
	    _cpdf_malloc_check((void *)pInf->compressedStream);
	    /* copy compressed data */
	    memcpy((void *)pInf->compressedStream, (void *)outputBuf, (size_t)comprLen);
	    free(outputBuf);			/* working buffer not needed now */
	    pInf->length_compressed = comprLen;	/* indicate compressed data length */
	    pInf->status |= 2;			/* set bit-1 to indicate compressed status */
	    /* fprintf(stderr, "input: %ld,  output: %ld\n", filesize, comprLen); */
	    if(pInf->pageMemStream) {
		cpdf_closeMemoryStream(pInf->pageMemStream);	/* free uncompressed memory stream */
		pInf->pageMemStream = NULL;
	    }
	}
	return(0);
}


int  _cpdf_freeAllPageInfos(CPDFdoc *pdf)
{
CPDFpageInfo *tpage;
CPDFintList *iL, *iLprev;
int i;
    for(i=1; i<= pdf->numPages; i++) {
	tpage = &pdf->pageInfos[i];
	if(tpage->pagenum == -1)
	    continue;		    /* skip unused pages */
	if(tpage->pageMemStream) {  /* _cpdf_WriteContentsFromMemory() should have freed it but... */
	    cpdf_closeMemoryStream(tpage->pageMemStream);
	    tpage->pageMemStream = NULL;
	}
	if(tpage->compressedStream) {  /* _cpdf_WriteContentsFromMemory() should have freed it but... */
	    free(tpage->compressedStream);
	    tpage->compressedStream = NULL;
	}
	if(tpage->defDomain) {
	    cpdf_freePlotDomain(tpage->defDomain);
	    tpage->defDomain = NULL;
	}
	if(tpage->contentfile) {
	    free(tpage->contentfile);
	    tpage->contentfile = NULL;
	}
 	iL = tpage->fontIdx;
	while(iL != NULL) {
	    iLprev = iL;
	    iL = iL->next;
	    free(iLprev);
	}
 	iL = tpage->imageIdx;
	while(iL != NULL) {
	    iLprev = iL;
	    iL = iL->next;
	    free(iLprev);
	}
 	iL = tpage->annotIdx;
	while(iL != NULL) {
	    iLprev = iL;
	    iL = iL->next;
	    free(iLprev);
	}
	if(tpage->mediaBox) {
	    free(tpage->mediaBox);
	    tpage->mediaBox = NULL;
	}
	if(tpage->cropBox) {
	    free(tpage->cropBox);
	    tpage->cropBox = NULL;
	}
	if(tpage->transition) {
	    free(tpage->transition);
	    tpage->transition = NULL;
	}
    } /* end for(i...) */
    return(0);
}


int _cpdf_file_open(CPDFdoc *pdf)
{
char pps[]=".pdf";
#if defined(_WIN32) || defined(WIN32)
extern int setmode(int filenum, int mode);
#endif
#if defined(UNIX) && !defined(SunOS5x)
extern pid_t getpid(void);
int mypid = getpid();
#else
int mypid = 123;
#endif

	if(!(pdf->filename_set)) {
	    strcpy(pdf->filenamepath, TEMP_DIRECTORY);
	    strcat(pdf->filenamepath, "_cpdf");
	    str_append_int(pdf->filenamepath, mypid);
	    strcat(pdf->filenamepath, "-");
	    str_append_int(pdf->filenamepath, pdf->docID);
	    strcat(pdf->filenamepath, pps);
	    pdf->filename_set = 1;
	}

	if(strcmp(pdf->filenamepath, "-") == 0)
	    pdf->useStandardOutput = 1;
	if(pdf->useStandardOutput) {
#if defined(_WIN32) || defined(WIN32)
	    pdf->fpcg = stdout;
	    setmode(fileno(pdf->fpcg), O_BINARY);
#else
	    pdf->fpcg = fdopen(1, BINARY_WRITE);	/* fd of 1 == stdout */
#endif
	}
	else
	    pdf->fpcg = fopen(pdf->filenamepath, BINARY_WRITE);	/* open for write */
	    /*  Open PDF file in binary mode to turn off LF -> CRLF translation.
		This is necessary because we are doing byte offset counting by strlen()
		function which assumes '\n' to be the end-of-line char.  Otherwise, the
		offsets in PDF xref table will be incorrect, and Acrobat Reader will give
		a warning about the damaged file (though it can read it fine).
	    */

	if(pdf->fpcg == NULL) {
		cpdf_Error(pdf, 2, "ClibPDF", "Cannot open %s for PDF output", pdf->filenamepath);
		return(1);
	}

	/* pdf->fncounter++; */		/* multiple files within a thread */
	return(0);			/* successful file/device open */
}

/* Close the PDF document file.  This gets called if PDF output is sent directly to a file
   without the use of memory stream for PDF.  */
void _cpdf_file_close(CPDFdoc *pdf)
{
    if(pdf->fpcg) {
	fclose(pdf->fpcg);
	pdf->fpcg = NULL;
    }
}



/* ---------------------------------------------------------------------------------------- */
/*   Private functions that write PDF objects */
/* ---------------------------------------------------------------------------------------- */
static char *pmStr[] = { "UseNone", "UseOutlines", "UseThumbs", "FullScreen", ""};
static char *plStr[] = { "SinglePage", "OneColumn", "TwoColumnLeft", "TwoColumnRight", ""};

long _cpdf_WriteCatalogObject(CPDFdoc *pdf, int objNumber)
{
	sprintf(pdf->spbuf, "%d 0 obj\n", objNumber); 				_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "<<\n"); 						_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Type /Catalog\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Pages %d 0 R\n", pdf->objIndex[CPDF_Pages]);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Outlines %d 0 R\n", pdf->objIndex[CPDF_Outlines]);
	_cpdf_pdfWrite(pdf, pdf->spbuf);
	if(pdf->viewerPrefs.pageMode) {
	    sprintf(pdf->spbuf, "/PageMode /%s\n", pmStr[pdf->viewerPrefs.pageMode & 3]);
	    _cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	/* There was a documentation bug in PDF Spec v1.2 about where /PageLayout key should
	   be present.  It should be in the Catalog object and not in the /ViewerPreferences dict.
	   It is fixed here, but Viewer Pref C struct is kept as before.  PDF Spec v1.3 had
	   this note.  Thanks to James W. Walker for pointing this out.
	*/
	if(pdf->viewerPrefs.pageLayout) {
	    sprintf(pdf->spbuf, "/PageLayout /%s\n", plStr[pdf->viewerPrefs.pageLayout & 3]);
	    _cpdf_pdfWrite(pdf, pdf->spbuf);
	}

	if(pdf->viewerPrefs.pageMode==3 || pdf->viewerPrefs.hideToolbar || pdf->viewerPrefs.hideMenubar ||
	   pdf->viewerPrefs.hideWindowUI || pdf->viewerPrefs.fitWindow || pdf->viewerPrefs.centerWindow) {
	    sprintf(pdf->spbuf, "/ViewerPreferences <<\n"); 			_cpdf_pdfWrite(pdf, pdf->spbuf);
	    if(pdf->viewerPrefs.hideToolbar) {
		sprintf(pdf->spbuf, "/HideToolbar true\n");			_cpdf_pdfWrite(pdf, pdf->spbuf);
	    }
	    if(pdf->viewerPrefs.hideMenubar) {
		sprintf(pdf->spbuf, "/HideMenubar true\n");			_cpdf_pdfWrite(pdf, pdf->spbuf);
	    }
	    if(pdf->viewerPrefs.hideWindowUI) {
		sprintf(pdf->spbuf, "/HideWindowUI true\n");			_cpdf_pdfWrite(pdf, pdf->spbuf);
	    }
	    if(pdf->viewerPrefs.fitWindow) {
		sprintf(pdf->spbuf, "/FitWindow true\n");			_cpdf_pdfWrite(pdf, pdf->spbuf);
	    }
	    if(pdf->viewerPrefs.centerWindow) {
		sprintf(pdf->spbuf, "/CenterWindow true\n");			_cpdf_pdfWrite(pdf, pdf->spbuf);
	    }
	    if(pdf->viewerPrefs.pageMode==3) {
		if(pdf->viewerPrefs.nonFSPageMode < 0 || pdf->viewerPrefs.nonFSPageMode > 2)
		    pdf->viewerPrefs.nonFSPageMode = 0;	/* full-screen mode is now allowed for this */
		sprintf(pdf->spbuf, "/NonFullScreenPageMode /%s\n", pmStr[pdf->viewerPrefs.nonFSPageMode]);
		_cpdf_pdfWrite(pdf, pdf->spbuf);
	    }
	    sprintf(pdf->spbuf, ">>\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	sprintf(pdf->spbuf, ">>\n");						_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "endobj\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	pdf->objByteOffset[objNumber+1] = pdf->currentByteCount;
	return(pdf->currentByteCount);
}

long _cpdf_WriteOutlinesObject(CPDFdoc *pdf, int objNumber)
{
	sprintf(pdf->spbuf, "%d 0 obj\n", objNumber);				_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "<<\n");						_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Type /Outlines\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Count %d\n", pdf->numOutlineEntries);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	if(pdf->numOutlineEntries) {
	    sprintf(pdf->spbuf, "/First %d 0 R\n", pdf->firstOLentry->objIndex);
	    _cpdf_pdfWrite(pdf, pdf->spbuf);
	    sprintf(pdf->spbuf, "/Last %d 0 R\n", pdf->lastOLentry->objIndex);
	    _cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	sprintf(pdf->spbuf, ">>\n");						_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "endobj\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	pdf->objByteOffset[objNumber+1] = pdf->currentByteCount;
	return(pdf->currentByteCount);
}

long _cpdf_WritePagesObject(CPDFdoc *pdf, int objNumber)
{
int i;
	sprintf(pdf->spbuf, "%d 0 obj\n", objNumber);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "<<\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Type /Pages\n");			_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Count %d\n", pdf->numKids);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Kids [ ");  			_cpdf_pdfWrite(pdf, pdf->spbuf);
	for(i=0; i< pdf->numKids; i++) {
	    sprintf(pdf->spbuf, "%d 0 R ", pdf->kidsIndex[i]);
	    _cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	sprintf(pdf->spbuf, "]\n>>\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "endobj\n");			_cpdf_pdfWrite(pdf, pdf->spbuf);
	pdf->objByteOffset[objNumber+1] = pdf->currentByteCount;
	return(pdf->currentByteCount);
}

long _cpdf_WritePageObject(CPDFdoc *pdf, CPDFpageInfo *pInf)
{
CPDFfontInfo *tfont;
CPDFimageInfo *timg;
CPDFannotInfo *tann;
CPDFintList *iL;

	sprintf(pdf->spbuf, "%d 0 obj\n", pInf->objIndex);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "<<\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Type /Page\n");			_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Parent %d 0 R\n", pInf->parent);  _cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Resources <<\n");			_cpdf_pdfWrite(pdf, pdf->spbuf);

	if(pInf->npFont) {
	    sprintf(pdf->spbuf, "/Font <<\n");	_cpdf_pdfWrite(pdf, pdf->spbuf);
	    iL = pInf->fontIdx;
	    while(iL != NULL) {
		if(iL->value >= 0) {
		    tfont = &pdf->fontInfos[iL->value];
		    sprintf(pdf->spbuf, "/%s %d 0 R\n", tfont->name, tfont->objIndex);
		    _cpdf_pdfWrite(pdf, pdf->spbuf);
		}
		iL = iL->next;
	    }
	    sprintf(pdf->spbuf, ">>\n"); _cpdf_pdfWrite(pdf, pdf->spbuf);
	}

	if(pInf->npImage) {
	    sprintf(pdf->spbuf, "/XObject <<\n");		_cpdf_pdfWrite(pdf, pdf->spbuf);
	    iL = pInf->imageIdx;
	    while(iL != NULL) {
	        if(iL->value >= 0) {
		    timg = &pdf->imageInfos[iL->value];
		    sprintf(pdf->spbuf, "/%s %d 0 R\n", timg->name, timg->objIndex);
		    _cpdf_pdfWrite(pdf, pdf->spbuf);
		}
		iL = iL->next;
	    }
	    sprintf(pdf->spbuf, ">>\n"); _cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	sprintf(pdf->spbuf, "/ProcSet %d 0 R >>\n", pdf->objIndex[CPDF_ProcSet]); _cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/MediaBox [%s]\n", pInf->mediaBox);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/CropBox [%s]\n", pInf->cropBox);			_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Rotate %d\n", pInf->orientation);			_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Contents %d 0 R\n", pInf->contents);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	if(pInf->npAnnot) {
	    sprintf(pdf->spbuf, "/Annots [ "); 					_cpdf_pdfWrite(pdf, pdf->spbuf);
	    iL = pInf->annotIdx;
	    while(iL != NULL) {
	        if(iL->value >= 0) {
		    tann = &pdf->annotInfos[iL->value];
		    sprintf(pdf->spbuf, "%d 0 R ", tann->objIndex);
		    _cpdf_pdfWrite(pdf, pdf->spbuf);
		}
		iL = iL->next;
	    }
	    sprintf(pdf->spbuf, "]\n"); 					_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	if(pInf->duration > 0.0) {
	    sprintf(pdf->spbuf, "/Dur %.3f\n", pInf->duration);			_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	if(pInf->transition) {
	    sprintf(pdf->spbuf, "/Trans << %s >>\n", pInf->transition);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	sprintf(pdf->spbuf, ">>\n");						_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "endobj\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	pdf->objByteOffset[pInf->objIndex+1] = pdf->currentByteCount;
	return(pdf->currentByteCount);
}

long _cpdf_WriteContentsFromFile(CPDFdoc *pdf, CPDFpageInfo *pInf)
{
long	filesize;
#ifdef EXT_ZLIBCOMP
int c;
char syscomdbuf[1024];
#else
Byte *inputBuf, *outputBuf;
uLong comprLen = 0;
int err;
#endif

	strcpy(pdf->contentfile, pInf->contentfile);		/* filaname for this page */

#ifdef EXT_ZLIBCOMP
	if(pdf->compressionON) {
	    sprintf(syscomdbuf, "%s %s", pdf->compress_command, pdf->contentfile);
	    system(syscomdbuf);
	    strcat(pdf->contentfile, ".zlib");
	}
#endif

	filesize = getFileSize(pdf->contentfile);
	sprintf(pdf->spbuf, "%d 0 obj\n", pInf->contents);		_cpdf_pdfWrite(pdf, pdf->spbuf);

#ifdef EXT_ZLIBCOMP
	if(pdf->compressionON) {
	    /* Here, we use external zlibcomp command for compression, so use filesize */
	    sprintf(pdf->spbuf, "<<\n/Length %ld\n/Filter [%s]\n>>\nstream\n", filesize, pdf->streamFilterList);
	    _cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	else {
	    sprintf(pdf->spbuf, "<</Length %ld>>\nstream\n", filesize);  _cpdf_pdfWrite(pdf, pdf->spbuf);
	}

	/* copy file from fpcontent to fpcg */
	if((pdf->fpcontent = fopen(pdf->contentfile, BINARY_READ)) != NULL) {
	    while((c = fgetc(pdf->fpcontent)) != EOF)
		    fputc(c, pdf->fpcg);
	    fclose(pdf->fpcontent);
	}
	pdf->currentByteCount += filesize;				/* add size of Contents file */
#else
/* else of #ifdef EXT_ZLIBCOMP -- So, here, we do zlib compression internally. */
	inputBuf = (Byte *)malloc((size_t)(filesize+16));
	_cpdf_malloc_check((void *)inputBuf);
	if((pdf->fpcontent = fopen(pdf->contentfile, BINARY_READ)) != NULL) {
	    fread(inputBuf, 1, (size_t)filesize, pdf->fpcontent);
	    fclose(pdf->fpcontent);
	}
	if(pdf->compressionON) {
	    comprLen = filesize+1024;		/* must set this for input */
	    outputBuf = (Byte *)malloc((size_t)comprLen);
	    _cpdf_malloc_check((void *)outputBuf);
	    err = compress(outputBuf, &comprLen, inputBuf, filesize);
	    if(err != Z_OK)
		cpdf_Error(pdf, -1, "ZLIB compress", "code=%d", err);
	    /* fprintf(stderr, "input: %ld,  output: %ld\n", filesize, comprLen); */
	    sprintf(pdf->spbuf, "<<\n/Length %ld\n/Filter [%s]\n>>\nstream\n", comprLen, pdf->streamFilterList);
	    _cpdf_pdfWrite(pdf, pdf->spbuf);
	    if(pdf->usePDFMemStream)
	        cpdf_writeMemoryStream( pdf->pdfMemStream, (char *)outputBuf, comprLen);
	    else
	        fwrite(outputBuf, 1, (size_t)comprLen, pdf->fpcg);
	    pdf->currentByteCount += comprLen;				/* add size of Contents file */
	    free(outputBuf);
	}
	else {
	    /* No compression */
	    sprintf(pdf->spbuf, "<</Length %ld>>\nstream\n", filesize); _cpdf_pdfWrite(pdf, pdf->spbuf);
	    if(pdf->usePDFMemStream)
	        cpdf_writeMemoryStream( pdf->pdfMemStream, (char *)inputBuf, filesize);
	    else
	        fwrite(inputBuf, 1, (size_t)filesize, pdf->fpcg);
	    pdf->currentByteCount += filesize;				/* add size of Contents file */
	}
	free(inputBuf);
#endif
/* end of #ifdef EXT_ZLIBCOMP */

#ifndef	CPDF_DEBUG
	remove(pdf->contentfile);				/* delete temporary content file */
#endif

	sprintf(pdf->spbuf, "\nendstream\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "endobj\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	pdf->objByteOffset[pInf->contents+1] = pdf->currentByteCount;
	return(pdf->currentByteCount);
}


long _cpdf_WriteContentsFromMemory(CPDFdoc *pdf, CPDFpageInfo *pInf)
{
Byte *inputBuf;
int filesize;
int bufSize;

	sprintf(pdf->spbuf, "%d 0 obj\n", pInf->contents);		_cpdf_pdfWrite(pdf, pdf->spbuf);

	if(pInf->status & 2) {
	    /* Compressed page content */
	    /* fprintf(stderr, "input: %ld,  output: %ld\n", filesize, pInf->length_compressed); */
	    sprintf(pdf->spbuf, "<<\n/Length %ld\n/Filter [%s]\n>>\nstream\n",
					pInf->length_compressed, pdf->streamFilterList);
	    _cpdf_pdfWrite(pdf, pdf->spbuf);
	    if(pdf->usePDFMemStream)
	        cpdf_writeMemoryStream( pdf->pdfMemStream, pInf->compressedStream, pInf->length_compressed);
	    else
	        fwrite((void *)pInf->compressedStream, 1, (size_t)(pInf->length_compressed), pdf->fpcg);
	    pdf->currentByteCount += pInf->length_compressed;	/* add size of Contents file */
	    free(pInf->compressedStream);
	    pInf->compressedStream = NULL;
	}
	else {
	    /* No compression */
	    cpdf_getMemoryBuffer(pInf->pageMemStream, (char **)&inputBuf, &filesize, &bufSize);
	    sprintf(pdf->spbuf, "<</Length %d>>\nstream\n", filesize); _cpdf_pdfWrite(pdf, pdf->spbuf);
	    if(pdf->usePDFMemStream)
	        cpdf_writeMemoryStream( pdf->pdfMemStream, (char *)inputBuf, filesize);
	    else
	        fwrite(inputBuf, 1, (size_t)filesize, pdf->fpcg);
	    pdf->currentByteCount += filesize;				/* add size of Contents file */
	    cpdf_closeMemoryStream(pInf->pageMemStream);
	    pInf->pageMemStream = NULL;
	}

	sprintf(pdf->spbuf, "\nendstream\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "endobj\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	pdf->objByteOffset[pInf->contents+1] = pdf->currentByteCount;
	return(pdf->currentByteCount);
}



long _cpdf_WriteProcSetArray(CPDFdoc *pdf, int objNumber)
{
	sprintf(pdf->spbuf, "%d 0 obj\n", objNumber);			_cpdf_pdfWrite(pdf, pdf->spbuf);
	if(pdf->numImages || pdf->inlineImages) {
	    sprintf(pdf->spbuf, "[/PDF /Text ");			_cpdf_pdfWrite(pdf, pdf->spbuf);
	    if(pdf->imageFlagBCI & 1) {
		sprintf(pdf->spbuf, "/ImageB ");		/* gray scale images */
		_cpdf_pdfWrite(pdf, pdf->spbuf);
	    }
	    if(pdf->imageFlagBCI & 2) {
		sprintf(pdf->spbuf, "/ImageC ");		/* color images */
		_cpdf_pdfWrite(pdf, pdf->spbuf);
	    }
	    if(pdf->imageFlagBCI & 4) {
		sprintf(pdf->spbuf, "/ImageI ");		/* indexed color images */
		_cpdf_pdfWrite(pdf, pdf->spbuf);
	    }
	    sprintf(pdf->spbuf, "]\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	else {
		sprintf(pdf->spbuf, "[/PDF /Text]\n");
		_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	sprintf(pdf->spbuf, "endobj\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	pdf->objByteOffset[objNumber+1] = pdf->currentByteCount;
	return(pdf->currentByteCount);
}

/* long _cpdf_WriteFont(int objNumber, char *fontName, char *baseFont, char *encoding) */
long _cpdf_WriteFont(CPDFdoc *pdf, CPDFfontInfo *fI)
{
int i, wcnt = 0;
	sprintf(pdf->spbuf, "%d 0 obj\n", fI->objIndex);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "<<\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Type /Font\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	if(fI->descLevel < 3)
	    sprintf(pdf->spbuf, "/Subtype /Type1\n");
	else
	    sprintf(pdf->spbuf, "/Subtype /Type0\n");		/* for CJK */
	_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Name /%s\n", fI->name);			_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/BaseFont /%s\n", fI->baseFont);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	/* Don't write encoding for Symbol, ZapfDingBats even if specified */
	if(fI->encoding != NULL && fI->afmIndex != 12 && fI->afmIndex != 13) {
	    sprintf(pdf->spbuf, "/Encoding /%s\n", fI->encoding);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
/* New addition starting with v1.10: 1999-03-07 */
	if(fI->descLevel > 0 && fI->descLevel < 3) {
	    /* encodeIndex: 0=WinAnsiEncoding, 1=MacRomanEncoding, 2=MacExpertEncoding, 3=StandardEncoding */
	    int enc = pdf->fontInfos[pdf->currentFont].encodeIndex;
	    /* For now, external fonts use width[] straight from PFM, and cannot use any encoding
	       other than the font's built-in encoding.
	    */
	    if(fI->descLevel == 2 || fI->afmIndex == 12 || fI->afmIndex == 13)
	        enc = 0;	/* Symbol, ZapfDingbats, or external fonts should not use a LUT. */

	    sprintf(pdf->spbuf, "/FirstChar %d\n", fI->firstchar);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	    sprintf(pdf->spbuf, "/LastChar %d\n", fI->lastchar);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	    sprintf(pdf->spbuf, "/Widths [\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	    for(i=fI->firstchar; i <= fI->lastchar; i++) {
		int ch = i;
	        if(enc)
		    ch = cpdf_encode2WinAnsi[enc-1][ch];	/* Look up, if encoding is not WinAnsiEncoding */
		wcnt++;
		if(wcnt%16 == 0) _cpdf_pdfWrite(pdf, "\n");
		sprintf(pdf->spbuf, "%d ", fI->width[ch]);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	    }
	    sprintf(pdf->spbuf, "]\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	    sprintf(pdf->spbuf, "/FontDescriptor %d 0 R\n", fI->fontDesc->objIndex);
	    _cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	/* For CJK CIDFontType2 */
	if(fI->descLevel == 3) {
	    sprintf(pdf->spbuf, "/DescendantFonts [ %d 0 R ]\n", fI->descendFont->objIndex);
	    _cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	sprintf(pdf->spbuf, ">>\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "endobj\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	pdf->objByteOffset[fI->objIndex+1] = pdf->currentByteCount;
	return(pdf->currentByteCount);
}

/* Added starting with v1.10: 1999-03-07
   This function is called only for non-base-14 fonts.
*/
long _cpdf_WriteFontDescriptor(CPDFdoc *pdf, CPDFfontInfo *fInf)
{
CPDFfontDescriptor *fD;
	fD = fInf->fontDesc;
	sprintf(pdf->spbuf, "%d 0 obj\n", fD->objIndex);			_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "<<\n");						_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Type /FontDescriptor\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Ascent %d\n", fD->ascent);			_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/CapHeight %d\n", fD->capHeight);			_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Descent %d\n", fD->descent);			_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Flags %ld\n", fD->flags);				_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/FontBBox [%d %d %d %d]\n",
	  fD->fontBBox[0], fD->fontBBox[1], fD->fontBBox[2], fD->fontBBox[3]);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/FontName /%s\n", fD->fontName);			_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/ItalicAngle %g\n", fD->italicAngle);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/StemV %d\n", fD->stemV);				_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/XHeight %d\n", fD->xHeight);			_cpdf_pdfWrite(pdf, pdf->spbuf);
	/* For CJK fonts */
	if(fInf->descLevel == 3) {
	    sprintf(pdf->spbuf, "/StemH %d\n", fD->stemH); 			_cpdf_pdfWrite(pdf, pdf->spbuf);
	    sprintf(pdf->spbuf, "/MissingWidth %d\n", fD->missingWidth);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	    sprintf(pdf->spbuf, "/Leading %d\n", fD->leading); 			_cpdf_pdfWrite(pdf, pdf->spbuf);
	    sprintf(pdf->spbuf, "/MaxWidth %d\n", fD->maxWidth); 		_cpdf_pdfWrite(pdf, pdf->spbuf);
	    sprintf(pdf->spbuf, "/AvgWidth %d\n", fD->avgWidth); 		_cpdf_pdfWrite(pdf, pdf->spbuf);
	    sprintf(pdf->spbuf, "/Style << %s >>\n", fD->style); 		_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	if(fD->extFont) {
	    sprintf(pdf->spbuf, "/FontFile %d 0 R\n", fD->extFont->objIndex);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	sprintf(pdf->spbuf, ">>\n");						_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "endobj\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	pdf->objByteOffset[fD->objIndex+1] = pdf->currentByteCount;
	return(pdf->currentByteCount);
}

/* For CJK CIDFontType2 Descendant font */
long _cpdf_WriteDescendantFont(CPDFdoc *pdf, CPDFfontInfo *fInf)
{
CPDFdescendantFont *fD;
	fD = fInf->descendFont;
	sprintf(pdf->spbuf, "%d 0 obj\n", fD->objIndex);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "<<\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Type /Font\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Subtype /CIDFontType2\n");		_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/BaseFont /%s\n", fInf->baseFont);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/FontDescriptor %d 0 R\n", fInf->fontDesc->objIndex);
		_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/CIDSystemInfo << %s >>\n", fD->cidSysInfo);  _cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/DW %d\n", fD->DW);			_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/W %s\n", fD->cidFontWidth);  		_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, ">>\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "endobj\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	pdf->objByteOffset[fD->objIndex+1] = pdf->currentByteCount;
	return(pdf->currentByteCount);
}


long _cpdf_WriteFontData(CPDFdoc *pdf, CPDFextFont *eFD)
{
char *memBuffer;
int memLen, bufSize;
	sprintf(pdf->spbuf, "%d 0 obj\n", eFD->objIndex);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "<<\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Length %ld\n", eFD->length);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Length1 %ld\n", eFD->length1);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Length2 %ld\n", eFD->length2);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Length3 %ld\n", eFD->length3);		_cpdf_pdfWrite(pdf, pdf->spbuf);

	sprintf(pdf->spbuf, ">>\nstream\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	cpdf_getMemoryBuffer(eFD->fontMemStream, &memBuffer, &memLen, &bufSize);
	if(pdf->usePDFMemStream)
	    cpdf_writeMemoryStream( pdf->pdfMemStream, (char *)memBuffer, eFD->length);
	else
	    fwrite(memBuffer, 1, (size_t)eFD->length, pdf->fpcg);
	pdf->currentByteCount += eFD->length;			/* add length of this stream */

	sprintf(pdf->spbuf, "\nendstream\nendobj\n");			_cpdf_pdfWrite(pdf, pdf->spbuf);
	pdf->objByteOffset[eFD->objIndex+1] = pdf->currentByteCount;
	return(pdf->currentByteCount);
}


/*                                 0               1            2           3              4 */
static char *colorspaces[] = { "DeviceGray", "DeviceGray", "DeviceGray", "DeviceRGB", "DeviceCMYK" };
static char *invarray[] =   {  "[1 0]",      "[1 0]",      "[1 0]",      "[1 0 1 0 1 0]", "[1 0 1 0 1 0]" };
/* colorspaces[2] (when number of components is 2) is probably GrayAlpha so DeviceGray */

long _cpdf_WriteImage(CPDFdoc *pdf, CPDFimageInfo *imgInf)
{
int ImDataOffset = -1;
int Invert = imgInf->invert;
FILE *fpimg;
char *inputBuf = NULL;
	sprintf(pdf->spbuf, "%d 0 obj\n", imgInf->objIndex);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	_cpdf_pdfWrite(pdf, "<<\n");
	_cpdf_pdfWrite(pdf, "/Type /XObject\n");
	_cpdf_pdfWrite(pdf, "/Subtype /Image\n");
	sprintf(pdf->spbuf, "/Name /%s\n", imgInf->name);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	if(imgInf->type != CPDF_IMG) {
	    sprintf(pdf->spbuf, "/Width %d\n", imgInf->width);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	    sprintf(pdf->spbuf, "/Height %d\n", imgInf->height);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	/* Reading of image data has been deferred until this function.
	   Image header had been read in cpdfImage.c: cpdf_importRawImage() to determine
	   basic image parameters as above.
	   We now read the image data.
	*/
	switch(imgInf->type) {
	    case JPEG_IMG:
		/* Now, copy the JPEG file to output, and adjust currentByteCount.
		   JPEG (baseline) is the only image file format which we can just copy
		   the file in as is.
		*/
		inputBuf = (char *)malloc((size_t)(imgInf->filesize+16));
		_cpdf_malloc_check((void *)inputBuf);
		if((fpimg = fopen(imgInf->filepath, BINARY_READ)) != NULL) {
		    fread(inputBuf, 1, (size_t)(imgInf->filesize), fpimg);
		    fclose(fpimg);
		}
		_cpdf_pdfWrite(pdf, "/Filter /DCTDecode\n");	
		break;

	    case TIFF_IMG:
		cpdf_readTIFFimageData(&inputBuf, imgInf, pdf->compressionON);
		/* imgInf->filesize is set now, and is the compressed image size. */
		if(imgInf->photometric == 0 && 
		    !(imgInf->process == COMPRESSION_CCITTFAX4 && imgInf->NumStrips == 1) ) {
		    /* min is white, then invert BW */
		    if(Invert) Invert = 0;
		    else       Invert = 1;
		}
		if(imgInf->process == COMPRESSION_CCITTFAX4 && imgInf->NumStrips == 1) {
		    _cpdf_pdfWrite(pdf, "/Filter /CCITTFaxDecode\n");
		    sprintf(pdf->spbuf, "/DecodeParms << /K -1 /Columns %d >>\n",
				imgInf->width);
		    _cpdf_pdfWrite(pdf, pdf->spbuf);
		}
		else if(pdf->compressionON) {
		    _cpdf_pdfWrite(pdf, "/Filter /FlateDecode\n");
		}
		break;

	    case GIF_IMG:	/* not implemented */
		break;

	    case CPDF_IMG:	/* PDFIMG - PDF Image XObject dump format */
		ImDataOffset = _cpdf_copyPDFimageHeader(pdf, imgInf);	/* copy the header part to PDF */
		if(ImDataOffset > 0)
		    cpdf_readPDFimageData(&inputBuf, imgInf, ImDataOffset);	/* binary image data */
		break;

	    case RAW_IMG:	/* RAW image data (kept at imgInf->data) */
		inputBuf = (char *) imgInf->data;
		if(pdf->compressionON) {
		    _cpdf_pdfWrite(pdf, "/Filter /FlateDecode\n");
		}
		break;

	    default:
		break;

	} /* end switch(imgInf->type) {} */

	/* filesize had been found earlier in cpdfImage.c: cpdf_importRawImage() */
	if(imgInf->type != CPDF_IMG) {
	    /* CPDF_IMG should have all option settings in it. */
	    if(Invert) {
		sprintf(pdf->spbuf, "/Decode %s\n", invarray[imgInf->ncomponents]);
		_cpdf_pdfWrite(pdf, pdf->spbuf);
	    }
	    sprintf(pdf->spbuf, "/BitsPerComponent %d\n", imgInf->bitspersample);
	    _cpdf_pdfWrite(pdf, pdf->spbuf);
	    if(imgInf->imagemask == 1 && imgInf->ncomponents == 1 && imgInf->bitspersample == 1)
	        sprintf(pdf->spbuf, "/ImageMask true\n");
	    else
	        sprintf(pdf->spbuf, "/ColorSpace /%s\n", colorspaces[imgInf->ncomponents]);
	    _cpdf_pdfWrite(pdf, pdf->spbuf);
	    sprintf(pdf->spbuf, "/Length %ld\n", imgInf->filesize);
	    _cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	_cpdf_pdfWrite(pdf, ">>\n");
	_cpdf_pdfWrite(pdf, "stream\n");

	if(pdf->usePDFMemStream)
	    cpdf_writeMemoryStream( pdf->pdfMemStream, inputBuf, imgInf->filesize);
	else
	    fwrite(inputBuf, 1, (size_t)(imgInf->filesize), pdf->fpcg);
	pdf->currentByteCount += (imgInf->filesize);		/* add size of Contents file */
	if(inputBuf && imgInf->type != RAW_IMG)
	    free(inputBuf);

	_cpdf_pdfWrite(pdf, "\nendstream\n");
	_cpdf_pdfWrite(pdf, "endobj\n");
	pdf->objByteOffset[imgInf->objIndex +1] = pdf->currentByteCount;
	return(pdf->currentByteCount);
}


long _cpdf_WriteAnnotation(CPDFdoc *pdf, CPDFannotInfo *aInf)
{
	sprintf(pdf->spbuf, "%d 0 obj\n", aInf->objIndex);			_cpdf_pdfWrite(pdf, pdf->spbuf);
	_cpdf_pdfWrite(pdf, "<<\n");
	_cpdf_pdfWrite(pdf, "/Type /Annot\n");
	if(aInf->type == ANNOT_TEXT) {
	    _cpdf_pdfWrite(pdf, "/Subtype /Text\n");
	    sprintf(pdf->spbuf, "/Rect [%.0f %.0f %.0f %.0f]\n",
			aInf->xLL, aInf->yLL, aInf->xUR, aInf->yUR);		_cpdf_pdfWrite(pdf, pdf->spbuf);

	    _cpdf_pdfWrite(pdf, "/T (");
	    _cpdf_pdfBinaryWrite(pdf, aInf->annot_title, aInf->title_len);	/* Unicode ready */
	    _cpdf_pdfWrite(pdf, ")\n");

	    _cpdf_pdfWrite(pdf, "/Contents (");
	    _cpdf_pdfBinaryWrite(pdf, aInf->content_link, aInf->content_len);	/* Unicode ready */
	    _cpdf_pdfWrite(pdf, ")\n");
	}
	else if(aInf->type == ANNOT_URL) {
	    _cpdf_pdfWrite(pdf, "/Subtype /Link\n");
	    sprintf(pdf->spbuf, "/Rect [%.0f %.0f %.0f %.0f]\n",
			aInf->xLL, aInf->yLL, aInf->xUR, aInf->yUR);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	    sprintf(pdf->spbuf, "/A << /S /URI\n/URI (%s)\n>>\n", aInf->content_link);
										_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	else if(aInf->type == ANNOT_GOTO) {
	    _cpdf_pdfWrite(pdf, "/Subtype /Link\n");
	    sprintf(pdf->spbuf, "/Rect [%.0f %.0f %.0f %.0f]\n",
			aInf->xLL, aInf->yLL, aInf->xUR, aInf->yUR);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	    sprintf(pdf->spbuf, "/Dest [%d 0 R %s]\n", pdf->pageInfos[aInf->page].objIndex, aInf->content_link);
										_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	else if(aInf->type == ANNOT_ACTION) {
	    _cpdf_pdfWrite(pdf, "/Subtype /Link\n");
	    sprintf(pdf->spbuf, "/Rect [%.0f %.0f %.0f %.0f]\n",
			aInf->xLL, aInf->yLL, aInf->xUR, aInf->yUR);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	    sprintf(pdf->spbuf, "/A <<\n%s\n>>\n", aInf->content_link);
										_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	
	if(aInf->BS != NULL) {
	    sprintf(pdf->spbuf, "/BS <<\n%s\n>>\n", aInf->BS);
										_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	else if(aInf->border_array != NULL) {
	    sprintf(pdf->spbuf, "/Border %s\n", aInf->border_array);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	else {
	    /* Both border_arrary or BS are not specified, put default */
	    sprintf(pdf->spbuf, "/Border [0 0 1]\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	}

	/* Other detail options added with ver 1.10 */
	sprintf(pdf->spbuf, "/C [%.4f %.4f %.4f]\n", aInf->r,  aInf->g,  aInf->b);
	_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/F %d\n", aInf->flags);
	_cpdf_pdfWrite(pdf, pdf->spbuf);

	_cpdf_pdfWrite(pdf, ">>\n");
	_cpdf_pdfWrite(pdf, "endobj\n");
	pdf->objByteOffset[aInf->objIndex+1] = pdf->currentByteCount;
	return(pdf->currentByteCount);
}

long _cpdf_WriteProducerDate(CPDFdoc *pdf, int objNumber)
{
char timebuf[64];
	sprintf(pdf->spbuf, "%d 0 obj\n", objNumber);				_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "<<\n");						_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Creator (%s)\n", pdf->creator_name);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/CreationDate (D:%s)\n", timestring(0, timebuf));
  	_cpdf_pdfWrite(pdf, pdf->spbuf);
/*
   As part of the license, the following conditions must be observed for the
   /Producer field: You may modify the /Producer field provided that the original
   LIBRARY_NAME and CPDF_VERSION are included intact formatted as [%s %s] in the
   modified field.  You may append any additional strings after the above.
   There are no restrictions on other fields.
*/
	sprintf(pdf->spbuf, "/Producer ([%s %s] %s)\n", LIBRARY_NAME, CPDF_VERSION, PLATFORM_NAME);
									_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Author (User: %s)\n", pdf->username);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Title (%s)\n", pdf->file_title);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Subject (%s)\n", pdf->file_subject);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Keywords (%s)\n", pdf->file_keywords);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, ">>\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "endobj\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	pdf->objByteOffset[objNumber+1] = pdf->currentByteCount;
	return(pdf->currentByteCount);
}


long _cpdf_WriteXrefTrailer(CPDFdoc *pdf, int objNumber)
{
int i;
	pdf->startXref = pdf->currentByteCount;
	sprintf(pdf->spbuf, "xref\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "0 %d\n", objNumber);			_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "0000000000 65535 f \n");			_cpdf_pdfWrite(pdf, pdf->spbuf);
	for(i=1; i< objNumber; i++) {
	     sprintf(pdf->spbuf, "%010ld 00000 n \n", pdf->objByteOffset[i]);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	sprintf(pdf->spbuf, "trailer\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "<<\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Size %d\n", objNumber);			_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Root %d 0 R\n", pdf->objIndex[CPDF_Catalog]);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "/Info %d 0 R\n", pdf->objIndex[CPDF_Info]);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, ">>\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "startxref\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "%ld\n", pdf->startXref);			_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "%%%%EOF\n");				_cpdf_pdfWrite(pdf, pdf->spbuf);

	return(pdf->currentByteCount);
}


void _cpdf_initDocumentStruct(CPDFdoc *pdf)
{
    pdf->docID = -1;			/* document ID (thread ID ?) */
    pdf->ps_pdf_mode = 0;		/* PDF=0, EPDF=1, PS=2, EPS=3, FDF=4 (not used) */
    pdf->pdfLevelMaj = 1;		/* PDF level, do not use operators beyond these */
    pdf->pdfLevelMin = 1;
    pdf->defdomain_unit = POINTSPERINCH;	/* unit for default domain */
    pdf->display_rotation = 270;
    pdf->useStandardOutput = 0;		/* send output to stdout if non-zero */
    pdf->linearizeON = 0;		/* linearized PDF */
    pdf->encryptionON = 0;		/* encryption to set security options */
    pdf->cryptDict.version = 1;		/* encryption dictionary stuff */
    pdf->cryptDict.revision = 2;
    pdf->cryptDict.permissions = 0xfffffffc;	/* p. 126 of the PDF Reference v1.3 */
    pdf->cryptDict.filter = "/Standard";	/* static and fixed for now */
    strcpy(pdf->cryptDict.owner, "");
    strcpy(pdf->cryptDict.user, "");

    pdf->compressionON = 0;		/* compress stream */
    pdf->compress_command = NULL;	/* command for LZW compression */
    pdf->streamFilterList = NULL;	/* for PDF stream /Filter spec */
    pdf->launchPreview = 1;		/* launch viewer application on the output file */
    pdf->filename_set = 0;		/* flag indicating if output filename is set explicitly */
    /* pdf->fncounter = 0; */			/* filename counter for a given process */
    pdf->inTextObj = 0;			/* flag indicating within Text block between BT ET */
    /* clear viewPrefs struct */
    pdf->viewerPrefs.pageMode =0;	/* This really belongs directly to Catalog obj, but here for convenience */
    pdf->viewerPrefs.hideToolbar =0;	/* when YES, tool bar in viewer will be hidden */
    pdf->viewerPrefs.hideMenubar =0;	/* when YES, menu bar in viewer will be hidden */
    pdf->viewerPrefs.hideWindowUI=0;	/* when YES, GUI elements in doc window will be hidden */
    pdf->viewerPrefs.fitWindow=0;	/* when YES, instructs the viewer to resize the doc window to the page size */
    pdf->viewerPrefs.centerWindow=0;	/* when YES, instructs the viewer to move the doc window to the screen's center */
    pdf->viewerPrefs.pageLayout=0;	/* Specifies 1-page display or two-page columnar displays */
    pdf->viewerPrefs.nonFSPageMode=0; 	/* Specifies pageMode coming out of the full-screen display mode */
    pdf->_cpdfErrorHandler = NULL;
    pdf->defaultDomain = NULL;		/* default plot domain */
    pdf->currentDomain = NULL;		/* current plot domain */
    pdf->xcurr=0.0;
    pdf->ycurr=0.0;			/* current raw coordinate.  Used in cpdfArc.c, cpdfRawPlot.c */
    pdf->x2points=1.0;
    pdf->y2points=1.0;		/* scaling factor for current domain */
    /* double xLlog, xHlog, yLlog, yHlog; */		/* scaling factor for current domain (logarithmic) */

    pdf->numFonts = 0;			/* number of fonts used */
    pdf->numExtFonts = 0;		/* number of external fonts embedded */
    /* extFontList = NULL; */		/* linked list of external fonts */
    /* CPDFfontInfo fontInfos[NMAXFONTS]; */		/* array of font infos */
    pdf->currentFont = 0;		/* current font index (into fontInfos[]) */
    pdf->inlineImages = 0;		/* in-line image count */
    pdf->numImages = 0;
    /* CPDFimageInfo imageInfos[NMAXIMAGES]; */
    pdf->imageFlagBCI = 0;		/* bit-0 (/ImageB), bit-1 (/ImageC), bit-2 (/ImageI) */
    pdf->imageSelector = 0;		/* image index (0 .. N-1) for multi-page TIFF file */
    pdf->numOutlineEntries = 0;		/* total # of outline (book mark) entries */
    pdf->firstOLentry = NULL;		/* pointer to first outline entry */
    pdf->lastOLentry = NULL;		/* pointer to last outline entry */

    pdf->hex_string = 0;		/* Current string mode: HEX or ASCII */
    pdf->font_size = 12.0;		/* current font size and info below */
    pdf->word_spacing = 0.0;
    pdf->char_spacing = 0.0;
    pdf->text_rise = 0.0;
    pdf->horiz_scaling = 100.0;		/* text horizontal scaling in percent */
    pdf->text_leading = 0.0;
    _cpdf_resetTextCTM(pdf);
/*
    pdf->textCTM.a = 1.0;
    pdf->textCTM.b = 0.0;
    pdf->textCTM.c = 0.0;
    pdf->textCTM.d = 1.0;
    pdf->textCTM.e = 0.0;
    pdf->textCTM.f = 0.0;
*/
    pdf->textClipMode = 0;

    pdf->usePDFMemStream = 1;		/* if non-zero use memory stream for PDF generation */
    pdf->pdfMemStream = NULL;		/* memory stream for PDF file that is currently active */
    pdf->useContentMemStream = 1;	/* if non-zero use memory stream for Content */
    pdf->currentMemStream = NULL;	/* memory stream for Content that is currently active */
    pdf->currentPage =1;		/* current page number that is being drawn */
    pdf->numPages =1;			/* number of pages - may be greater than actual # of pages */
    /* CPDFpageInfo pageInfos[NMAXPAGES+1]; */	/* array of pageInfo structure for all pages */
    pdf->numKids = 0;			/* actual # of pages counted for Pages object */
    /* int kidsIndex[NMAXPAGES]; */	/* object index list for kids to be written to Pages object */
    pdf->scratchMem = NULL;		/* use this as non-overflowing scratch pad */
    pdf->fpcg = NULL; 			/* Output file */
    pdf->fpcontent = NULL;		/* Content stream (need length) */
    pdf->numAnnots = 0;			/* count of annotations */
    /* CPDFannotInfo annotInfos[NMAXANNOTS]; */	/* array of annotInfo structure for all annotations */
    /* char mediaBox[64]; */		/* MediaBox for current page*/
    /* cropBox[64]; */			/* CropBox for current page */
    pdf->currentByteCount = 0;		/* # of bytes written, or offset of next object */
    /* char creator_name[64]; */	/* Info: set it by cpdf_setCreator() */
    /* char file_title[64]; */		/* Info: title of PDF file */
    /* char file_subject[64]; */	/* Info: subject of PDF file */
    /* char file_keywords[128]; */	/* Info: keywords */
    /* char username[64]; */		/* user name */
    /* char filenamepath[1024]; */
    /* char contentfile[1024]; */
    /* long objByteOffset[NMAXOBJECTS]; */	/* offset into object number N */
    /* int  objIndex[NMAXOBJECTS]; */	/* object index for selected objects */
    pdf->startXref = 0;			/* offset of xref */
    /* char spbuf[2048]; */		/* scratch buffer for sprintf */

}


/* ---------------------------------------------------------------------------------------- */
#ifdef MacOS8

/* MacOS function to set info fork items */
/* Many thanks to James W. Walker <http://www.jwwalker.com/> for Mac related fixes. */
void SetFileInfo(const char *fileName, OSType fileType, OSType fileCreator)
{
    FInfo outFileInfo;
    Str31 nameCopy;
	strcpy( (char*)nameCopy, fileName );
	c2pstr( (char*)nameCopy );
	GetFInfo(nameCopy, 0, &outFileInfo);
	outFileInfo.fdType = fileType;
	outFileInfo.fdCreator = fileCreator;
	SetFInfo(nameCopy, 0, &outFileInfo);
}
#endif

