/* cpdfAnnot.c -- Text annotations and URL actions.
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf.
 *  or http://www.fastio.com/LICENSE.pdf


## MULTI_THREADING_OK  ##
1999-11-27 [io]
	Unicode ready.  Unocode must be passed as HEX, but written into PDF as binary.

 FIXME:  There are a lot of duplications of code here.  Consolidate and clean up.
	 [1999-04-03]
*/

#include "config.h"
#include "version.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "cpdflib.h"
#include "cpdf_private.h"
#include "cglobals.h"

/****************************************************************************/
/* This function will enable and set Text Annotation for the PDF file */
/* Only the first 64 chars of title string will be used, the remainder will be truncated */
/* Call after cpdf_open(), but before cpdf_init() */

void cpdf_setAnnotation(CPDFdoc *pdf, float xll, float yll, float xur, float yur,
	const char *title, const char *annot, CPDFannotAttrib *attr) 
{
	cpdf_rawSetAnnotation( pdf, x_Domain2Points(pdf,  xll), y_Domain2Points(pdf, yll),
			       x_Domain2Points(pdf, xur), y_Domain2Points(pdf, yur),
			       title, annot, attr);
}


void cpdf_setActionURL(CPDFdoc *pdf, float xll, float yll, float xur, float yur,
	const char *linkspec, CPDFannotAttrib *attr)
{
	cpdf_rawSetActionURL(pdf, x_Domain2Points(pdf, xll), y_Domain2Points(pdf, yll),
			      x_Domain2Points(pdf, xur), y_Domain2Points(pdf, yur),
			      linkspec, attr);
}


void cpdf_setLinkGoToPage(CPDFdoc *pdf, float xll, float yll, float xur, float yur, int page,
	const char *fitmode, CPDFannotAttrib *attr)
{
	cpdf_rawSetLinkGoToPage(pdf, x_Domain2Points(pdf, xll), y_Domain2Points(pdf, yll),
			      x_Domain2Points(pdf, xur), y_Domain2Points(pdf, yur),
			      page, fitmode, attr);
}


void cpdf_setLinkAction(CPDFdoc *pdf, float xll, float yll, float xur, float yur,
	const char *action_dict, CPDFannotAttrib *attr)
{
	cpdf_rawSetLinkAction(pdf, x_Domain2Points(pdf, xll), y_Domain2Points(pdf, yll),
			      x_Domain2Points(pdf, xur), y_Domain2Points(pdf, yur),
			      action_dict, attr);
}



int cpdf_includeTextFileAsAnnotation(CPDFdoc *pdf, float xll, float yll, float xur, float yur,
	const char *title, const char *filename, CPDFannotAttrib *attr)
{
	return( cpdf_rawIncludeTextFileAsAnnotation(pdf,
				x_Domain2Points(pdf, xll), y_Domain2Points(pdf, yll),
			      	x_Domain2Points(pdf, xur), y_Domain2Points(pdf, yur),
			      	title, filename, attr) );
}

/* Raw coordinate versions ===================================================== */

void cpdf_rawSetAnnotation(CPDFdoc *pdf, float xll, float yll, float xur, float yur,
	const char *title, const char *annot, CPDFannotAttrib *attr) 
{
int foundInPageList = 0;	/* in pageInfos[].annotIdx[] */
char *ptr, *ptr2, ch;
int  escapecount = 0, newlncount = 0;
CPDFannotInfo *aI;
	if(pdf->numAnnots >= pdf->nMaxAnnots) {
	    cpdf_Error(pdf, 1, "ClibPDF", "too many annotations in this PDF: %d", pdf->nMaxAnnots);
	    return;
	}
	aI = &pdf->annotInfos[pdf->numAnnots];
	aI->type = ANNOT_TEXT;
	aI->xLL = xll;
	aI->yLL = yll;
	aI->xUR = xur;
	aI->yUR = yur;

	/* -------- Added for ver 1.10 and later ----------- */
	aI->border_array = NULL;
	aI->BS = NULL;
	if(attr != NULL) {
	    aI->flags = attr->flags;
	    aI->r = attr->r;
	    aI->g = attr->g;
	    aI->b = attr->b;
	    if(attr->BS != NULL) {
		aI->BS = (char *)malloc((size_t)(strlen(attr->BS) + 1));
		_cpdf_malloc_check((void *)aI->BS);
		strcpy(aI->BS, attr->BS);
	    }
	    else if(attr->border_array != NULL) {
		aI->border_array = (char *)malloc((size_t)(strlen(attr->border_array) + 1));
		_cpdf_malloc_check((void *)aI->border_array);
		strcpy(aI->border_array, attr->border_array);
	    }
	    else {
	        aI->border_array = NULL;
	        aI->BS = NULL;
	    }

	}
	else {
	    /* NULL spec for detailed annotation attributes, so give defaults */
	    aI->flags = 0;
	    aI->r = 1.0;
	    aI->g = 1.0;
	    aI->b = 0.0;
	    aI->BS = NULL;
	    aI->border_array = NULL;
	}
	/* -------- Above added for ver 1.10 and later ----------- */

	if(pdf->hex_string) {
	    /* Annotation is in HEX encoding */
	    long length;
	    char *tempbuf = (char *)malloc((size_t)(strlen(annot)/2 + 3));
	    _cpdf_malloc_check((void *)tempbuf);
	    cpdf_convertHexToBinary(annot, (unsigned char *)tempbuf, &length);
	    aI->content_link = cpdf_escapeSpecialCharsBinary(tempbuf, length, &(aI->content_len));
	    free(tempbuf);

	    tempbuf = (char *)malloc((size_t)(strlen(title)/2 + 3));
	    _cpdf_malloc_check((void *)tempbuf);
	    cpdf_convertHexToBinary(title, (unsigned char *)tempbuf, &length);
	    aI->annot_title = cpdf_escapeSpecialCharsBinary(tempbuf, length, &(aI->title_len));
	    free(tempbuf);
	}
	else {
	    /* Annotation is not in HEX mode */
	    int str_len;
	    /* Scan annotation string and determine how long the escaped string will be
		by counting chars that require escaping, and # of newline chars. */
	    ptr = (char *)annot;
	    while( (ch = *ptr++) != '\0') {
		if( (ch == '(') || (ch == ')') || (ch == '\\') ) escapecount++;
		if( (ch == '\n') || (ch == '\r') ) newlncount++;
	    }
    
	    /* Allocate sufficient string space */
	    aI->content_link = (char *)malloc((size_t)(strlen(annot) + escapecount + newlncount*3 + 16));
	    _cpdf_malloc_check((void *)aI->content_link);
	    ptr = (char *)annot;
	    ptr2 = aI->content_link;
	    /* escape special chars \ and () for PS/PDF string */
	    while( (ch = *ptr++) != '\0') {
		if( (ch == '\\') || (ch == '(') || (ch == ')') ) {
		    *ptr2++ = '\\';
		    *ptr2++ = ch;
		}
		else if((ch == '\n') || (ch == '\r')) {
		    *ptr2++ = '\\';
		    *ptr2++ = 'r';
		    *ptr2++ = '\\';
		    *ptr2++ = 'n';		/* convert LF to CRLF for maximum compatibility */
		}
		else
		    *ptr2++ = ch;
	    }
	    *ptr2 = '\0';
	    aI->content_len = strlen(aI->content_link);		/* length of annoation */

	    str_len = strlen(title);
	    aI->annot_title = (char *) malloc((size_t)(str_len +1));
	    _cpdf_malloc_check((void *)aI->annot_title);
	    strcpy(aI->annot_title, title);
	    aI->title_len = str_len;
	}

	/*  Store the index of this annotation in the list for the current page in pageInfos[]. */
	{   /* braces for local variable */
	    CPDFintList *iLprev, *iL = pdf->pageInfos[pdf->currentPage].annotIdx;
	    iLprev = iL;
	    foundInPageList = 0;
	    while(iL != NULL) {
		if(iL->value == pdf->numAnnots) {
		    foundInPageList = 1;
		    break;
		}
		iLprev = iL;	/* save the last one */
		iL = iL->next;
	    }
	    if(!foundInPageList) {
		iL = (CPDFintList *)malloc((size_t)sizeof(CPDFintList));	/* alloc the next one */
		_cpdf_malloc_check((void *)iL);
		iLprev->next = iL;			/* this will be the new last one */
		iLprev->value = pdf->numAnnots;	/* store the index to current font */
		iL->next = NULL;
		iL->value = -1;
		pdf->pageInfos[pdf->currentPage].npAnnot++;
	    }
	}

	pdf->numAnnots++;
}



void cpdf_rawSetActionURL(CPDFdoc *pdf, float xll, float yll, float xur, float yur,
	const char *linkspec, CPDFannotAttrib *attr)
{
int foundInPageList = 0;
CPDFannotInfo *aI;
char *s;
	if(pdf->numAnnots >= pdf->nMaxAnnots) {
	    cpdf_Error(pdf, 1, "ClibPDF", "too many annotations in this PDF: %d", pdf->nMaxAnnots);
	    return;
	}
	aI = &pdf->annotInfos[pdf->numAnnots];
	aI->type = ANNOT_URL;
	aI->xLL = xll;
	aI->yLL = yll;
	aI->xUR = xur;
	aI->yUR = yur;

	/* -------- Added for ver 1.10 and later ----------- */
	aI->border_array = NULL;
	aI->BS = NULL;
	if(attr != NULL) {
	    aI->flags = attr->flags;
	    aI->r = attr->r;
	    aI->g = attr->g;
	    aI->b = attr->b;
	    if(attr->BS != NULL) {
		aI->BS = (char *)malloc((size_t)(strlen(attr->BS) + 1));
		_cpdf_malloc_check((void *)aI->BS);
		strcpy(aI->BS, attr->BS);
	    }
	    else if(attr->border_array != NULL) {
		aI->border_array = (char *)malloc((size_t)(strlen(attr->border_array) + 1));
		_cpdf_malloc_check((void *)aI->border_array);
		strcpy(aI->border_array, attr->border_array);
	    }
	    else {
	        aI->border_array = NULL;
	        aI->BS = NULL;
	    }

	}
	else {
	    /* NULL spec for detailed annotation attributes, so give defaults */
	    aI->flags = 0;
	    aI->r = 0.0;
	    aI->g = 0.0;
	    aI->b = 1.0;
	    aI->BS = NULL;
	    s = "[0 0 0]";
	    aI->border_array = (char *)malloc((size_t)(strlen(s) + 1));
	    _cpdf_malloc_check((void *)aI->border_array);
	    strcpy(aI->border_array, s);
	}
	/* -------- Above added for ver 1.10 and later ----------- */


	/* Allocate sufficient string space */
	aI->content_link = (char *)malloc((size_t)(strlen(linkspec) +1));
	_cpdf_malloc_check((void *)aI->content_link);
	strcpy(aI->content_link, linkspec);
	aI->content_len = strlen(aI->content_link);		/* length of annoation */
	aI->annot_title = NULL;
	aI->title_len = 0;
	/*  Store the index of this link in the list for the current page in pageInfos[]. */
	{   /* braces for local variable */
	    CPDFintList *iLprev, *iL = pdf->pageInfos[pdf->currentPage].annotIdx;
	    iLprev = iL;
	    foundInPageList = 0;
	    while(iL != NULL) {
		if(iL->value == pdf->numAnnots) {
		    foundInPageList = 1;
		    break;
		}
		iLprev = iL;	/* save the last one */
		iL = iL->next;
	    }
	    if(!foundInPageList) {
		iL = (CPDFintList *)malloc((size_t)sizeof(CPDFintList));	/* alloc the next one */
		_cpdf_malloc_check((void *)iL);
		iLprev->next = iL;			/* this will be the new last one */
		iLprev->value = pdf->numAnnots;	/* store the index to current font */
		iL->next = NULL;
		iL->value = -1;
		pdf->pageInfos[pdf->currentPage].npAnnot++;
	    }
	}

	pdf->numAnnots++;
}

/* Goto another page within the same PDF file (we are generating). */

void cpdf_rawSetLinkGoToPage(CPDFdoc *pdf, float xll, float yll, float xur, float yur, int page,
	const char *fitmode, CPDFannotAttrib *attr)
{
int foundInPageList = 0;
CPDFannotInfo *aI;

	if(pdf->numAnnots >= pdf->nMaxAnnots) {
	    cpdf_Error(pdf, 1, "ClibPDF", "too many annotations in this PDF: %d", pdf->nMaxAnnots);
	    return;
	}
	aI = &pdf->annotInfos[pdf->numAnnots];
	aI->type = ANNOT_GOTO;
	aI->xLL = xll;
	aI->yLL = yll;
	aI->xUR = xur;
	aI->yUR = yur;
	aI->page = page;

	/* -------- Added for ver 1.10 and later ----------- */
	aI->border_array = NULL;
	aI->BS = NULL;
	if(attr != NULL) {
	    aI->flags = attr->flags;
	    aI->r = attr->r;
	    aI->g = attr->g;
	    aI->b = attr->b;
	    if(attr->BS != NULL) {
		aI->BS = (char *)malloc((size_t)(strlen(attr->BS) + 1));
		_cpdf_malloc_check((void *)aI->BS);
		strcpy(aI->BS, attr->BS);
	    }
	    else if(attr->border_array != NULL) {
		aI->border_array = (char *)malloc((size_t)(strlen(attr->border_array) + 1));
		_cpdf_malloc_check((void *)aI->border_array);
		strcpy(aI->border_array, attr->border_array);
	    }
	    else {
	        aI->border_array = NULL;
	        aI->BS = NULL;
	    }

	}
	else {
	    /* NULL spec for detailed annotation attributes, so give defaults */
	    aI->flags = 0;
	    aI->r = 0.0;
	    aI->g = 0.0;
	    aI->b = 1.0;
	    aI->BS = NULL;
	    aI->border_array = NULL;
	}
	/* -------- Above added for ver 1.10 and later ----------- */


	/* Allocate sufficient string space */
	aI->content_link = (char *)malloc((size_t)(strlen(fitmode) +1));
	_cpdf_malloc_check((void *)aI->content_link);
	strcpy(aI->content_link, fitmode);
	aI->content_len = strlen(aI->content_link);		/* length of annoation */
	aI->annot_title = NULL;
	aI->title_len = 0;
	/*  Store the index of this link in the list for the current page in pageInfos[]. */
	{   /* braces for local variable */
	    CPDFintList *iLprev, *iL = pdf->pageInfos[pdf->currentPage].annotIdx;
	    iLprev = iL;
	    foundInPageList = 0;
	    while(iL != NULL) {
		if(iL->value == pdf->numAnnots) {
		    foundInPageList = 1;
		    break;
		}
		iLprev = iL;	/* save the last one */
		iL = iL->next;
	    }
	    if(!foundInPageList) {
		iL = (CPDFintList *)malloc((size_t)sizeof(CPDFintList));	/* alloc the next one */
		_cpdf_malloc_check((void *)iL);
		iLprev->next = iL;			/* this will be the new last one */
		iLprev->value = pdf->numAnnots;	/* store the index to current font */
		iL->next = NULL;
		iL->value = -1;
		pdf->pageInfos[pdf->currentPage].npAnnot++;
	    }
	}

	pdf->numAnnots++;
}


/*  This allows general link action to be taken by specifying an Action Dictionary.
    This allows opening another PDF file, and even URL action done by the function
    cpdf_rawSetActionURL() above.
*/

void cpdf_rawSetLinkAction(CPDFdoc *pdf, float xll, float yll, float xur, float yur,
	const char *action_dict, CPDFannotAttrib *attr)
{
int foundInPageList = 0;
CPDFannotInfo *aI;
char *s;
	if(pdf->numAnnots >= pdf->nMaxAnnots) {
	    cpdf_Error(pdf, 1, "ClibPDF", "too many annotations in this PDF: %d", pdf->nMaxAnnots);
	    return;
	}
	aI = &pdf->annotInfos[pdf->numAnnots];
	aI->type = ANNOT_ACTION;
	aI->xLL = xll;
	aI->yLL = yll;
	aI->xUR = xur;
	aI->yUR = yur;

	/* -------- Added for ver 1.10 and later ----------- */
	aI->border_array = NULL;
	aI->BS = NULL;
	if(attr != NULL) {
	    aI->flags = attr->flags;
	    aI->r = attr->r;
	    aI->g = attr->g;
	    aI->b = attr->b;
	    if(attr->BS != NULL) {
		aI->BS = (char *)malloc((size_t)(strlen(attr->BS) + 1));
		_cpdf_malloc_check((void *)aI->BS);
		strcpy(aI->BS, attr->BS);
	    }
	    else if(attr->border_array != NULL) {
		aI->border_array = (char *)malloc((size_t)(strlen(attr->border_array) + 1));
		_cpdf_malloc_check((void *)aI->border_array);
		strcpy(aI->border_array, attr->border_array);
	    }
	    else {
	        aI->border_array = NULL;
	        aI->BS = NULL;
	    }

	}
	else {
	    /* NULL spec for detailed annotation attributes, so give defaults */
	    aI->flags = 0;
	    aI->r = 0.0;
	    aI->g = 0.0;
	    aI->b = 1.0;
	    aI->BS = NULL;
	    s = "[0 0 0]";
	    aI->border_array = (char *)malloc((size_t)(strlen(s) + 1));
	    _cpdf_malloc_check((void *)aI->border_array);
	    strcpy(aI->border_array, s);
	}
	/* -------- Above added for ver 1.10 and later ----------- */


	/* Allocate sufficient string space */
	aI->content_link = (char *)malloc((size_t)(strlen(action_dict) +1));
	_cpdf_malloc_check((void *)aI->content_link);
	strcpy(aI->content_link, action_dict);
	aI->content_len = strlen(aI->content_link);		/* length of annoation */
	aI->annot_title = NULL;
	aI->title_len = 0;
	/*  Store the index of this link in the list for the current page in pageInfos[]. */
	{   /* braces for local variable */
	    CPDFintList *iLprev, *iL = pdf->pageInfos[pdf->currentPage].annotIdx;
	    iLprev = iL;
	    foundInPageList = 0;
	    while(iL != NULL) {
		if(iL->value == pdf->numAnnots) {
		    foundInPageList = 1;
		    break;
		}
		iLprev = iL;	/* save the last one */
		iL = iL->next;
	    }
	    if(!foundInPageList) {
		iL = (CPDFintList *)malloc((size_t)sizeof(CPDFintList));	/* alloc the next one */
		_cpdf_malloc_check((void *)iL);
		iLprev->next = iL;			/* this will be the new last one */
		iLprev->value = pdf->numAnnots;	/* store the index to current font */
		iL->next = NULL;
		iL->value = -1;
		pdf->pageInfos[pdf->currentPage].npAnnot++;
	    }
	}

	pdf->numAnnots++;
}


int cpdf_rawIncludeTextFileAsAnnotation(CPDFdoc *pdf, float xll, float yll, float xur, float yur,
	const char *title, const char *filename, CPDFannotAttrib *attr)
{
int ch;
FILE *fp;
CPDFmemStream *memStream;
char *memBuffer;
int memLen, bufSize;
	if((fp = fopen(filename, "r")) == NULL) {
	    cpdf_Error(pdf, 1, "ClibPDF", "cannot open annotation text file: %s", filename);
	    return(-1);
	}
	memStream = cpdf_openMemoryStream();
	/* This thing should not hurt HEX mode file */
	while( (ch = fgetc(fp)) != EOF) {
	    if(ch == '\\' || ch == '(' || ch == ')')
		cpdf_memPutc('\\', memStream);		/* escape special chars */
            cpdf_memPutc(ch, memStream);
	}
	fclose(fp);
	cpdf_getMemoryBuffer(memStream, &memBuffer, &memLen, &bufSize);
	cpdf_rawSetAnnotation(pdf, xll, yll, xur, yur, title, memBuffer, attr);
	cpdf_closeMemoryStream(memStream);
	return(0);
}


int  _cpdf_freeAllAnnotInfos(CPDFdoc *pdf)
{
int i;
CPDFannotInfo *aI;
	for(i=0; i< pdf->numAnnots; i++) {
	    aI = &pdf->annotInfos[i];
	    if(aI->content_link)
		free(aI->content_link);
	    if(aI->annot_title)
		free(aI->annot_title);
	    if(aI->border_array)
		free(aI->border_array);
	    if(aI->BS)
		free(aI->BS);
	}
	return(0);
}

