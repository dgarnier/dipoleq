/* cpdfOutline.c
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

1999-11-27 [io]
	Unicode ready.  Unocode must be passed as HEX, but written into PDF as binary.

1999-10-15 [IO] 
	Added DEST_Y to allow specification of Y position only, keeping all others.

1999-08-17 [IO] 
	## MULTI_THREADING_OK ##

1999-04-03 [IO]
	Support for Action dictionary in addition to goto a page.

1998-11-26 [IO]
	Multi-level Outline (Book Mark) support.
	The only limitation of this version is that outline entries must be added
	in the top-to-bottom order, (i.e., in the order they will appear in the outline view).
	This version does not support insertion of an outline entry in the middle of an outline
	tree, or deletion of entries.  This will be fixed at some point.  It is a matter of
	adding a few more functions that perform editing of a linked list.
*/


/* #define CPDF_DEBUG 1 */

#include "config.h"
#include "version.h"

#include <string.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "cpdflib.h"		/* This must be included before all other local include files */
#include "cpdf_private.h"
#include "cglobals.h"

/* Caution:
   For now, only adding at the end is allowed.  I.e., you can't insert in the middle
   of an existing outline tree.

-------------------------------------------------------------------------------------------------------
## PARAMETER DEFINITIONS ##

   "afterThis" is another outline entry below which the new entry will be added.
	For the first every outline entry for a document use afterThis = NULL.
   If "sublevel" is non-zero (OL_SUBENT), the new entry will be added as a sub-entry under afterThis.
	If it is zero (OL_SAME), the new entry will be added at the same level as "afterThis."
   If "open" is non-zero (OL_OPEN), subentries under this entry will be open and visible.
	If it is zero (OL_CLOSED), subentries will NOT be visible upon opening the document
	(you can manually open).
   "page" is the destination page number for this outline entry as specified in
	cpdf_pageInit().
   "title" is a string for the outline.  It will automatically be escaped of special PDF chars.
   "mode" should be one of:
	DEST_NULL	0	/XYZ null null null (keep current position and zoom)
	DEST_XYZ	1	/XYZ left top zoom (zoom=1.0 is normal size) 
	DEST_FIT	2	/Fit
	DEST_FITH	3	/FitH top
	DEST_FITV	4	/FitV left
	DEST_FITR	5	/FitR left bottom right top
	DEST_FITB	6	/FitB		(fit bounding box to window) PDF-1.1
	DEST_FITBH	7	/FitBH top	(fit width of bounding box to window) PDF-1.1
	DEST_FITBV	8	/FitBV left	(fit height of bounding box to window) PDF-1.1
        DEST_Y		9	/XYZ null y null

    "p1 .. p4" are parameters for various "modes" above.  Not all parameters are used.
	pass a dummy float variable (the value doesn't matter) for those parameters unused for
	a given "mode."
-------------------------------------------------------------------------------------------------------
*/

CPDFoutlineEntry *cpdf_addOutlineEntry(CPDFdoc *pdf, CPDFoutlineEntry *afterThis, int sublevel,
	int open, int page, const char *title, int mode, float p1, float p2, float p3, float p4)
{
CPDFoutlineEntry *newOL, *parentptr;
int openp=1;
    if(page > pdf->numPages) {
	cpdf_Error(pdf, 1, "ClibPDF", "Page %d has not yet been created for outline (book mark).", page);
	return(NULL);
    }
    newOL = (CPDFoutlineEntry *)malloc((size_t)sizeof(CPDFoutlineEntry));
    _cpdf_malloc_check((void *)newOL);
    newOL->dest_page = page;
    newOL->open = open;
    newOL->dest_attr_act_dict = _cpdf_dest_attribute(mode, p1, p2, p3, p4);

    if(pdf->hex_string) {
	/* Outline title is in HEX encoding */
	long length;
	char *tempbuf = (char *)malloc((size_t)(strlen(title)/2 + 3));
	_cpdf_malloc_check((void *)tempbuf);
	cpdf_convertHexToBinary(title, (unsigned char *)tempbuf, &length);
	newOL->title = cpdf_escapeSpecialCharsBinary(tempbuf, length, &(newOL->title_len));
	free(tempbuf);
    }
    else {
	/* Outline title is in standard ASCII */
        newOL->title = cpdf_escapeSpecialChars(title);	/* this function allocates storage */
	newOL->title_len = strlen(newOL->title);;
    }

    if(sublevel && afterThis != NULL) {
	/* Add the new entry at one level below as "afterThis", i.e., as a sub-section. */
	newOL->parent = afterThis;		/* that is our parent */
	newOL->prev = NULL;
	newOL->first = NULL;			/* first, last are only needed if kids are added */
	newOL->last = NULL;
	newOL->next = NULL;
	newOL->count = 0;			/* will become non-zero only if kids are added later */
        /* now adding a sub-section requires modification to parent entry */
	if(afterThis->count == 0)
	    afterThis->first = newOL;		/* this is the first kid of afterThis */
	afterThis->last = newOL;
	/* Increment "count" of all open ancestors */
	parentptr = newOL->parent;			/* initialize to my parent */
	if(parentptr == NULL)
	    pdf->numOutlineEntries++;		/* incremeent global outline entry count */
	else if(parentptr->open)
	    pdf->numOutlineEntries++;		/* incremeent global outline entry count */
	while((parentptr != NULL) && openp) {
	    parentptr->count++;			/* increment kids count of ancestors */
	    openp = parentptr->open;		/* if parent is open, then continue on */
	    parentptr = parentptr->parent;
	}
    }
    else {
	/* Add the new entry at the same level as "afterThis" */
	if(afterThis == NULL) {
	    /* This is the first ever outline entry for this document. */
	    newOL->parent = NULL;		/* eventuall will point to global Outlines object */
	    newOL->prev = NULL;
	    newOL->first = NULL;		/* first, last are only needed if kids are added */
	    newOL->last = NULL;
	    newOL->next = NULL;
	    newOL->count = 0;			/* will become non-zero only if kids are added later */
	    /* These 3 are for the global Outlines object */
	    pdf->numOutlineEntries++;		/* incremeent global outline entry count */
	    pdf->firstOLentry = newOL;
	    pdf->lastOLentry = newOL;
	}
	else {
	    /* We are adding an entry at the same level as "afterThis" */
	    newOL->parent = afterThis->parent;	/* common parent */
	    newOL->prev = afterThis;
	    newOL->first = NULL;		/* first, last are only needed if kids are added */
	    newOL->last = NULL;
	    newOL->next = NULL;
	    newOL->count = 0;			/* will become non-zero only if kids are added later */
	    /* new, set fields of the previous entry */
	    afterThis->next = newOL;
	    if(afterThis->parent)
		afterThis->parent->last = newOL;

	    /* Increment "count" of all open ancestors */
	    parentptr = newOL->parent;			/* initialize to my parent */
	    if(parentptr == NULL)
		pdf->numOutlineEntries++;		/* incremeent global outline entry count */
	    else if(parentptr->open)
		pdf->numOutlineEntries++;		/* incremeent global outline entry count */
	    while((parentptr != NULL) && openp) {
		parentptr->count++;			/* increment kids count of ancestors */
		openp = parentptr->open;		/* if parent is open, then continue on */
		parentptr = parentptr->parent;
	    }
	}
    }
 
    if(newOL->parent == NULL)
        pdf->lastOLentry = newOL;		/* if our direct parent is the global Outline object */
    return(newOL);
}

/* FIXME: This is nearly a copy of the above cpdf_addOutlineEntry() to quickly add Action support.
   Remove duplication and consolidate at some point.
*/
CPDFoutlineEntry *cpdf_addOutlineAction(CPDFdoc *pdf, CPDFoutlineEntry *afterThis, int sublevel, int open,
		const char *action_dict, const char *title)
{
CPDFoutlineEntry *newOL, *parentptr;
int openp=1, page= -1;		/* nagative page indicates action rather than page destination. */
    newOL = (CPDFoutlineEntry *)malloc((size_t)sizeof(CPDFoutlineEntry));
    _cpdf_malloc_check((void *)newOL);
    newOL->dest_page = page;
    newOL->open = open;
    newOL->dest_attr_act_dict = _cpdf_copy_string_with_malloc(action_dict);

    if(pdf->hex_string) {
	/* Outline title is in HEX encoding */
	long length;
	char *tempbuf = (char *)malloc((size_t)(strlen(title)/2 + 3));
	_cpdf_malloc_check((void *)tempbuf);
	cpdf_convertHexToBinary(title, (unsigned char *)tempbuf, &length);
	newOL->title = cpdf_escapeSpecialCharsBinary(tempbuf, length, &(newOL->title_len));
	free(tempbuf);
    }
    else {
	/* Outline title is in standard ASCII */
        newOL->title = cpdf_escapeSpecialChars(title);	/* this function allocates storage */
	newOL->title_len = strlen(newOL->title);;
    }

    if(sublevel && afterThis != NULL) {
	/* Add the new entry at one level below as "afterThis", i.e., as a sub-section. */
	newOL->parent = afterThis;		/* that is our parent */
	newOL->prev = NULL;
	newOL->first = NULL;			/* first, last are only needed if kids are added */
	newOL->last = NULL;
	newOL->next = NULL;
	newOL->count = 0;			/* will become non-zero only if kids are added later */
        /* now adding a sub-section requires modification to parent entry */
	if(afterThis->count == 0)
	    afterThis->first = newOL;		/* this is the first kid of afterThis */
	afterThis->last = newOL;
	/* Increment "count" of all open ancestors */
	parentptr = newOL->parent;			/* initialize to my parent */
	if(parentptr == NULL)
	    pdf->numOutlineEntries++;		/* incremeent global outline entry count */
	else if(parentptr->open)
	    pdf->numOutlineEntries++;		/* incremeent global outline entry count */
	while((parentptr != NULL) && openp) {
	    parentptr->count++;			/* increment kids count of ancestors */
	    openp = parentptr->open;		/* if parent is open, then continue on */
	    parentptr = parentptr->parent;
	}
    }
    else {
	/* Add the new entry at the same level as "afterThis" */
	if(afterThis == NULL) {
	    /* This is the first ever outline entry for this document. */
	    newOL->parent = NULL;		/* eventuall will point to global Outlines object */
	    newOL->prev = NULL;
	    newOL->first = NULL;		/* first, last are only needed if kids are added */
	    newOL->last = NULL;
	    newOL->next = NULL;
	    newOL->count = 0;			/* will become non-zero only if kids are added later */
	    /* These 3 are for the global Outlines object */
	    pdf->numOutlineEntries++;		/* incremeent global outline entry count */
	    pdf->firstOLentry = newOL;
	    pdf->lastOLentry = newOL;
	}
	else {
	    /* We are adding an entry at the same level as "afterThis" */
	    newOL->parent = afterThis->parent;	/* common parent */
	    newOL->prev = afterThis;
	    newOL->first = NULL;		/* first, last are only needed if kids are added */
	    newOL->last = NULL;
	    newOL->next = NULL;
	    newOL->count = 0;			/* will become non-zero only if kids are added later */
	    /* new, set fields of the previous entry */
	    afterThis->next = newOL;
	    if(afterThis->parent)
		afterThis->parent->last = newOL;

	    /* Increment "count" of all open ancestors */
	    parentptr = newOL->parent;			/* initialize to my parent */
	    if(parentptr == NULL)
		pdf->numOutlineEntries++;		/* incremeent global outline entry count */
	    else if(parentptr->open)
		pdf->numOutlineEntries++;		/* incremeent global outline entry count */
	    while((parentptr != NULL) && openp) {
		parentptr->count++;			/* increment kids count of ancestors */
		openp = parentptr->open;		/* if parent is open, then continue on */
		parentptr = parentptr->parent;
	    }
	}
    }
 
    if(newOL->parent == NULL)
        pdf->lastOLentry = newOL;		/* if our direct parent is the global Outline object */
    return(newOL);
}



char *_cpdf_copy_string_with_malloc(const char *action_dict)
{
char *sptr;
    sptr = (char *)malloc((size_t)(strlen(action_dict)+1));
    _cpdf_malloc_check((void *)(sptr));
    strcpy(sptr, action_dict);
    return(sptr);
}

/* This function will return a malloced memory, which must eventually be freed.*/
char *_cpdf_dest_attribute(int mode, float p1, float p2, float p3, float p4)
{
char buff[128], *sptr;
    switch(mode) {
	default:
	case DEST_NULL:	/* keep current display location and zoom */
			sprintf(buff, "/XYZ null null null");
			break;
	case DEST_XYZ:	/* 1 - /XYZ left top zoom (equivalent to above; zoom=1.0 is normal size) */
			sprintf(buff, "/XYZ %.0f %.0f %.3f", p1, p2, p3);
			break;
	case DEST_FIT:	/* 2 - /Fit  */
			sprintf(buff, "/Fit");
			break;
	case DEST_FITH:	/* 3 - /FitH top */
			sprintf(buff, "/FitH %.0f", p1);
			break;
	case DEST_FITV:	/* 4 - /FitV left */
			sprintf(buff, "/FitV %.0f", p1);
			break;
	case DEST_FITR:	/* 5 - /FitR left bottom right top */
			sprintf(buff, "/FitR %.0f %.0f %.0f %.0f", p1, p2, p3, p4);
			break;
	case DEST_FITB:	/* 6 - /FitB   (fit bounding box to window) PDF-1.1 */
			sprintf(buff, "/FitB");
			break;
	case DEST_FITBH: /* 7 - /FitBH top   (fit width of bounding box to window) PDF-1.1 */
			sprintf(buff, "/FitBH %.0f", p1);
			break;
	case DEST_FITBV: /* 8 - /FitBV left   (fit height of bounding box to window) PDF-1.1 */
			sprintf(buff, "/FitBV %.0f", p1);
			break;
	case DEST_Y:	 /* 9 - /XYZ null y null (positions only Y, leave others unchanged) */
			sprintf(buff, "/XYZ null %.0f null", p1);
			break;
    }
    sptr = (char *)malloc((size_t)(strlen(buff)+1));
    _cpdf_malloc_check((void *)(sptr));
    strcpy(sptr, buff);
    return(sptr);
}


/* Call in cpdfInit.c: cpdf_finalizeAll() as
   _cpdf_serializeOutlineEntries(&objectCount, firstOLentry, lastOLentry);
*/

void _cpdf_serializeOutlineEntries(int *objCount, CPDFoutlineEntry *first, CPDFoutlineEntry *last)
{
CPDFoutlineEntry *olent;
    olent = first;
    do {
	olent->objIndex = (*objCount)++;	/* set serialized index, and increment count */
	/* fprintf(stderr,"OutlineEnt: page=%d, object=%d (%s)\n",
		olent->dest_page, olent->objIndex, olent->title);
	*/
	/* serialize kids first */
	if(olent->first != NULL)
	    _cpdf_serializeOutlineEntries(objCount, olent->first, olent->last);	/* recurse */
	olent = olent->next;		/* do the next one on the same level */
    } while(olent != NULL);
}

void _cpdf_WriteOutlineEntries(CPDFdoc *pdf, CPDFoutlineEntry *first, CPDFoutlineEntry *last)
{
CPDFoutlineEntry *olent;
    olent = first;
    do {
	_cpdf_WriteOneOutlineEntry(pdf, olent);	/* write this one */
	/* write kids */
	if(olent->first != NULL)
	    _cpdf_WriteOutlineEntries(pdf, olent->first, olent->last);	/* recurse */
	olent = olent->next;		/* do the next one on the same level */
    } while(olent != NULL);
}


long _cpdf_WriteOneOutlineEntry(CPDFdoc *pdf, CPDFoutlineEntry *olent)
{
	sprintf(pdf->spbuf, "%d 0 obj\n", olent->objIndex);		_cpdf_pdfWrite(pdf, pdf->spbuf);
	sprintf(pdf->spbuf, "<<\n");					_cpdf_pdfWrite(pdf, pdf->spbuf);
	/* NULL parent means the parent is the global Outlines object */
	if(olent->parent)
	    sprintf(pdf->spbuf, "/Parent %d 0 R\n", olent->parent->objIndex);
	else
	    sprintf(pdf->spbuf, "/Parent %d 0 R\n", pdf->objIndex[CPDF_Outlines]);
	_cpdf_pdfWrite(pdf, pdf->spbuf);

	if(olent->dest_page >= 0) {
	    /* goto a specified page in the current document */
	    sprintf(pdf->spbuf, "/Dest [%d 0 R %s]\n",
		pdf->pageInfos[olent->dest_page].objIndex, olent->dest_attr_act_dict);
	    _cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	else {
	    /* Action Outline, we just copy content of dictionary as this is quite complex.
		No checking is performed on the correctness of syntax.
	    */
	    sprintf(pdf->spbuf, "/A <<\n%s\n>>\n", olent->dest_attr_act_dict);
	    _cpdf_pdfWrite(pdf, pdf->spbuf);
	}

	_cpdf_pdfWrite(pdf, "/Title (");
	_cpdf_pdfBinaryWrite(pdf, olent->title, olent->title_len);	/* Unicode ready */
	_cpdf_pdfWrite(pdf, ")\n");

	/* Now, write out, first, last, next, prev if they are not NULL */
	if(olent->next) {
	    sprintf(pdf->spbuf, "/Next %d 0 R\n", olent->next->objIndex);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	if(olent->prev) {
	    sprintf(pdf->spbuf, "/Prev %d 0 R\n", olent->prev->objIndex);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	if(olent->first) {
	    sprintf(pdf->spbuf, "/First %d 0 R\n", olent->first->objIndex);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	if(olent->last) {
	    sprintf(pdf->spbuf, "/Last %d 0 R\n", olent->last->objIndex);	_cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	if(olent->count != 0) {
	    if(olent->open)
	        sprintf(pdf->spbuf, "/Count %d\n", olent->count);
	    else
	        sprintf(pdf->spbuf, "/Count -%d\n", olent->count);
	    _cpdf_pdfWrite(pdf, pdf->spbuf);
	}
	_cpdf_pdfWrite(pdf, ">>\n");
	_cpdf_pdfWrite(pdf, "endobj\n");
	pdf->objByteOffset[olent->objIndex + 1] = pdf->currentByteCount;
	return(pdf->currentByteCount);
}

/* static int freecnt = 0; */

/* Free entries from the end */
void _cpdf_freeAllOutlineEntries(CPDFoutlineEntry *firstOL, CPDFoutlineEntry *lastOL)
{
CPDFoutlineEntry *olent, *olent2;

    olent = firstOL;
    do {
	/* free kids first */
	if(olent->first != NULL)
	    _cpdf_freeAllOutlineEntries(olent->first, olent->last);	/* recurse */
	if(olent->title) free(olent->title);
	if(olent->dest_attr_act_dict) free(olent->dest_attr_act_dict);
	olent2 = olent->next;	/* save the next one at the same level */
	free(olent);
	olent = olent2;
	/* freecnt++;
	fprintf(stderr,"freecount=%d\n", freecnt);
	*/
    } while(olent != NULL);
}



