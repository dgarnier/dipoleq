/* outline.c -- Test program for multi-level outline (book mark) in a large PDF file (200 pages).
 * Copyright (C) 1998,1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

NEXTSTEP/OPENSTEP
cc -Wall -o outline outline.c -I/usr/local/include -lcpdfpm

BSDI/gcc
gcc -Wall -o outline outline.c -L. -lcpdfpm -lm

SunOS5.4/CC-4.0
cc -o outline outline.c -Xa -L. -lcpdfpm -lm

1999-08-24  [io] for v2.00 API
1999-04-05  [IO] Mods for adding cpdf_addOutlineAction().
1998-11-27  [IO] Begin
*/

#define	NPAGES		200		/* number of pages to generate */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cpdflib.h"

void show_textStuff(CPDFdoc *pdf, float yloc, int page, int chapter, int section, int subsection)
{
float ybump = 0.3;
char strbuf[128];
    /* == Text examples == */
    /* cpdf_setgrayFill(pdf, 0.0); */				/* Black */
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Times-Roman", "MacRomanEncoding", 24.0);
    sprintf(strbuf, "Page: %d", page);
    cpdf_text(pdf, 1.0, yloc, 0.0, strbuf);	/* cpdf_text() may be repeatedly used */
    yloc -= ybump;
    sprintf(strbuf, "Chapter: %d", chapter);
    cpdf_text(pdf, 1.0, yloc, 0.0, strbuf);	/* cpdf_text() may be repeatedly used */
    yloc -= ybump;
    sprintf(strbuf, "Section: %d.%d", chapter, section);
    cpdf_text(pdf, 1.0, yloc, 0.0, strbuf);	/* cpdf_text() may be repeatedly used */
    yloc -= ybump;
    sprintf(strbuf, "Subsection: %d.%d.%d", chapter, section, subsection);
    cpdf_text(pdf, 1.0, yloc, 0.0, strbuf);	/* cpdf_text() may be repeatedly used */
    yloc -= ybump;

    cpdf_endText(pdf);
}

int main(int argc, char *argv[])
{
CPDFdoc *pdf;
int page;
int chapter, lastchap=0;
int section, lastsect=0;
int subsection, lastsubsect=0;
int same_sub = 0;
CPDFoutlineEntry *lastChapOL, *lastSectOL, *lastSubsectOL;
float p1=0.0, p2=0.0, p3=0.0, p4=0.0;
char strbuf[128], *act_dict;
CPDFdocLimits dL = {NPAGES+50, -1, -1, -1, 7000};
CPDFviewerPrefs vP = { PM_OUTLINES, 0, 0, 0, 0, 0, 0, 0 };

    /* == Initialization == */
    /* cpdf_setDocumentLimits(int maxPages, int maxFonts, int maxImages, int maxAnnots, int maxObjects) */
    /* Default limit is 100 pages, so we increase it for our large document */
    /* OBSOLETE: cpdf_setDocumentLimits(NPAGES+50, -1, -1, -1, 7000); */
    pdf = cpdf_open(0, &dL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_setOutputFilename(pdf, "outline.pdf");	/* specifying it here saves memory */
    cpdf_init(pdf);

    /* This will have outline (book marks) visible upon document open */
    /* API changed: cpdf_setViewerPreferences(PM_OUTLINES, 0, 0, 0, 0, 0, 0, 0); */
    cpdf_setViewerPreferences(pdf, &vP);

    /* Add 3 action outlines using new Action Dictionary API (ClibPDF v1.10) */
    /* ### NOTE: The very FIRST cpdf_addOutline____() must use NULL as its first argument */
    act_dict = "/Type /Action /S /URI /URI (http://www.fastio.com/example/arctest.pdf)";
    lastChapOL = cpdf_addOutlineAction(pdf, NULL, OL_SAME, OL_OPEN, act_dict, "FastIO Arc Example (URI)");

    /* page number for /GoToR starts at 0 */
    act_dict = "/Type /Action /S /GoToR /D [0 /XYZ null null null] /F (./domain.pdf) /NewWindow true";
    lastChapOL = cpdf_addOutlineAction(pdf, lastChapOL, OL_SAME, OL_OPEN, act_dict, "Domain Example (local)");

    act_dict = "/Type /Action /S /URI /URI (http://partners.adobe.com/asn/developer/PDFS/TN/PDFSPEC.PDF)";
    lastChapOL = cpdf_addOutlineAction(pdf, lastChapOL, OL_SAME, OL_OPEN, act_dict, "PDF Reference v1.3 at Adobe");


    /* == Generate many pages, and add outline (book mark) entries */
    for(page=1; page <= NPAGES; page++) {
	fprintf(stderr, "Page %d  \r", page);  fflush(stderr);
        cpdf_pageInit(pdf, page, PORTRAIT, LETTER, LETTER);	/* create page */
	chapter = (page-1) / 20 + 1;			/* 10 chapters: new chapter every 20 pages */
	section = ((page-1) % 20) / 4 + 1;		/* 5 sections per chapter: new section every 4 pages */
	subsection = ((page-1) % 4) * 2 + 1;		/* 8 subsections per section: 2 subsections per page */

	/* For this example, we have just minimal page content */
	show_textStuff(pdf, 9.5, page, chapter, section, subsection);	/* draw text on page */

	p1 = 0.0;		/* left */
	p2 = 11.0*72.0;		/* top */
	p3 = 0.0;		/* zoom (== 0 means keep the current zoom) */

	if(chapter != lastchap) {
	    /* new chapter */
	    sprintf(strbuf, "Chapter %d", chapter);
	    lastChapOL = cpdf_addOutlineEntry(pdf, lastChapOL, OL_SAME, OL_OPEN, page, strbuf,
						DEST_XYZ, p1, p2, p3, p4);
	    lastchap = chapter;
	    lastSectOL = lastChapOL;	/* first section of a given chapter is added to chapter outline */
	}
	if(section != lastsect) {
	    /* new section */
	    sprintf(strbuf, "Section %d.%d", chapter, section);
	    if(section == 1) same_sub = OL_SUBENT;	/* sect 1 is a subentry of chapter */
	    else             same_sub = OL_SAME;	/* otherwise, add it to prev section at the same level */
	    lastSectOL = cpdf_addOutlineEntry(pdf, lastSectOL, same_sub, OL_CLOSED, page, strbuf,
						DEST_XYZ, p1, p2, p3, p4);
	    lastsect = section;
	    lastSubsectOL = lastSectOL;	/* first subsection of a given section is added to section outline */
	}
	if(subsection != lastsubsect) {
	    /* new subsection */
	    sprintf(strbuf, "Subsection %d.%d.%d", chapter, section, subsection);
	    if(subsection == 1) same_sub = OL_SUBENT;	/* subsect 1 is a subentry of section */
	    else                same_sub = OL_SAME;	/* otherwise, add it to prev subsection at the same level */
	    /* All subsections for this section are closed upon display.  OL_CLOSED does that. */
	    lastSubsectOL = cpdf_addOutlineEntry(pdf, lastSubsectOL, same_sub, OL_CLOSED, page, strbuf,
						DEST_XYZ, p1, p2, p3, p4);
	    lastsubsect = subsection;
	}

	/* Another section on the bottom half of page, so we have 2 subsections per page */
	subsection++;
	show_textStuff(pdf, 4.5, page, chapter, section, subsection);	/* draw text on page */

	/* The second subsection is half way down the page.  So, when the outline for
	   this subsection is clicked and the view is smaller than the page, it will
	   show the lower half of the page by specifying p2 (top) as half way down, and
	   using "/XYZ left top zoom" destination specification.
	   There are many destination options,  For example, the destination can be
	   "/FitR left top right bottom" which will fit this rectangle within the view.
	*/
	p2 = 5.0*72.0;		/* top is now moved down half way */
	if(subsection != lastsubsect) {
	    /* new subsection */
	    sprintf(strbuf, "Subsection %d.%d.%d", chapter, section, subsection);
	    if(subsection == 1) same_sub = OL_SUBENT;	/* subsect 1 is a subentry of section */
	    else                same_sub = OL_SAME;	/* otherwise, add it to prev subsection at the same level */
	    lastSubsectOL = cpdf_addOutlineEntry(pdf, lastSubsectOL, same_sub, OL_OPEN, page, strbuf,
						DEST_XYZ, p1, p2, p3, p4);
	    lastsubsect = subsection;
	}

	/* Calling cpdf_finalizePage() will close the page and compress the page
	   content if compression is enabled.  This will conserve memory and will reduce
	   VM swapping if your PDF document has a large number of pages.
	*/
	cpdf_finalizePage(pdf, page);

    } /* end of for(page....) */

    /* == generate output and clean up == */
    fprintf(stderr, "\nNow writing PDF...\n");
    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */

    /* == Clean up == */
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
    cpdf_close(pdf);			/* shut down */
    return(0);
}

