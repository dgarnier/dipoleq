/*  ClibPDF library globals -- This file must be imported after cpdflib.h in each module.
 * Copyright (C) 1998-1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf.
 *  or http://www.fastio.com/LICENSE.pdf

 1999-08-12 [IO]
	Multi-threaded version.  Most of stuff here is for backward compatibility or for
	really constant data.
*/


/* MAINDEF is defined in file cpdfInit.c. In all other files, it is not defined. */
#if defined(MAINDEF)
    int global_docID = 0;		/* storage for doc ID serial number */

    int g_nMaxPages = NMAXPAGES;
    int g_nMaxFonts = NMAXFONTS;
    int g_nMaxImages = NMAXIMAGES;
    int g_nMaxAnnots = NMAXANNOTS;
    int g_nMaxObjects = NMAXOBJECTS;
    /* CPDFdoc *g_pdf_doc; */		/* Global PDF document struct for backward compatibility */
    /* Don't change these, use cpdf_setMonthNames(char *mnArray[]) for other languages. */
    char *monthNameEnglish[] = { "Jan", "Feb", "Mar", "Apr", "May", "Jun",
			     	 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
    CPDFglobalErrorHandler _cpdfGlobalErrorHandler = NULL;	/* global error handler */
    
#else
    extern int global_docID;		/* storage for doc ID serial number */

    extern int g_nMaxPages;
    extern int g_nMaxFonts;
    extern int g_nMaxImages;
    extern int g_nMaxAnnots;
    extern int g_nMaxObjects;
    /* extern CPDFdoc *g_pdf_doc; */		/* Global PDF document struct for backward compatibility */
    extern char *monthNameEnglish[];
    extern CPDFglobalErrorHandler _cpdfGlobalErrorHandler;	/* global error handler */


#endif


