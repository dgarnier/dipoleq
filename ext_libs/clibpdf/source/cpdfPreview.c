/* cpdfPreview.c -- launch PDF/Acrobat viewer on the output file.
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

    This module is highly platform dependent.  Please modify the function
    cpdf_openPDFinViewer() in this file and "config.h" to make it work on
    a given platform.

    If you have added PDF-opening capability for a platform not supported
    in the original distribution, please submit your changes, i.e.,
	* Modified cpdfPreview.c, config.h (any other files changed)
	* OS/platform name and version
	* Compiler/Development Environment name and version.
	* Name of PDF viewer application and version
	* Your name and e-mail etc. for acknowledgement

============ MacOS Programmers -- Do you know the answer to the following?  =========

	Now, programmatically opening a PDF file works fine for MacOS.
	I still have this question.
	
	[1] I was told that SIZE resource (see cpdfMacOS.r) must be linked in,
	      but it seems to work without. Is that really needed, if so,
	      where and how do I add it?  What's Rez ?

1999-08-17 [IO] 
	## MULTI_THREADING_OK ##

1998-10-16  FastIO Systems, [IO]
-----------------------------------------------------------------------------------------------
*/

#include "config.h"
#include "version.h"

#include <stdio.h>
#include <string.h>
#if defined(_WIN32) || defined(WIN32)
    #include <direct.h>		/* for _getcwd() */
    #include <windows.h>
#endif

#ifdef MacOS8
    #include <Types.h>
    #include <QuickDraw.h>
    /* #include <Fonts.h> */
    /* #include <Windows.h> */
    /* #include <TextEdit.h> */
    /* #include <Dialogs.h> */
    /* #include <Menus.h> */
    #include <AppleEvents.h>
    #include <StandardFile.h>
    #include <Files.h>
    #include <Errors.h>
    #include <Strings.h>	/* c2pstr() */
    /* #include <Aliases.h> */
    typedef void (*MapCreatorProc)(FSSpec *document, OSType type, OSType *creator);
    extern OSErr LaunchTheDocument(FSSpec *document, Boolean use_remote_apps, MapCreatorProc remap);
    extern void SetFrontCreator(OSType  inProcessSignature);
#endif

#include "cpdflib.h"
#include "cpdf_private.h"
#include "cglobals.h"

/****************************************************************************/

int cpdf_openPDFfileInViewer(CPDFdoc *pdf, const char *pdffilepath)
{
#ifdef MacOS8
FSSpec fsSpec;
OSErr  err;
OSType creator='CARO';
static int firstTime = 1;
#endif
char strbuf[1024], *filep;
extern int system(const char *str);

    if(pdffilepath != NULL)
	filep = (char *)pdffilepath;	/* file is specified */
    else {
	/* NULL has been passed */
	filep = pdf->filenamepath;	/* PDF file just saved. */
	if(pdf->useStandardOutput || !pdf->launchPreview)
	    return(-1);
    }

/* ----- Platform-specific section --------------------------------------------------- */

#if defined(BSDI) || defined(FreeBSD)  || defined(Linux) || defined(SunOS5x) || defined(SolarisGCC)
    sprintf(strbuf, "%s %s &", PDF_VIEWER_COMMAND, filep);
    system(strbuf);
#endif

#if defined(NEXTSTEP) || defined(MacOSX)
    sprintf(strbuf, "%s %s", PDF_VIEWER_COMMAND, filep);	/* we use "open" */
    system(strbuf);
#endif

#ifdef MacOS8
    if(firstTime) {
        InitGraf(&qd.thePort);		/* apparently, this is needed */
	firstTime = 0;			/* only once, in case this is a persistent program */
    }
    strcpy(strbuf, filep);
    c2pstr(strbuf);
    err = FSMakeFSSpec(0, 0, (StringPtr)strbuf, &fsSpec);
    if(err != noErr) {
        cpdf_Error(pdf, 1, "ClibPDF", "cannot find PDF file to open: %s", filep);
        return(-2);
    }
    LaunchTheDocument(&fsSpec, false, NULL);
    SetFrontCreator(creator);	/* bring the window to front */
#endif

#if defined(_WIN32) || defined(WIN32)
    #ifndef __BORLANDC__
        ShellExecute(GetActiveWindow(), NULL, filep,
    		NULL, _getcwd(strbuf, 1023), SW_SHOWNORMAL);
    #else
        ShellExecute(GetActiveWindow(), NULL, filep,
    		NULL, getcwd(strbuf, 1023), SW_SHOWNORMAL);
    #endif
#endif

    return(0);
}


int cpdf_launchPreview(CPDFdoc *pdf)
{
    return( cpdf_openPDFfileInViewer(pdf, NULL) );
}



/* ==== MacOS 8 or 7 specific stuff below ============================= */

#ifdef MacOS8

/* SendOpen.c -- MacOS 8 or 7 opening a document programmatically.
   Original source and documentation for this part of code is in
   "Technote 1002: On Launching an App with a Document" (Oct. 95)
   from Apple Developer Technical Support (DTS) at:

    http://developer.apple.com/technotes/tn/tn1002.html

Code snippet illustrating how to launch an application and send it
an open document Apple event.
by John Montbriand, 1995.
Copyright (C) 1995 by Apple Computer, Inc.  All Rights Reserved.

The author (John Montbriand) has compiled and run the following with
the MPW C, MPW PPC, SymantecC, MrC, and CodeWarrior C compilers and
is satisfied that it is bug free for System 7.0 and onward.

*/

/*
#include <Types.h>
#include <QuickDraw.h>
#include <Fonts.h>
#include <Windows.h>
#include <TextEdit.h>
#include <Dialogs.h>
#include <Menus.h>
#include <AppleEvents.h>
#include <StandardFile.h>
#include <Files.h>
#include <Errors.h>
#include <Aliases.h>
*/

/* MapCreatorProcs, when provided as a parameter to LaunchTheDocument,
will be called to allow your app to re-map the creator type for a file
to a different creator type. Your function should either modify *creator
to refer to a preferred application by replacing its value with the
value of the preferred application's creator bytes, or leave *creator
untouched if you want the creator bytes from the document to be used.
*/

/* typedef void (*MapCreatorProc)(FSSpec *document, OSType type, OSType *creator); */

/* LaunchTheDocument attempts to open the document with either an 
application that is already loaded in memory and running, or by finding 
and launching the appropriate application with the document as a 
parameter. If no application can be found with the creator bytes found 
in the document, fnfErr is returned. Local Volumes are searched before 
remote server volumes. Remote server volumes are only searched if 
use_remote_apps is true. If the optional creator remapping procedure is 
provided (remap != NULL), then it will be called to allow the caller to 
re-map the creator bytes to a preferred application for the document's 
type.
*/

OSErr LaunchTheDocument(FSSpec *document, Boolean use_remote_apps, MapCreatorProc remap)
{
    OSErr err;
    AEAddressDesc target_desc; /* address desc for target application */
    AEDescList files_list;
    AEDesc file_desc, parameter_desc;
    AppleEvent the_apple_event;
    AppleEvent the_reply;
    AliasHandle the_alias;
    FInfo fndrInfo;
    OSType app_creator_bytes;
    LaunchParamBlockRec launch_pb;

    /* initialize our records to NULL descriptors */
    AECreateDesc(typeNull, NULL, 0, &target_desc);
    AECreateDesc(typeNull, NULL, 0, &files_list);
    AECreateDesc(typeNull, NULL, 0, &file_desc);
    AECreateDesc(typeNull, NULL, 0, &the_apple_event);
    AECreateDesc(typeNull, NULL, 0, &the_reply);
    AECreateDesc(typeNull, NULL, 0, &parameter_desc);
    the_alias = NULL;

    /* get the target application's creator */
    err = FSpGetFInfo(document, &fndrInfo);
    if (err != noErr) goto launch_the_document_termination;

    /* call the remap function if applicable */
    app_creator_bytes = fndrInfo.fdCreator;
    if (remap != NULL)
	    remap(document, fndrInfo.fdType, &app_creator_bytes);

    /* create an open documents Apple event */
    err = AECreateDesc(typeApplSignature, (Ptr)&app_creator_bytes,
	    sizeof(OSType), &target_desc);
    if (err != noErr) goto launch_the_document_termination;
    err = AECreateAppleEvent(kCoreEventClass, kAEOpenDocuments,
	    &target_desc, kAutoGenerateReturnID, kAnyTransactionID,
		    &the_apple_event);
    if (err != noErr) goto launch_the_document_termination;

    /* create a one element list of files to send in the event */
    err = AECreateList(NULL, 0, false, &files_list);
    if (err != noErr) goto launch_the_document_termination;
    err = NewAlias(NULL, document, &the_alias);
    if (err != noErr) goto launch_the_document_termination;
    HLock((Handle) the_alias);
    err = AECreateDesc(typeAlias, (Ptr) (*the_alias),
	    GetHandleSize((Handle) the_alias), &file_desc);
    HUnlock((Handle) the_alias);
    if (err != noErr) goto launch_the_document_termination;
    err = AEPutDesc(&files_list, 0, &file_desc);
    if (err != noErr) goto launch_the_document_termination; 

    /* add the file list to the open documents event */
    err = AEPutParamDesc(&the_apple_event, keyDirectObject,
	    &files_list);
    if (err != noErr) goto launch_the_document_termination;
    /* send the Apple event */
    err = AESend(&the_apple_event, &the_reply, kAENoReply,
		    kAENormalPriority, kNoTimeOut, NULL, NULL);

    /* if the target could not be found */
    if (err == connectionInvalid) { /* no such app running? */
        HVolumeParam vol_pb;
        DTPBRec desktop_pb;
        FSSpec application;
        GetVolParmsInfoBuffer volinfo;
        HIOParam param_pb;

        /* search desktop files on local volumes first */
        param_pb.ioNamePtr = NULL;
        param_pb.ioBuffer = (Ptr) &volinfo;
        param_pb.ioReqCount = sizeof(volinfo);
        err = fnfErr; /* default return value */
        vol_pb.ioNamePtr = NULL;

        for(vol_pb.ioVolIndex = 1;
		PBHGetVInfoSync((HParmBlkPtr) &vol_pb) == noErr;
		vol_pb.ioVolIndex++) {
	    param_pb.ioVRefNum = vol_pb.ioVRefNum;
	    if (PBHGetVolParmsSync((HParmBlkPtr)&param_pb) == noErr
	    && volinfo.vMServerAdr == 0) {
		desktop_pb.ioCompletion = NULL;
		desktop_pb.ioVRefNum = vol_pb.ioVRefNum;
		desktop_pb.ioNamePtr = NULL;
		desktop_pb.ioIndex = 0;
		if (PBDTGetPath(&desktop_pb) == noErr) {
		    desktop_pb.ioFileCreator = app_creator_bytes;
		    desktop_pb.ioNamePtr = application.name;
		    if (PBDTGetAPPLSync(&desktop_pb) == noErr) {
			application.vRefNum = vol_pb.ioVRefNum;
			application.parID = desktop_pb.ioAPPLParID;
			err = noErr;
			break;
		    }
		}
	    }
        }

	/* if this fails, search remove volumes if allowed... */
	if (err != noErr && use_remote_apps) {
	    for (vol_pb.ioVolIndex = 1;
		    PBHGetVInfoSync((HParmBlkPtr) &vol_pb) == noErr;
		    vol_pb.ioVolIndex++) {
		param_pb.ioVRefNum = vol_pb.ioVRefNum;
		if (PBHGetVolParmsSync((HParmBlkPtr) &param_pb) == noErr
		&& volinfo.vMServerAdr != 0) {
		    desktop_pb.ioCompletion = NULL;
		    desktop_pb.ioVRefNum = vol_pb.ioVRefNum;
		    desktop_pb.ioNamePtr = NULL;
		    desktop_pb.ioIndex = 0;
		    if (PBDTGetPath(&desktop_pb) == noErr) {
			desktop_pb.ioFileCreator = app_creator_bytes;
			desktop_pb.ioNamePtr = application.name;
			if (PBDTGetAPPLSync(&desktop_pb) == noErr) {
			    application.vRefNum = vol_pb.ioVRefNum;
			    application.parID = desktop_pb.ioAPPLParID;
			    err = noErr;
			    break;
			}
		    }
		}
	    }  /* end for(vol_pb.ioVolIndex = 1; ...) */
	}  /* end if(err != noErr && use_remote_apps) */
	if (err != noErr)
	    goto launch_the_document_termination;

	/* coerce the apple event to app parameters */
	err = AECoerceDesc(&the_apple_event, typeAppParameters, &parameter_desc);
	if (err != noErr)
	    goto launch_the_document_termination;

	/* launch the application */
	launch_pb.launchBlockID = extendedBlock;
	launch_pb.launchEPBLength = extendedBlockLen;
	launch_pb.launchFileFlags = 0;
	launch_pb.launchControlFlags = launchContinue + launchNoFileFlags;
	launch_pb.launchAppSpec = &application;
	HLock((Handle) parameter_desc.dataHandle);
	launch_pb.launchAppParameters = (AppParametersPtr) (*parameter_desc.dataHandle);
	err = LaunchApplication(&launch_pb);
	HUnlock((Handle) parameter_desc.dataHandle);
    } /* end if (err == connectionInvalid) */

launch_the_document_termination:
	/* clean up, and go.. */
	if (the_alias != NULL) DisposeHandle((Handle) the_alias);
	AEDisposeDesc(&parameter_desc);
	AEDisposeDesc(&target_desc);
	AEDisposeDesc(&file_desc);
	AEDisposeDesc(&files_list);
	AEDisposeDesc(&the_apple_event);
	AEDisposeDesc(&the_reply);
	return err;
}
/* end of LaunchTheDocument() */

/* This function brings the opened document window to front
   --- based on the following tip from USENET:
    Date: Sat, 17 Oct 1998 21:35:34 -0600
    From: dtp@pluto.njcc.com (David T. Pierson)
    Subject: Re: Apple Event for bringing remote window to front?
    Newsgroups: comp.sys.mac.programmer.help
*/
void SetFrontCreator(OSType inProcessSignature)
{
ProcessSerialNumber  thePSN;
ProcessInfoRec       thePRec;
Boolean              isDone=false;

    thePSN.highLongOfPSN = kNoProcess;
    thePSN.lowLongOfPSN = kNoProcess;
    thePRec.processAppSpec = nil;
    thePRec.processInfoLength = sizeof(ProcessInfoRec);
    thePRec.processName = nil;
    
    while (! isDone) {  
	if(( GetNextProcess(&thePSN)==noErr ) && ( GetProcessInformation(&thePSN, &thePRec)==noErr )) {  
	    if ( thePRec.processSignature== inProcessSignature ) {
		SetFrontProcess(&thePSN);	/* found the app */
		isDone=true;
	    }
	}
	else
	    isDone=true;	/* App not found, but just give up. */
    }
}

/*  end of SetFrontCreator() */


/* sample creator mapping function. The function simply tells
        LaunchTheDocument that every document of type 'TEXT' should be
        opened using MPW Shell
*/

/* --- rest of file commented out ------------------------------------
void MyMapCreator(FSSpec *document, OSType type, OSType *creator) {
        if (type == 'TEXT')
                *creator = 'BOBO';
}


#ifdef powerc
QDGlobals qd;
#endif

void main(void) {
        SFTypeList typeList;
        StandardFileReply reply;

        InitGraf(&qd.thePort);
        InitFonts();
        InitWindows();
        TEInit();
        InitMenus();
        InitDialogs(0);
        InitCursor();
        StandardGetFile(NULL, -1, typeList, &reply);

        if (reply.sfGood)
                LaunchTheDocument(&reply.sfFile, false, MyMapCreator);
}

File SendOpen.r containing Rez source for the 'SIZE' resource.
Because you're working with Apple events, you must
include a SIZE resource with the isHighLevelEventAware bit set: 

resource 'SIZE' (-1, purgeable) {
                reserved,
                acceptSuspendResumeEvents,
                reserved,
                canBackground,
                multiFinderAware,
                backgroundAndForeground,
                dontGetFrontClicks,
                ignoreChildDiedEvents,
                is32BitCompatible,
                isHighLevelEventAware,
                localAndRemoteHLEvents,
                isStationeryAware,
                dontUseTextEditServices,
                reserved,
                reserved,
                reserved,
                524288,
                524288
};

Summary

Launching an application with a document under System 7 requires
the use of the AppleEvent Manager to either directly send an open
document event to an existing application or indirectly send an
open document event to an application as part of its launch sequence.
This Technote describes how you can do this in your applications
and provides a simple example of how you can use the technique.

Further References

"Creating And Sending Apple Events" (Chapter 5), Inside Macintosh: Interapplication Communication 
"Handling the Required Apple Events" (4-11), Inside Macintosh: Interapplication Communication 
"Desktop Manager" (Chapter 9), Inside Macintosh: More Macintosh Toolbox 
"Obtaining Volume Information" (2-145), Inside Macintosh: Files 
"Alias Manager" (Chapter 4), Inside Macintosh: Files 

Acknowledgements

Thanks to Eric Anderson, Nitin Ganatra, and Bob Wambaugh.
------------------------------------------------------------------ */
#endif		/* MacOS8 */

