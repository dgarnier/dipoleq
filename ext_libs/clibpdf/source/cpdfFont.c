/* cpdfFont.c  -- Font management and string width etc.
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf.
 *  or http://www.fastio.com/LICENSE.pdf

## MULTI_THREADING_OK ## 

1999-11-25 [io]  Bad case of open file leak in _cpdf_readFontPathFromMapFile() fixed.
		 Thanks to Pierre-Yves Monnet.
1999-10-14 [io]  Optimization of cpdf_stringWidth().
		 Added cpdf_ascent();
		 Fixed cpdf_capHeight() to return point size, rather than fraction of font size.
1999-09-30 [io]  Correct char width look up mechanism in place.
		 WinAnsiEncoding char width table completed.
1999-09-26 [io]  Added CJK Font Support.
1998-07-11 [io]
*/

/* #define CPDF_DEBUG 	1 */

#include "config.h"
#include "version.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "cpdflib.h"		/* This must be included before all other local include files */
#include "cpdf_private.h"
#include "cglobals.h"
#include "cpdfAFM.h"		/* font width array */

/* ================================================================================== */
/*  Change NUMPSFONTS to reflect the actual number of fonts defined in cpdfAFM.h.     */
/* ================================================================================== */
#ifdef  PS_FONTS39
#define NUMPSFONTS 48		/* including CJK */
#define NUMROMANFONTS 41
#else
#define NUMPSFONTS 14
#define NUMROMANFONTS 14
#endif

#ifdef MacOS8
#include <stat.mac.h>
#include <unix.mac.h>
#include <StandardFile.h>
#include <Memory.h>
#include <Files.h>
#include <Strings.h>	/* c2pstr() */
#endif

/* This won't work for Unicode and all Asian char encodings. */
#define	iskanji(c)	((c)>0x7f)

/* Returns actual captal letter height in points.
   Used mainly for centering string in the height dimension.
*/

/* ## MULTI_THREADING_OK ## */
float cpdf_capHeight(CPDFdoc *pdf)
{
int idxAFM;
float capH;
CPDFfontInfo *fInfo = &(pdf->fontInfos[pdf->currentFont]);
    idxAFM = fInfo->afmIndex;		/* index into AFM data arrays */
    if(fInfo->descLevel < 2 || fInfo->descLevel == 3)	/* Base-14, Metrics embedded or CJK */
        capH = (float)_cpdfCapHeight[idxAFM]/1000.0;		/* get it from static array compiled in */
    else
        capH = (float)(fInfo->fontDesc->capHeight)/1000.0;
    return(capH * pdf->font_size);
}

float cpdf_ascent(CPDFdoc *pdf)
{
int idxAFM;
float ascent;
CPDFfontInfo *fInfo = &(pdf->fontInfos[pdf->currentFont]);
    idxAFM = fInfo->afmIndex;		/* index into AFM data arrays */
    if(fInfo->descLevel < 2 || fInfo->descLevel == 3)	/* Base-14, Metrics embedded or CJK */
        ascent = (float)_cpdfAscender[idxAFM]/1000.0;		/* get it from static array compiled in */
    else
        ascent = (float)(fInfo->fontDesc->ascent)/1000.0;
    return(ascent * pdf->font_size);
}

float cpdf_descent(CPDFdoc *pdf)
{
int idxAFM;
float descent;
CPDFfontInfo *fInfo = &(pdf->fontInfos[pdf->currentFont]);
    idxAFM = fInfo->afmIndex;		/* index into AFM data arrays */
    if(fInfo->descLevel < 2 || fInfo->descLevel == 3)	/* Base-14, Metrics embedded or CJK */
        descent = (float)_cpdfDescender[idxAFM]/1000.0;		/* get it from static array compiled in */
    else
        descent = (float)(fInfo->fontDesc->descent)/1000.0;
    return(descent * pdf->font_size);
}

/* ## MULTI_THREADING_OK ## */
float cpdf_stringWidth(CPDFdoc *pdf, const unsigned char *str)
{
float swid = 0.0;
float fsize = pdf->font_size;
float chspacing = pdf->char_spacing;
float wdspacing = pdf->word_spacing;
unsigned char *p;
int c1, c2, ch;
CPDFfontInfo *fInfo = &(pdf->fontInfos[pdf->currentFont]);

    p = (unsigned char *)str;
    if(fInfo->descLevel == 3) {
	/* Current font is a CJK font */
	unsigned char *endptr;
        endptr = p + strlen((const char *)str);
	while( p < endptr) {
	    c1 = *p++;
	    c2 = (iskanji(c1) ? *p++ :  *p);	/* if kanji 2-nd byte, if ASCII c2 == c1 */
	    if( iskanji(c1) ) {
		ch = (c1 << 8) + c2;
		swid += fsize + chspacing;		/* assume width = 1000 */
	    }
	    else {
		/* roman chars in CJK */
		swid += fsize * 0.5 + chspacing;	/* assume 500 for Roman, HW kana */
		if(c1 == ' ')
		    swid += wdspacing;
	    }
	}
	swid -= chspacing;		/* spacing between chars is one less */
	swid *= pdf->horiz_scaling/100.0;
    }
    else {
	/* Current font is roman or symbol font */
	/* 0=WinAnsiEncoding, 1=MacRomanEncoding, 2=MacExpertEncoding, 3=StandardEncoding */
	int enc = fInfo->encodeIndex; 
	short *encodeLUT, *chWidth;
        if(fInfo->afmIndex == 12 || fInfo->afmIndex == 13)
	    enc = 0;	/* Symbol or ZapfDingbats, just turn off the extra lookup for these */
	if(enc)
	    encodeLUT = cpdf_encode2WinAnsi[enc-1];	/* encoding LUT if not WinAnsi */
	chWidth = fInfo->width;
	fsize *= 0.001;				/* normalize here for 1000 unit grid */
	while((ch = *p++) != '\0') {
	    if(enc)
		ch = encodeLUT[ch];	/* Look up, if encoding is not WinAnsiEncoding */
	    swid += fsize * (float)chWidth[ch] + chspacing;
	    if(ch == ' ')
		swid += wdspacing;
	}
	swid -= chspacing;		/* spacing between chars is one less */
	swid *= pdf->horiz_scaling/100.0;
    }
#ifdef CPDF_DEBUG
    printf("%g inches: %s\n", swid/72.0, str);
#endif
    return(swid);
}


int cpdf_setFont(CPDFdoc *pdf, const char *basefontname, const char *encodename, float size)
{
int i, fontOK= 0, reUse = -1, isNewVal = 0;
int idxAFM = -1;		/* index into AFM metrics/descriptor arrays in cpdfAFM.h */
int foundInPageList = 0;	/* in pageInfos[].fontIdx[] */
int fontIndexFound = 0;		/* in fontInfos[] */
int baseFontFound = 0;
char *realbasefont;
CPDFfontInfo *newFont;
CPDFfontDescriptor *extFontDesc;
CPDFextFont *extFontData, *efd, *extFontReuse;
char fontname[32], *font;

    /* First, look in the static (compiled in) list of internal fonts */
    realbasefont = (char *)basefontname;
    for(i=0; i<NUMPSFONTS; i++) {
        int cmp;
	if(i < NUMROMANFONTS)
	    cmp = strcmp(basefontname, cpdf_fontnamelist[i]);  /* roman/symbol */
	else
	    cmp = strncmp(basefontname, cpdf_fontnamelist[i], strlen(cpdf_fontnamelist[i]));  /* CJK */
	if( cmp == 0) {
	    fontOK = 1;
	    idxAFM = i;		/* index into arrays in cpdfAFM.h */
	    break;
	}
    }

    extFontReuse = NULL;

    /* Look this font name in the list of already loaded external fonts */
    for(i=0; i< pdf->numExtFonts; i++) {
	efd = &pdf->extFontList[i];
	if( strcmp(basefontname, efd->basefont) == 0) {
	    fontOK = 1;
	    reUse = i;
	    extFontReuse = efd;
	    idxAFM = i+NUMPSFONTS;		/* index into arrays in cpdfAFM.h */
#ifdef CPDF_DEBUG
	    fprintf(stderr, "Already loaded: %s\n", basefontname);
#endif
	    break;
	}
    }

    if(!fontOK) {
	/* Font is not in the internal fonts built into ClibPDF or the list of loaded external fonts.
	   Check font directories for PFM/PFB files.
	*/
	newFont = &pdf->fontInfos[pdf->numFonts];
	extFontDesc = (CPDFfontDescriptor *)malloc(sizeof(CPDFfontDescriptor));
	_cpdf_malloc_check((void *)extFontDesc);
	extFontData = &pdf->extFontList[pdf->numExtFonts];

	/* The function below will look for font in the following places in the order listed:
	   1. Current directory.
	   2. Directries specified by cpdf_setFontDirectories();
	   3. Font files specified by a fontmap file set by cpdf_setFontMapFile();
		Note: For cases 1 and 2, the PFM and PFB files must bear exactly the
			same base name (including case) as font name specified in
			cpdf_setFont();  i.e.,  if you use: cpdf_setFont("Tekton", ..);
			Then, files "Tekton.pfm" and "Tekton.pfb" must exist.
		      For case 3, font filenames may differ from that specified in
			cpdf_setFont(), becuase the mapfile provides the mapping information.
	*/
	if( _cpdf_loadExternalFont(pdf, basefontname, newFont, extFontDesc, extFontData) == 0) {
	    /* External font successfuly loaded */
	    extFontReuse = extFontData;		/* we need it here too */
	    extFontData->parentDesc = extFontDesc;
	    if(pdf->numExtFonts == 0)
	        extFontData->prev = NULL;	/* set "prev" to point to the previous item in the stack */
	    else
		extFontData->prev = &pdf->extFontList[pdf->numExtFonts-1];
	    extFontData->basefont = (char *)malloc(strlen(basefontname)+1);
	    _cpdf_malloc_check((void *)extFontData->basefont);
	    strcpy(extFontData->basefont, basefontname);	/* This is stub name as used by cpdf_setFont() */
	    extFontData->doneAlready = 0;
	    extFontDesc->doneAlready = 0;
	    idxAFM = pdf->numExtFonts + NUMPSFONTS;
	    pdf->numExtFonts++;
	}
	else {
	    /* Font not found anywhere in the 3 places, so substitute something from base-14. */
	    cpdf_Error(pdf, 1, "ClibPDF", "Font not found: %s >> substituting Times-Roman.", basefontname);
	    realbasefont = "Times-Roman";
	    idxAFM = 4;		/* index of Times-Roman */
	    free(extFontDesc);
	    /* free(extFontData); */
	}
    }

    /* check if this font has already been included or new font */
    if( (isNewVal= _isNewFont(pdf, realbasefont, encodename, &fontIndexFound, &baseFontFound)) ) {
	/* This section is entered if either realbasefont or encode name or both differ from
	   existing font in the list
	*/
#ifdef CPDF_DEBUG
	fprintf(stderr, "%s (%s): FI=%d, baseFI=%d\n", realbasefont, encodename, fontIndexFound, baseFontFound);
#endif
	if(pdf->numFonts >= pdf->nMaxFonts) {
	    cpdf_Error(pdf, 1, "ClibPDF", "Too many fonts in this PDF: %d. Increase limit by cpdf_open()",
			pdf->numFonts);
	    return(1);
	}
        sprintf(fontname, "Fcpdf%d", pdf->numFonts);
	font = fontname;
	newFont = &pdf->fontInfos[pdf->numFonts];
	if(idxAFM >= NUMPSFONTS ) {
 	    if(extFontReuse != NULL && baseFontFound >= 0) {
		/* copy from existing instance of the same base font, but with different encoding */
		memcpy(newFont, &pdf->fontInfos[baseFontFound], sizeof(CPDFfontInfo));
	    }
	}
	pdf->currentFont = pdf->numFonts;			/* index of curren font */
	newFont->afmIndex = idxAFM;		/* index into the AFM array in cpdfAFM.h */
	if(idxAFM < 14)
	    newFont->descLevel = 0;		/* base 14 font */
	else if(idxAFM < NUMROMANFONTS)
	    newFont->descLevel = 1;		/* metrics/descriptors from cpdfAFM.h */
	else if(idxAFM < NUMPSFONTS) {
	    newFont->descLevel = 3;		/* CJK font */
	    cpdf_setPDFLevel(pdf, 1, 2);	/* CJK requires PDF-1.2 */
	}
	else
	    newFont->descLevel = 2;		/* external font */

	/* This name is like "Fcpdf1" that appears in content */
	newFont->name = (char *)malloc((size_t)(strlen(fontname) + 1));
	_cpdf_malloc_check((void *)newFont->name);
	strcpy(newFont->name, fontname);
	/* name passed in cpdf_setFont() */
	newFont->stubFont = (char *)malloc((size_t)(strlen(realbasefont) + 1));
	_cpdf_malloc_check((void *)newFont->stubFont);
	strcpy(newFont->stubFont, realbasefont);

	/* Only for internal fonts.  Embedded fonts use real name from PFM file */
	if(newFont->descLevel < 2 || newFont->descLevel == 3) {
	    newFont->baseFont = (char *)malloc((size_t)(strlen(realbasefont) + 1));
	    _cpdf_malloc_check((void *)newFont->baseFont);
	    strcpy(newFont->baseFont, realbasefont);
	}
	else {
	    newFont->baseFont = (char *)malloc((size_t)(strlen(newFont->fontDesc->fontName) + 1));
	    _cpdf_malloc_check((void *)newFont->baseFont);
	    strcpy(newFont->baseFont, newFont->fontDesc->fontName);
	}

	if(encodename != NULL) {
	    newFont->encoding = (char *)malloc((size_t)(strlen(encodename) + 1));
	    _cpdf_malloc_check((void *)newFont->encoding);
	    strcpy(newFont->encoding, encodename);
	    newFont->encodeIndex = 0;		/* default to WinAnsiEncoding */
	    for(i=0; i<4; i++) {
		if(strcmp(encodename, cpdf_fontEncodings[i]) == 0) {
		    newFont->encodeIndex = i;
		    break;
		}
	    }
	}
	else {
	    newFont->encoding = NULL;
	    newFont->encodeIndex = 3;	/* StandardEncoding */
	}

	/* --- Extra font metrics/descriptor info added starting with v1.10 --- */
	if(newFont->descLevel < 2) {
	    newFont->firstchar = 32;
	    newFont->lastchar = 255;
	    newFont->width = charWidth[idxAFM];	/* in cpdfAFM.h */
	}

	if(newFont->descLevel == 1 || newFont->descLevel == 3) {
	    /* allocate memory for Font Descriptor for non-base 14 and CJK fonts */
	    newFont->fontDesc = (CPDFfontDescriptor *)malloc((size_t)sizeof(CPDFfontDescriptor));
	    _cpdf_malloc_check((void *)newFont->fontDesc);
	    newFont->fontDesc->ascent = _cpdfAscender[idxAFM];
	    newFont->fontDesc->capHeight = _cpdfCapHeight[idxAFM];
	    newFont->fontDesc->descent = _cpdfDescender[idxAFM];
	    newFont->fontDesc->flags = _cpdfFontFlags[idxAFM];
	    newFont->fontDesc->italicAngle = _cpdfItalicAngle[idxAFM];
	    newFont->fontDesc->stemV = _cpdfStemV[idxAFM];
	    newFont->fontDesc->xHeight = _cpdfXHeight[idxAFM];

	    /* Next 5 items are valid only for CJK fonts */
	    newFont->fontDesc->stemH = _cpdf_stemH[idxAFM];
	    newFont->fontDesc->missingWidth = _cpdf_missingWidth[idxAFM];
	    newFont->fontDesc->leading = _cpdf_leading[idxAFM];
	    newFont->fontDesc->maxWidth = _cpdf_maxWidth[idxAFM];
	    newFont->fontDesc->avgWidth = _cpdf_avgWidth[idxAFM];
	    newFont->fontDesc->style = cpdf_fontStyleList[idxAFM];
	    /* for CJK font only */
	    if(newFont->descLevel == 3) {
	        newFont->descendFont = (CPDFdescendantFont *)malloc((size_t)sizeof(CPDFdescendantFont));
	        _cpdf_malloc_check((void *)newFont->descendFont);
		newFont->descendFont->doneAlready = 0;
		newFont->descendFont->DW = 1000;
		newFont->descendFont->cidSysInfo = cpdf_cidSysInfo[idxAFM];
		newFont->descendFont->cidFontWidth = cpdf_cidFontWidth[idxAFM];
	    }
	    else
		newFont->descendFont = NULL;

	    newFont->fontDesc->doneAlready = 0;
	    for(i=0; i<4; i++)
		newFont->fontDesc->fontBBox[i] = _cpdfFontBBox[idxAFM][i];
	    /* Point fontName in descriptor to the baseFont name in the Font object above.
	       This is OK for PS Type1 fonts.  It will have to be modified once
	       other font types are supported.
	    */
	    newFont->fontDesc->fontName = newFont->baseFont;
	    newFont->fontDesc->extFont = NULL;
	}
	else if(newFont->descLevel == 2) {
	    /* newFont->fontDesc = extFontReuse->parentDesc;
	       newFont->fontDesc->extFont = extFontReuse;
	    */
	}
	else {
	    newFont->fontDesc = NULL;		/* No font descriptor is needed for base-14 fonts */
	    newFont->descendFont = NULL;	/* No descendant font is needed for base-14 fonts */
	}

	pdf->numFonts++;	/* increment font count for next font */
    }
    else {
	/* This font+encoding has been used already, and is in fontInfos[fontIndexFound] */
	pdf->currentFont = fontIndexFound;		/* index of curren font */
	font = pdf->fontInfos[fontIndexFound].name;
    }

    /*  Store the index of this font in the used font list in pageInfos[].
	Multiple pages can share fonts.
    */
    {   /* braces for local variable */
	CPDFintList *iLprev, *iL = pdf->pageInfos[pdf->currentPage].fontIdx;
	iLprev = iL;
	foundInPageList = 0;
	while(iL != NULL) {
	    if(iL->value == pdf->currentFont) {
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
	    iLprev->value = pdf->currentFont;	/* store the index to current font */
	    iL->next = NULL;
	    iL->value = -1;
	    pdf->pageInfos[pdf->currentPage].npFont++;
	}
    }
    /* Write font spec line to Contents stream */
    pdf->inTextObj = 1;
    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "/%s %.3f Tf\n", font, size);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
        fprintf(pdf->fpcontent, "/%s %.3f Tf\n", font, size);
    pdf->font_size = size;
    
    return(0);
}

/* Returns 0 if font is loaded successfully. If font not found or error, returns non-zero */
int _cpdf_loadExternalFont(CPDFdoc *pdf, const char *fontname,
	CPDFfontInfo *eFI, CPDFfontDescriptor *eFDesc, CPDFextFont *eFData)
{
int retval = 0;
int fontlocation = 0;	/* 0=current dir, 1=in specified font dir, 2=fontmap file */
char filePFM[1024], filePFB[1024];
volatile int fnsize = sizeof(filePFM);

    /* int _cpdf_readPFM(char *pfmfile, CPDFfontInfo *fI, CPDFfontDescriptor *fD); */
    strcpy(filePFM, fontname);
    strcat(filePFM, ".pfm");
    if((retval = _cpdf_readPFM(filePFM, eFI, eFDesc))) {	/* trying in current dir */
	/* not in current dir */
	snprintf(filePFM, fnsize, "%s%c%s.pfm", pdf->pfm_directory, DIR_SEPARATOR, fontname);	/* prepend possibly set dir */
	if((retval = _cpdf_readPFM(filePFM, eFI, eFDesc))) {
	    /* Not in specified directory path. So try finding it in Fontmap file. */
	    if((retval =_cpdf_readFontPathFromMapFile(pdf, fontname, filePFM, filePFB))) {
		if(retval == -1)
		    cpdf_Error(pdf, 1, "ClibPDF", "_cpdf_loadExternalFont() Font not found: %s", fontname);
		else if(retval == 1)
		    cpdf_Error(pdf, 1, "ClibPDF", "_cpdf_loadExternalFont() Fontmap file not found: %s",
				pdf->fontmapfile);
		else if(retval == 2)
		    cpdf_Error(pdf, 1, "ClibPDF", "_cpdf_loadExternalFont() Syntax error in fontmap file: %s",
				pdf->fontmapfile);
		return(retval);
	    }
	    /* There was an entry in map file, so try with that */
	    if((retval = _cpdf_readPFM(filePFM, eFI, eFDesc))) {
	        /* Could not find font anywhere. */
		cpdf_Error(pdf, 1, "ClibPDF",
			"_cpdf_loadExternalFont() Fontmap file %s or PFM file %s for %s is invalid.",
			pdf->fontmapfile, filePFM, fontname);
		return(retval);
	    }
	    else
	        fontlocation = 2;
	}
	else {
	  /* PFM was in specified dir */
	  fontlocation = 1;
	}
    }
    /* Now, we are going to read the PFB file.  We assume PFB font will be found following
	the method by which we found the PFM file.
    */
    /* int _cpdf_readPFB(const char *pfbfile, CPDFextFont *eFI); */
    if(fontlocation == 0) {
	strcpy(filePFB, fontname);
	strcat(filePFB, ".pfb");
	if((retval = _cpdf_readPFB(filePFB,  eFData)))
	    return(retval);
    }
    else if(fontlocation == 1) {
	snprintf(filePFB, fnsize, "%s%c%s.pfb", pdf->pfb_directory, DIR_SEPARATOR, fontname);
	if((retval = _cpdf_readPFB(filePFB, eFData)))
	    return(retval);
    }
    else { /* fontlocation == 2 */
	/* from fontmap file.  filePFB should have been read when looking for PFM file above. */
	if((retval = _cpdf_readPFB(filePFB,  eFData)))
	    return(retval);
    }
    /* All OK */
    eFI->fontDesc = eFDesc;
    eFDesc->extFont = eFData;
    return(0);
}

/* if both base and encodename are different returns 2 */
/* if only encode name is different return 1 */
/* if both are the same, returns 0 */
int _isNewFont(CPDFdoc *pdf, const char *basefontname, const char *encodename, int *fontFound, int *baseFound)
{
int i, isNew = 2;
    *baseFound = -1;
    /* look in the list */
    for(i=0; i < pdf->numFonts; i++) {
	if( (strcmp(basefontname, pdf->fontInfos[i].stubFont) == 0) ) {
	    if((encodename == NULL && pdf->fontInfos[i].encoding == NULL) ||
		(encodename != NULL && pdf->fontInfos[i].encoding != NULL
		   && strcmp(encodename, pdf->fontInfos[i].encoding) == 0) ) {
		isNew = 0;		/* already defined: basefont AND encoding are same */
		*fontFound = i;	/* return found font index */
		*baseFound = i;
		break;
	    }
	}
    }
    if(isNew) {
	for(i=0; i < pdf->numFonts; i++) {
	    if( (strcmp(basefontname, pdf->fontInfos[i].stubFont) == 0) ) {
		    isNew = 1;		/* already defined: basefont name only (encoding differs) */
		    *baseFound = i;
		    break;
	    }
	}
	*fontFound = pdf->numFonts;	/* index of the new font */
    }

    return(isNew);
}

/*  Font map file format:
# comment line
% comment line
Font-Name-Italic	(/usr/local/fonts/pfm/FNI___.pfm)	(/usr/local/fonts/pfb/FNI___.pfb)
 ...
*/

#define	MAX_FDIR_PREFIX	10

int _cpdf_readFontPathFromMapFile(CPDFdoc *pdf, const char *fontname, char *filePFM, char *filePFB)
{
FILE *fp;
int i, foundFlag = -1, wasAbsPath = 0;
int foundFontDirSpec = 0;
char lbuf[2100], keybuf[128], pfmdir[1024], pfbdir[1024], temp_buf[1024];
int Nabspref = 0, len=0;
char *p, abspref[MAX_FDIR_PREFIX][8];
char dir_sep[2];

    if((fp = fopen(pdf->fontmapfile, "r")) == NULL)
	return(1);

    /* Some initializations */
    dir_sep[0] = DIR_SEPARATOR;		/* direcotry separator char (from config.h) */
    dir_sep[1] = 0;
    strcpy(pfmdir, "");
    strcpy(pfbdir, "");

    /* while #1 */
    while( ( fgets(lbuf, 2098, fp) ) != NULL ) {
	/* Skip comment or blank lines */
	if(lbuf[0] == '#' || lbuf[0] == '%' || strlen(lbuf) < 3)
	    continue;		/* pretty much blank line */

	/* Handle special ones: $Font_Directories and $Abs_Path_Prefix */
	if(lbuf[0] == '$') {
	    if(strncmp(lbuf, "$Font_Directories", 17 ) == 0) {
		/* Read font directories to be prepended to each font map entries */
		if((p = strchr(lbuf, '(')) == NULL) {
		    fclose(fp);
		    return(2);
		}
		p++;	/* char after first ( */
		sscanf(p, "%[^)]", pfmdir);
		if((p = strchr(p, '(')) == NULL) {
		    fclose(fp);
		    return(2);
		}
		p++;	/* char after 2-nd ( */
		sscanf(p, "%[^)]", pfbdir);
#ifdef CPDF_DEBUG
			fprintf(stderr, "PFM dir path in fontmap= %s\n", pfmdir);
			fprintf(stderr, "PFB dir path in fontmap= %s\n", pfbdir);
#endif
		foundFontDirSpec = 1;
	    }
	    else if(strncmp(lbuf, "$Abs_Path_Prefix", 16) == 0) {
		/* Now, we get absolut path prexif specs. */
		p = lbuf;
		/* while #2 */
		while((p = strchr(p, '(')) != NULL && Nabspref < MAX_FDIR_PREFIX) {
		    p++;	/* char after first ( */
		    strcpy(keybuf, "");
		    sscanf(p, "%[^)]", keybuf);
		    len = strlen(keybuf);
		    if(len > 0 && len < 8) {
		        strcpy(abspref[Nabspref], keybuf);
#ifdef CPDF_DEBUG
			fprintf(stderr, "Abs_PAthPrefix[%d]= %s\n", Nabspref, abspref[Nabspref]);
#endif
			Nabspref++;
		    }
		} /* end while #2 ((p= strchr(p, '(')) .... */
	    }
	    continue;		/* while loop */
	}  /* end if(lbuf[0] == '#' ...) */

	sscanf( lbuf, "%s", keybuf);	/* first word */
	if(strcmp(fontname, keybuf) == 0) {
	    /* Now we have what looks like a fontname map spec */
	    /* First get filePFM */
	    if((p = strchr(lbuf, '(')) == NULL) {
		fclose(fp);
		return(2);
	    }
	    p++;	/* char after first ( */
	    sscanf(p, "%[^)]", temp_buf);	/* PFM file/path filePFM */
	    /* OK, prepend pfmdir unless spec is an absolute path */
	    wasAbsPath = 0;
	    if(Nabspref && foundFontDirSpec) {
#ifdef CPDF_DEBUG
			fprintf(stderr, "Checking for absolute path\n");
#endif
		for(i=0; i<Nabspref; i++) {
		    len = strlen(abspref[i]);
		    if(strncmp(temp_buf, abspref[i], len) == 0) {
			/* absolute path, so return spec in () as is */
		        strcpy(filePFM, temp_buf);
			wasAbsPath = 1;
#ifdef CPDF_DEBUG
			fprintf(stderr, "Absolute path: %s, %d\n", temp_buf, len);
#endif
			break;	/* for loop */
		    }
		} /* end for(i=..) */
		if( !wasAbsPath) {
		    /* not absolute path, so prepend pfmdir */
#ifdef CPDF_DEBUG
		    fprintf(stderr, "Not absolute path, so prepending %s%s\n", pfmdir, dir_sep);
#endif
		    strcpy(filePFM, pfmdir);
		    strcat(filePFM, dir_sep);	/* add dire sep */
		    strcat(filePFM, temp_buf);
		}
	    }
	    else { /* No $* entries found in fontmap file */
		strcat(filePFM, temp_buf);
	    }

	    /* Next, get filePFB */
	    if((p = strchr(p, '(')) == NULL) {
		fclose(fp);
		return(2);
	    }
	    p++;	/* char after 2-nd ( */
	    sscanf(p, "%[^)]", temp_buf);	/* PFB file/path filePFB */
	    /* OK, prepend pfbdir unless spec is an absolute path */
	    wasAbsPath = 0;
	    if(Nabspref && foundFontDirSpec) {
		for(i=0; i<Nabspref; i++) {
		    len = strlen(abspref[i]);
		    if(strncmp(temp_buf, abspref[i], len) == 0) {
			/* absolute path, so return spec in () as is */
		        strcpy(filePFB, temp_buf);
			wasAbsPath = 1;
			break;	/* for loop */
		    }
	        } /* end for(i=..) */
		if( !wasAbsPath) {
		    /* not absolute path, so prepend pfbdir */
		    strcpy(filePFB, pfbdir);
		    strcat(filePFB, dir_sep);	/* add dire sep */
		    strcat(filePFB, temp_buf);
		}
	    }
	    else { /* No $* entries found in fontmap file */
		strcat(filePFB, temp_buf);
	    }

#ifdef CPDF_DEBUG
	    fprintf(stderr, "PFM file from fontmap= %s\n", filePFM);
	    fprintf(stderr, "PFB file from fontmap= %s\n", filePFB);
#endif
	    foundFlag = 0;		/* we found the entry */
	    break;			/* break while() loop */
	} /* end if(strcmp(fontname, keybuf) == 0)  .. */
    } /* end while #1 */
    fclose(fp);
    return(foundFlag);
}

void cpdf_setFontMapFile(CPDFdoc *pdf, const char *mapfile)
{
	strncpy(pdf->fontmapfile, mapfile, 1023);
}

void cpdf_setFontDirectories(CPDFdoc *pdf, const char *pfmdir, const char *pfbdir)
{
	strncpy(pdf->pfm_directory, pfmdir, 1020);
	strncpy(pdf->pfb_directory, pfbdir, 1020);
}

int  _cpdf_freeAllFontInfos(CPDFdoc *pdf)
{
CPDFfontInfo *tfont;
CPDFextFont *efd;
int i;

    for(i=0; i< pdf->numExtFonts; i++) {
	efd = &pdf->extFontList[i];
	if(efd->basefont) {
	    free(efd->basefont);
	    efd->basefont = NULL;
	}
	if(efd->fontMemStream) {
	    cpdf_closeMemoryStream(efd->fontMemStream);
	}
    }

    for(i=0; i< pdf->numFonts; i++) {
#ifdef CPDF_DEBUG
	fprintf(stderr, "Freeing font: %d\n", i);
#endif
	tfont = &pdf->fontInfos[i];
	if(tfont->name) {
	    free(tfont->name);
	    tfont->name = NULL;
	}
	if(tfont->baseFont) {
	    free(tfont->baseFont);
	    tfont->baseFont = NULL;
	}
	if(tfont->stubFont) {
	    free(tfont->stubFont);
	    tfont->stubFont = NULL;
	}
	if(tfont->encoding) {
	    free(tfont->encoding);
	    tfont->encoding = NULL;
	}

	/* Free descendant font space for CJK fonts */
	if(tfont->descLevel == 3 && tfont->descendFont) {
	    free(tfont->descendFont);
	    tfont->descendFont = NULL;
	}

	if(tfont->descLevel && tfont->fontDesc) {
	    /* Free font name string memory for external font */
	    if(tfont->descLevel == 2 && tfont->fontDesc->fontName) {
		free(tfont->fontDesc->fontName);
	    }
	    free(tfont->fontDesc);
	    tfont->fontDesc = NULL;
	}
	if(tfont->descLevel == 2 && tfont->width) {
	    free(tfont->width);
	    tfont->width = NULL;
	}
    }
    return(0);
}


