/* cpdf_private.h -- C language PRIVATE API definitions for ClibPDF library
 * ---- Private Functions and Macros: DO NOT CALL THESE FUNCTIONS.
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf
------------------------------------------------------------------------------------
1999-08-17 [IO]
	Private functions separated from cpdflib.h.
	## MULTI_THREADING_OK ##
*/

#ifndef __CLIBPDF_PRIV_H__
#define __CLIBPDF_PRIV_H__

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

/* =============================================================================== */

void _cpdf_defaultErrorHandler(CPDFdoc *pdf, int level, const char* module, const char* fmt, va_list ap);
void _cpdf_defaultGlobalErrorHandler(int level, const char* module, const char* fmt, va_list ap);
int _Big_Endian(void);
int _cpdf_inc_docID(CPDFdoc *pdf);
unsigned char _cpdf_nibbleValue(unsigned char hexchar);

void _cpdf_pdfWrite(CPDFdoc *pdf, char *s);
void _cpdf_pdfBinaryWrite(CPDFdoc *pdf, char *s, long length);
void _cpdf_initDocumentStruct(CPDFdoc *pdf);
void _cpdf_finalizeLinearAll(CPDFdoc *pdf);

long _cpdf_WriteCatalogObject(CPDFdoc *pdf, int objNumber);
long _cpdf_WriteOutlinesObject(CPDFdoc *pdf, int objNumber);
long _cpdf_WritePagesObject(CPDFdoc *pdf, int objNumber);
long _cpdf_WritePageObject(CPDFdoc *pdf, CPDFpageInfo *pInf);
long _cpdf_WriteContentsFromFile(CPDFdoc *pdf, CPDFpageInfo *pInf);
long _cpdf_WriteContentsFromMemory(CPDFdoc *pdf, CPDFpageInfo *pInf);
long _cpdf_WriteProcSetArray(CPDFdoc *pdf, int objNumber);
long _cpdf_WriteFont(CPDFdoc *pdf, CPDFfontInfo *fI);
long _cpdf_WriteFontDescriptor(CPDFdoc *pdf, CPDFfontInfo *fInf);
long _cpdf_WriteDescendantFont(CPDFdoc *pdf, CPDFfontInfo *fInf);
long _cpdf_WriteFontData(CPDFdoc *pdf, CPDFextFont *eFD);
long _cpdf_WriteImage(CPDFdoc *pdf, CPDFimageInfo *imgInf);
long _cpdf_WriteAnnotation(CPDFdoc *pdf, CPDFannotInfo *aInf);
long _cpdf_WriteProducerDate(CPDFdoc *pdf, int objNumber);
long _cpdf_WriteXrefTrailer(CPDFdoc *pdf, int objNumber);

/* Page related private functions */
int  _cpdf_closeContentMemStreamForPage(CPDFdoc *pdf, int page);
int  _cpdf_freeAllPageInfos(CPDFdoc *pdf);

/* Annotation and hyper link private functions */
int  _cpdf_freeAllAnnotInfos(CPDFdoc *pdf);

/* Text related private functions */
void _cpdf_resetTextCTM(CPDFdoc *pdf);

/* Font related private functions */
int _cpdf_loadExternalFont(CPDFdoc *pdf, const char *fontname,
	CPDFfontInfo *eFI, CPDFfontDescriptor *eFDesc, CPDFextFont *eFData);
int _cpdf_readPFB(const char *pfbfile, CPDFextFont *eFI);
int _cpdf_readPFM(const char *pfmfile, CPDFfontInfo *fI, CPDFfontDescriptor *fD);
int _cpdf_readFontPathFromMapFile(CPDFdoc *pdf, const char *fontname, char *filePFM, char *filePFB);
int  _cpdf_freeAllFontInfos(CPDFdoc *pdf);
int _isNewFont(CPDFdoc *pdf, const char *basefontname, const char *encodename, int *fontFound, int *baseFound);


/* Image related private functions */
int _isNewImage(CPDFdoc *pdf, const char *filepath, int imgSel, int *imageFound);
int  _cpdf_freeAllImageInfos(CPDFdoc *pdf);
int _cpdf_copyPDFimageHeader(CPDFdoc *pdf, CPDFimageInfo *timgInfo);


void str_append_int(char *buf, int num);
int _cpdf_file_open(CPDFdoc *pdf);
void _cpdf_file_close(CPDFdoc *pdf);

/* Domain related (private) */
void _do_meshLines_X(CPDFplotDomain *aDomain);
void _do_meshLines_Y(CPDFplotDomain *aDomain);

/* Axis related (private) */
void _do_oneTick(CPDFdoc *pdf, CPDFaxis *anAx, float vt, float ticlen);
void _do_linearTics(CPDFdoc *pdf, CPDFaxis *anAx);
void _do_logTics(CPDFdoc *pdf, CPDFaxis *anAx);
void _do_timeTics(CPDFdoc *pdf, CPDFaxis *anAx);
char *fix_trailingZeros(char *sstr);
void _do_oneNumber(CPDFdoc *pdf, CPDFaxis *anAx, float v, float ticlen);
void _do_oneTimeNumber(CPDFdoc *pdf, CPDFaxis *anAx, float v, struct tm *vtm, int majorBumpVar, float ticlen);

void _do_linearNumbers(CPDFdoc *pdf, CPDFaxis *anAx);
void _do_logNumbers(CPDFdoc *pdf, CPDFaxis *anAx);
void _do_timeNumbers(CPDFdoc *pdf, CPDFaxis *anAx);
void _do_axisLabel(CPDFdoc *pdf, CPDFaxis *anAx);
int _bittest(int aNumber, int bitpos);
void _setDefaultTimeBumpVar(float fndays, int *minorBumpVar, int *majorBumpVar, int *minorBump, int *majorBump);
float _bump_tm_Time(struct tm *rT, struct tm *vT, int bumpVar, int bump);
void _printfTime(struct tm *vtm);
char *_yearFormat(int year, int flag, char *stbuf);

/* memory stream debug function */

void _checkMemMagic(char *idstr, CPDFmemStream *memStream);
void _cpdf_malloc_check(void *buf);

/* Private outline (book mark) functions */
char *_cpdf_copy_string_with_malloc(const char *action_dict);
char *_cpdf_dest_attribute(int mode, float p1, float p2, float p3, float p4);
void _cpdf_serializeOutlineEntries(int *objCount, CPDFoutlineEntry *first, CPDFoutlineEntry *last);
void _cpdf_WriteOutlineEntries(CPDFdoc *pdf, CPDFoutlineEntry *first, CPDFoutlineEntry *last);
long _cpdf_WriteOneOutlineEntry(CPDFdoc *pdf, CPDFoutlineEntry *olent);
void _cpdf_freeAllOutlineEntries(CPDFoutlineEntry *first, CPDFoutlineEntry *last);


/* Private functions for arc drawing */
void _cpdf_arc_small(CPDFdoc *pdf, float x, float y, float r, float midtheta, float hangle, int mvlnto0, int ccwcw);


#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif  /*  __CLIBPDF_PRIV_H__  */

