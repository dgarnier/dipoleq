/* cpdflib.h -- C language API definitions for ClibPDF library
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf
------------------------------------------------------------------------------------
*/


#ifndef __CLIBPDF_H__
#define __CLIBPDF_H__

#include <stdio.h>
#include <time.h>
#include <stdarg.h>

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef __CYGWIN__
#undef _WIN32
#endif

#ifdef _WIN32
/* 4305 is loss of precision, e.g., by assigning double to float.
   but I am tired of VC++ generating errors on OK things like:
   float a = 1.23;
   Somehow, VC++ thinks 1.23 is 'const double'
   So, here is a pragma to turn off that waring until I find a
   switch that makes it treat 1.23 as float as well.
*/
#pragma warning (disable: 4305 4244)
/* #pragma warning (disable: 4244) */

#endif

/* Change these parameters if you are using ClibPDF to create REALLY large PDF files
   all the time.  These are default values, but may be changed by calling
   cpdf_setGlobalDocumentLimits() before cpdf_open().
   These 5 items may be set PER PDF file.
*/
#define NMAXOBJECTS	5000		/* maximum number of objects for XREF */
#define NMAXPAGES	100		/* maximum number of pages */
#define NMAXFONTS	100		/* maximum number of fonts */
#define NMAXIMAGES	100		/* maximum number of different images */
#define NMAXANNOTS	100		/* maximum number of annotations and links */

/* These three may be set on a per page basis via cpdf_setPerPageLimits() */
#define NMAXFONTS_PP	50		/* maximum number of fonts per page */
#define NMAXIMAGES_PP	50		/* maximum number of different images per page */
#define NMAXANNOTS_PP	50		/* maximum number of annotations and links per page */

/* --------------------------------------------------------------------------------- */
#define PI		3.141592654
#ifndef YES
  #define YES		1
#endif
#ifndef NO
  #define NO		0
#endif

/* for use with cpdf_createTimeAxis(), and  cpdf_createTimeAxis() */
#define LINEAR		0
#define LOGARITHMIC	1
#define TIME		2

/* for use with cpdf_setTimeAxisNumberFormat() */
#define MONTH_NUMBER	0
#define MONTH_NAME	1
#define YEAR_FULL	0
#define YEAR_4DIGIT	0
#define YEAR_2DIGIT	1

/* mesh interval and offset parameters for cpdf_setLinearMeshParams() */
#define X_MESH		0
#define Y_MESH		1

/* for cpdf_pageInit() */
#define PORTRAIT	0
#define LANDSCAPE	1
#define ALTLANDSCAPE	2

#define inch		72.0
#define cm		28.3464566929
#define POINTSPERINCH	72.0   		/* number of points per inch */
#define POINTSPERCM	28.3464566929	/* number of points per cm */
#define POINTSPERMM	2.83464566929	/* number of points per mm */

/* Conv factor to get char height of char '0' from font size.
   0.676 for Times-Roman, 0.688 for Times-Bold, 0.703 for Helvetica, 0.71 for Helvetica-Bold.
   This does vary from one char to another.  The number below is a compromize value.
*/
#define FONTSIZE2HEIGHT		0.7

/* standard page sizes in points */
#define LETTER			"0 0 612 792"
#define LEGAL			"0 0 612 1008"
#define A4			"0 0 595 842"
#define B5			"0 0 499 708"
#define C5			"0 0 459 649"
#define DL			"0 0 312 624"
#define EXECUTIVE		"0 0 522 756"
#define COMM10			"0 0 297 684"
#define MONARCH			"0 0 279 540"
#define FILM35MM		"0 0 528 792"

#define DEFAULT_PAGESIZE	LETTER

/* Log axis tick/number selector masks */
#define LOGAXSEL_1		0x0002
#define LOGAXSEL_13		0x000A
#define LOGAXSEL_125		0x0026
#define LOGAXSEL_12468		0x0156
#define LOGAXSEL_12357		0x00AE
#define LOGAXSEL_123456789	0x03FE
#define LOGAXSEL_MIN		0x0001
#define LOGAXSEL_MAX		0x0400

/* Text Rendering Mode: cpdf_setTextRenderingMode() */
#define	TEXT_FILL		0
#define	TEXT_STROKE		1
#define	TEXT_FILL_STROKE	2
#define	TEXT_INVISIBLE		3
#define	TEXT_FILL_CLIP		4
#define	TEXT_STROKE_CLIP	5
#define	TEXT_FILL_STROKE_CLIP	6
#define	TEXT_CLIP		7

/* Text centering mode: cpdf_rawTextAligned(), cpdf_textAligned() */
#define	TEXTPOS_LL	0	/* lower left */
#define	TEXTPOS_LM	1	/* lower middle */
#define	TEXTPOS_LR	2	/* lower right */
#define	TEXTPOS_ML	3	/* middle left */
#define	TEXTPOS_MM	4	/* middle middle */
#define	TEXTPOS_MR	5	/* middle right */
#define	TEXTPOS_UL	6	/* upper left */
#define	TEXTPOS_UM	7	/* upper middle */
#define	TEXTPOS_UR	8	/* upper right */

/* For in-line image placement function: cpdf_placeInLineImage() */
#define IMAGE_MASK	0
#define	CS_GRAY		1
#define	CS_RGB		2
#define	CS_CMYK		3

/* For cpdf_pointer() */
#define PTR_RIGHT	0
#define PTR_DOWN	1
#define PTR_LEFT	2
#define PTR_UP		3

/* Page transition types: cpdf_setPageTransition() */
#define	TRANS_NONE	0
#define	TRANS_SPLIT	1
#define	TRANS_BLINDS	2
#define	TRANS_BOX	3
#define	TRANS_WIPE	4
#define	TRANS_DISSOLVE	5
#define	TRANS_GLITTER	6


typedef enum {
	SECOND,
	MINUTE,
	HOUR,
	DAY,
	MONTH,
	YEAR
} CPDFtimeTypes;

typedef enum {
	CPDF_Root,
	CPDF_Catalog,
	CPDF_Outlines,
	CPDF_Pages,
	CPDF_Page,
	CPDF_Contents,
	CPDF_ProcSet,
	CPDF_Annots,
	CPDF_Info
} CPDFobjTypes;


/* 2 x 3 matrix for CTM */
typedef struct {
    float a; float b;		/*  a   b   0  */
    float c; float d;		/*  c   d   0  */
    float x; float y;		/*  x   y   1  */
} CPDFctm;

typedef struct {
	unsigned long   magic_number;	/* to check stream validity */
	char  *buffer;    	/* pointer to buffer's beginning */
	/* char  *buf_ptr; */		/* current buffer pointer */
	int   count;  		/* # of bytes currently in buffer */
	int   bufSize;   	/* Total size of buffer -- bufSize expands as needed. */
} CPDFmemStream;

typedef struct _cpdf_intlist CPDFintList;
struct _cpdf_intlist {
    int	        value;		/* integer value */
    CPDFintList *next;		/* NULL here terminates the list */
};

/* For passing as "align" to cpdf_textBox() and cpdf_rawTextBox() */
#define	TBOX_LEFT		0
#define	TBOX_CENTER		1
#define	TBOX_RIGHT		2
#define	TBOX_JUSTIFY		3

/* TextBox attribute structure: for cpdf_textBox() and cpdf_rawTextBox() */
typedef struct {
	int alignmode;		/* one of the modes above (default TBOX_LEFT) */
	int NLmode;		/* if non-zero, NL is is a line break, if 0 reformatted (default 0) */
	float paragraphSpacing;	/* extra space between paragraphs (default = 0.0) */
	int noMark;		/* if non-zero, text is actually not printed (for fitting text to box) */
} CPDFtboxAttr;

/* External font for embedding: Holds font data stream (essentialy content of PFB file) */
typedef struct extfont_  CPDFextFont;
typedef struct fontinfo_ CPDFfontInfo;

typedef struct fontdesc_ {
    int objIndex;
    int doneAlready;
    int ascent;
    int capHeight;
    int descent;
    unsigned long flags;
    int fontBBox[4];
    char *fontName;		/* For external font, this will contain read font name from PFM */
    float italicAngle;		/* Some fonts uses fractional angle. Distiller uses floats */
    int stemV;
    /* ------------ Everything above required ---------------------- */
    int stemH;			/* optional - needed for CJK */
    int xHeight;		/* AFM files have this, and Distiller puts this */
    int missingWidth;		/* for CJK */
    int leading;		/* for CJK */
    int maxWidth;		/* for CJK */
    int avgWidth;		/* for CJK */
    char *style;		/* for CJK */
    CPDFextFont *extFont;	/* for later font embedding */
} CPDFfontDescriptor;

/* for CIDkeyed Font for CJK */
typedef struct fontdescendant_ {
    int objIndex;
    int doneAlready;
    int DW;
    char *cidSysInfo;
    char *cidFontWidth;
} CPDFdescendantFont;

/* This is CPDFextFont struct */
struct extfont_ {
    int objIndex;
    int doneAlready;
    CPDFextFont *prev;	/* */
    CPDFfontDescriptor *parentDesc;
    CPDFfontInfo *parentFInfo;
    int type;			/* Font type: 1=Type1 */
    int subtype;
    char *basefont;		/* Stub fontname specified in cpdf_setFont() and fontmap.lst file. */
    long length;		/* length of font stream */
    long length1;		/* ASCII portion of Type1 */
    long length2;		/* BINARY portion of Type1 */
    long length3;		/* ASCII portion of Type1 */
    CPDFmemStream *fontMemStream;	/* memory stream containing font stream */
};	/* This is CPDFextFont */


/* This is CPDFfontInfo  struct */
struct fontinfo_ {
    int objIndex;
    char *name;			/* This is the name like "Fcpdf1" used in Page Content */
    char *baseFont;		/* For external font, this will contain real font name read from PFM */
    char *stubFont;		/* font name as specified in cpdf_setFont() */
    char *encoding;
    /* --- added stuff below starting with version ClibPDF v1.10 --- */
    int  afmIndex;
    int  encodeIndex;	/* 0=WinAnsiEncoding, 1=MacRomanEncoding, 2=MacExpertEncoding, 3=StandardEncoding */
    int  descLevel;	/* 0=Base14, 1=Metrics, 2=Embedded, 3=CJK built-in */
    /* char *subtype; */	/* Type1 only for now -- hard coded in cpdfInit.c: _cpdf_WriteFont() */
    short  firstchar;
    short  lastchar;
    short  *width;	/* array of char widths: always of 256 elements [0..255] */
    CPDFfontDescriptor *fontDesc;	/* font descriptor */
    CPDFdescendantFont *descendFont;	/* for CJK font */
};	/* This is CPDFfontInfo */


typedef struct {
    int pageMode;	/* This really belongs directly to Catalog obj, but here for convenience */
    int hideToolbar;	/* when YES, tool bar in viewer will be hidden */
    int hideMenubar;	/* when YES, menu bar in viewer will be hidden */
    int hideWindowUI;	/* when YES, GUI elements in doc window will be hidden */
    int fitWindow;	/* when YES, instructs the viewer to resize the doc window to the page size */
    int centerWindow;	/* when YES, instructs the viewer to move the doc window to the screen's center */
    int pageLayout;	/* Specifies 1-page display or two-page columnar displays */
    int nonFSPageMode;	/* Specifies pageMode coming out of the full-screen display mode */
} CPDFviewerPrefs;

/* Values for pageMode and nonFSPageMode */
#define	PM_NONE		0	/* default - neither outline nor thumbnails will be visible */
#define PM_OUTLINES	1	/* open the doc with outline visible */
#define PM_THUMBS	2	/* open the doc with thumbnails visible */
#define PM_FULLSCREEN	3	/* open the doc in full screen mode */

#define PL_SINGLE	0	/* default - one page at a time */
#define PL_1COLUMN	1	/* display pages in one column */
#define PL_2LCOLUMN	2	/* 2-column display, with odd pages on the left */
#define PL_2RCOLUMN	3	/* 2-column display, with odd pages on the right */


/* Image file type */
#define JPEG_IMG		0
#define G4FAX_IMG		1
#define G3FAX_IMG		2
#define TIFF_IMG		3
#define GIF_IMG			4
#define CPDF_IMG		5	/* our custom image format for fast loading into PDF */
#define RAW_IMG			6	/* 1999-11-19 */

/* Flags for cpdf_importImage() and cpdf_placeImageData() */
#define IM_GSAVE		1
#define IM_IMAGEMASK		2
#define IM_INVERT		4

typedef struct {
	int objIndex;
	char *name;		/* Im0, Im1, etc */
	int type;		/* image file type */
	int imagemask;		/* 1 if this image is /ImageMask, 0 if not */
	int invert;		/* invert B/W or gray */
	int process;		/* M_SOF# -- jpeg process or TIFF compression type  */
	int NumStrips;		/* number of strips (only for TIFF) */
	int imgSelector;	/* index (0..N-1) to an image in multi-page file (TIFF) */
	int photometric;	/* photometric interpretation for tiff */
	int orientation;	/* orientation tag in TIFF file (not always present, normally 1) */
	int width;		/* # of pixels horizontal */
	int height;		/* # of pixels vertical */
	int resUnit;		/* resolution unit, e.g., DPI */
	float resX;		/* resolution for horizontal dimension */
	float resY;		/* resolution for height dimension */
	int ncomponents;	/* # of color components */
	int bitspersample;	/* bits per sample */
	long filesize;		/* # of bytes in file or data */
	char *filepath;		/* path to image file */
	void *data;		/* image data for non-JPEG images */
} CPDFimageInfo;

/* For future use in cpdf_placeImageData() and cpdf_rawPlaceImageData() */
typedef struct {
	int reserved;
} CPDFimgAttr;

/* Annotation, hyperlink info object */
#define	ANNOT_TEXT	0
#define	ANNOT_URL	1
#define ANNOT_GOTO	2
#define ANNOT_ACTION	3

/* Annotation flags, may be OR'ed together */
#define	AF_INVISIBLE	0x0001
#define	AF_HIDDEN	0x0002
#define AF_PRINT	0x0004
#define AF_NOZOOM	0x0008
#define	AF_NOROTATE	0x0010
#define AF_NOVIEW	0x0020
#define	AF_READONLY	0x0040

typedef struct {
    int objIndex;
    int type;		/* annotation or link type, see above defines */
    int page;
    int flags;		/* see PDF spec 1.3, page 83-84 */
    float xLL;
    float yLL;
    float xUR;
    float yUR;			/* annotation box */
    float r;			/* RGB color components */
    float g;
    float b;
    char *border_array;		/* array of 3 numbers plus optional dash array */
    char *BS;			/* content of BS (border style dict) overrides border_array if not NULL */
    char *content_link;		/* annotation text content, link URI specification, or action dictionary */
    long content_len;
    char *annot_title;		/* annotation box title */
    long title_len;
} CPDFannotInfo;

typedef struct {
    int flags;
    char *border_array;		/* as above */
    char *BS;			/* content of BS (border style dict) overrides border_array if not NULL */
    float r;			/* RGB color components */
    float g;
    float b;
} CPDFannotAttrib;

/* Structure for outline (bookmark) linked list */
typedef struct _cpdf_outline CPDFoutlineEntry;
struct _cpdf_outline {
    int objIndex;		/* serialized object index (for xref) */
    int count;			/* total number of sub entries (descendants) */
    int dest_page;		/* page number set by cpdf_pageInit(), negative for Action outline */
    int open;			/* if zero, all subsections under this entry will be closed */
    char *dest_attr_act_dict;	/* Destination spec after "3 0 R" part, or content of Action dictionary */
    char *title;		/* title string */
    long  title_len;		/* byte length of title above */
    CPDFoutlineEntry *parent;	/* pointer to parent outline entry */
    CPDFoutlineEntry *prev;	/* pointer to previous outline entry */
    CPDFoutlineEntry *next;	/* pointer to next outline entry */
    CPDFoutlineEntry *first;	/* pointer to first outline entry */
    CPDFoutlineEntry *last;	/* pointer to last outline entry */
};

/* Outline (book mark) destination modes (see page 95, Table 6.20 of THE PDF Reference) */
#define DEST_NULL		0	/* keep current display location and zoom */
#define DEST_XYZ		1	/* /XYZ left top zoom (equivalent to above) */
#define DEST_FIT		2	/* /Fit  */
#define DEST_FITH		3	/* /FitH top */
#define DEST_FITV		4	/* /FitV left */
#define DEST_FITR		5	/* /FitR left bottom right top */
#define DEST_FITB		6	/* /FitB   (fit bounding box to window) PDF-1.1 */
#define DEST_FITBH		7	/* /FitBH top   (fit width of bounding box to window) PDF-1.1 */
#define DEST_FITBV		8	/* /FitBV left   (fit height of bounding box to window) PDF-1.1 */
#define DEST_Y			9	/* /XYZ null y null */

#define OL_SUBENT		1	/* add as sub-entry */
#define OL_SAME			0	/* add outline at the same level */
#define OL_OPEN			1	/* outline open */
#define OL_CLOSED		0	/* outline closed (all subentries under this) */

#define  DOMAIN_MAGIC_NUMBER	0xdada3333

typedef struct _cpdf_doc CPDFdoc;	/* see blow for struct _cpdf_doc */

typedef struct _cpdf_domain {
	unsigned long magic;		/* domain magic number */
	CPDFdoc *pdfdoc;		/* PDF document this domain belongs to */
	/* struct _cpdf_domain *parent; */ /* pointer to parent domain, null if the parent is the entire page. */
	float xloc, yloc;		/* coordinate of lower-left corner of this domain in points */
	float width, height;		/* width and height of this domain in points */
	float xvalL, xvalH;		/* low- and high-limit values of the X axis in domain unit */
	float yvalL, yvalH;		/* low- and high-limit values of the Y axis of the domain */
	struct tm xvTL, xvTH;		/* low- and high-limit values for time X axis of the domain */
	int xtype, ytype;		/* axis flags: 0=linear, 1=log, 2=time */
	int polar;			/* reserved */
	int enableMeshMajor;
	int enableMeshMinor;
	char *meshDashMajor;		/* dash array spec for major mesh lines */
	char *meshDashMinor;		/* dash array spec for minor mesh lines */
	float meshLineWidthMajor;
	float meshLineWidthMinor;
	float meshLineColorMajor[3];
	float meshLineColorMinor[3];
	/* int numChildren; */		/* # of child domains */
	/* struct _cpdf_domain **children; */	/* array of pointers to child plot domains */
	/* for linear X axis */
	float xvalFirstMeshLinMajor;	/* value of first major mesh line */
	float xvalFirstMeshLinMinor;	/* value of first minor mesh line */
	float xmeshIntervalLinMajor;	/* mesh interval for linear axis */
	float xmeshIntervalLinMinor;
	/* for linear Y axis */
	float yvalFirstMeshLinMajor;	/* value of first major mesh line */
	float yvalFirstMeshLinMinor;	/* value of first minor mesh line */
	float ymeshIntervalLinMajor;	/* mesh interval for linear axis */
	float ymeshIntervalLinMinor;
} CPDFplotDomain;

#define AXIS_MAGIC_NUMBER   0xafafafaf

typedef struct {
	unsigned long magic;		/* axis magic number */
	CPDFplotDomain *plotDomain;	/* pointer to parent domain, null if the parent is the entire page. */
	float angle;			/* angle of axis, 0.0 for X-axis, 90.0 for Y-axis */
	int   type;			/* 0=linear, 1=logarithmic, 2=time */
	float xloc, yloc;		/* location of the beginning of axis relative to domain's xloc, yloc */
	float length;			/* length of axis in points */
	float axisLineWidth;		/* width of axis line in points */
	float valL, valH;		/* high and low values of the axis (for numbering and ticks) */
	struct tm vTL, vTH;		/* low- and high-limit values for time X axis of the domain */
	int   ticEnableMajor;		/* 0=No tics, 1=Enabled (regular), 2=Free style (list provided) */
	int   ticEnableMinor;		/* 0=No tics, 1=Enabled (regular), 2=Free style (list provided) */
	float ticLenMajor;		/* length of major ticks in points */
	float ticLenMinor;		/* length of minor ticks in points */
	float tickWidthMajor;		/* linewidth for major ticks */
	float tickWidthMinor;
	int   ticPosition;		/* tick position: 0=CWside (below X), 1=Middle, 2=CCWside (above X) */
	int   numPosition;		/* number (label) position: 0=CWside (below X),  2=CCWside (above X) */

	int   numEnable;		/* 0=No #s, 1=Enabled (regular), 2=Free style (list provided) */
	float ticNumGap;		/* gap (in points) between tic end and number */
	float numFontSize;
	int   useMonthName;		/* non-zero for using month names rather than numbers */
	int   use2DigitYear;		/* non-zero value will use 2-digit year for display */
	int   horizNumber;		/* number text is horizontal if non-zero */
	int   numStyle;			/* axis number style regular, exponent etc. */
	/* int   numPrecision; */	/* # of digits after decimal point */
	char  *numFormat;		/* set axis number format */
	char  *numFontName;		/* Font name for numbers */
	char  *numEncoding;		/* Font encoding for numbers */

	float numLabelGap;		/* gap (in points) between number and axis label */
	float labelFontSize;
	int   horizLabel;
	char  *labelFontName;		/* Font name for axis label */
	char  *labelEncoding;		/* Font encoding for axis label */
	char  *axisLabel;		/* Axis label string, if NULL, no label is shown */

	/* for linear axis */
	float valFirstTicLinMajor;	/* value of first major tick */
	float valFirstTicLinMinor;	/* value of first minor tick */
	float ticIntervalLinMajor;	/* tick interval for linear axis, Major ticks will be numbered. */
	float ticIntervalLinMinor;
	/* for log axis */
	int   ticSelectorLog;		/* log axis tick enable mask */
	int   numSelectorLog;		/* log axis number enable mask */
	/* for time axis */
	int   lastMin;
	int   lastHour;
	int   lastDay;
	int   lastMonth;
	int   lastYear;
} CPDFaxis;


/* Describes attributes and resources used on each page.
   Each one of these corresponds to one Page object.
*/
typedef struct {
	int pagenum;			/* page number from cpdf_pageInit() call */
	int objIndex;			/* obj index of itself */
	int parent;			/* obj index of parent Pages object */
	int contents;			/* obj index of Contents stream for this page */
	int status;			/* bit flags: b-1=compressed, b-0=closed */
	long length_compressed;		/* length (bytes) of compressed page content stream */
	char *compressedStream;		/* compressed page content stream */
	CPDFmemStream  *pageMemStream;	/* uncompressed page content memory stream */
	CPDFplotDomain *defDomain;	/* default domain for this page */
	int orientation;		/* page orientation */
	int npFont;			/* # of fonts used on this page */
	int npImage;			/* # of images used on this page */
	int npAnnot;			/* # of annotations on for this page */
	CPDFintList *fontIdx;		/* linked list of fonts containing indices into fontInfos[] */
	CPDFintList *imageIdx;		/* linked list of images containing indices into imageInfos[] */
	CPDFintList *annotIdx;		/* linked list of annotations and links as above */
        char *mediaBox;			/* MediaBox () */
        char *cropBox;			/* CropBox () */
	FILE *fppage;			/* file stream pointer for this page */
	char *contentfile;		/* file for Content stream -- only when memory strem is not used */
	float duration;			/* if > 0.0 in seconds, the page will be displayed for that period */
	char *transition;		/* transition effects */
} CPDFpageInfo;


typedef struct {
	int nMaxPages;
	int nMaxFonts;
	int nMaxImages;
	int nMaxAnnots;
	int nMaxObjects;
} CPDFdocLimits;


typedef struct {
	int objIndex;			/* obj index of itself */
	int version;			/* currently 1 */
	int revision;			/* currently 2 */
	unsigned long permissions;
	char *filter;			/* currently points to static string "/Standard" */
	char owner[32];			/* data for owner password */
	char user[32];			/* data for user password */
} CPDFencryptDict;

/* Error handler */
typedef void (*CPDFglobalErrorHandler)(int level, const char *module, const char *fmt, va_list ap);
typedef void (*CPDFerrorHandler)(CPDFdoc *pdf, int level, const char *module, const char *fmt, va_list ap);

/* PDF document structure for allowing simultaneous multiple document creation */ 
struct _cpdf_doc {
    int docID;					/* document ID (thread ID ?) */
    int ps_pdf_mode;				/* PDF=0, EPDF=1, PS=2, EPS=3 (not used) */
    int pdfLevelMaj;				/* PDF level, do not use operators beyond these */
    int pdfLevelMin;
    int init_check;				/* used in cpdfInit.c */
    char **monthName;				/* see cpdfAxis.c */
    char *monthNameEncoding;
    float defdomain_unit;			/* unit for default domain */
    int display_rotation;
    int useStandardOutput;			/* send output to stdout if non-zero */
    int linearizeON;				/* linearized PDF */
    int encryptionON;				/* encryption to set security options */
    CPDFencryptDict cryptDict;			/* encryption dictionary */
    int compressionON;				/* compress stream */
    char *compress_command;			/* command for LZW compression */
    char *streamFilterList;			/* for PDF stream /Filter spec */
    int launchPreview;				/* launch viewer application on the output file */
    int filename_set;				/* flag indicating if output filename is set explicitly */
    int inTextObj;				/* flag indicating within Text block between BT ET */
    CPDFviewerPrefs viewerPrefs;		/* viewer preferences such as whether outline should be open */
    CPDFerrorHandler _cpdfErrorHandler;		/* error handler */
    CPDFplotDomain *defaultDomain;		/* default plot domain */
    CPDFplotDomain *currentDomain;		/* current plot domain */
    float xcurr, ycurr;				/* current raw coordinate.  Used in cpdfArc.c, cpdfRawPlot.c */
    float x2points, y2points;			/* scaling factor for current domain */
    double xLlog, xHlog, yLlog, yHlog;		/* scaling factor for current domain (logarithmic) */
    int nMaxFonts;				/* max number of fonts as a variable - can't be changed */
    int numFonts;				/* number of fonts used: Unique for each fontname+encoding */
    int numExtFonts;				/* number of external fonts embedded */
    CPDFextFont *extFontList;			/* point to linked list of external fonts */
    CPDFfontInfo *fontInfos;			/* array of font infos */
    int currentFont;				/* current font index (into fontInfos[]) */
    int inlineImages;				/* in-line image count */
    int nMaxImages;				/* maximum number of unique images  - can't be changed */
    int numImages;
    CPDFimageInfo *imageInfos;
    int imageFlagBCI;				/* bit-0 (/ImageB), bit-1 (/ImageC), bit-2 (/ImageI) */
    int imageSelector;				/* image index (0 .. N-1) for multi-page TIFF file */
    int numOutlineEntries;			/* total # of outline (book mark) entries */
    CPDFoutlineEntry *firstOLentry;		/* pointer to first outline entry */
    CPDFoutlineEntry *lastOLentry;		/* pointer to last outline entry */

    int   hex_string;				/* current string mode: If non-zero HEX */
    float font_size;				/* current font size and info below */
    float word_spacing;
    float char_spacing;
    float text_rise;
    float horiz_scaling;			/* text horizontal scaling in percent */
    float text_leading;
    CPDFctm textCTM; 				/* { 1.0, 0.0, 0.0, 1.0, 0.0, 0.0 }; */
    int textClipMode;

    int usePDFMemStream;			/* if non-zero use memory stream for PDF generation */
    CPDFmemStream *pdfMemStream;		/* memory stream for PDF file that is currently active */
    int useContentMemStream;			/* if non-zero use memory stream for Content */
    CPDFmemStream *currentMemStream;		/* memory stream for Content that is currently active */
    int currentPage;				/* current page number that is being drawn */
    int nMaxPages;				/* maximum number of pages */
    int numPages;				/* number of pages - may be greater than actual # of pages */
    CPDFpageInfo *pageInfos;			/* array of pageInfo structure for all pages (alloc nMaxPages+1) */
    int numKids;				/* actual # of pages counted for Pages object */
    int *kidsIndex;				/* object index list for kids to be written to Pages object */
    CPDFmemStream *scratchMem;			/* use this as non-overflowing scratch pad */
    FILE *fpcg; 				/* Output file */
    FILE *fpcontent;				/* Content stream (need length) */
    int  nMaxAnnots;				/* maximum number of annotations  - can't be changed */
    int  numAnnots;				/* count of annotations */
    CPDFannotInfo *annotInfos;			/* array of annotInfo structure for all annotations */
    char mediaBox[64];				/* MediaBox for current page*/
    char cropBox[64];				/* CropBox for current page */
    long currentByteCount;			/* # of bytes written, or offset of next object */
    char creator_name[64];			/* Info: set it by cpdf_setCreator() */
    char file_title[64];			/* Info: title of PDF file */
    char file_subject[64];			/* Info: subject of PDF file */
    char file_keywords[128];			/* Info: keywords */
    char username[64];				/* user name */
    char filenamepath[1024];
    char contentfile[1024];
    char pfm_directory[1024];
    char pfb_directory[1024];
    char fontmapfile[1024];
    int  nMaxObjects;				/* maximum number of objects for xref */
    long *objByteOffset;			/* offset into object number N */
    int  *objIndex;				/* object index for selected objects */
    long startXref;				/* offset of xref */
    char spbuf[8192];				/* scratch buffer for sprintf */
/* These are from cpdfAxis.c */
/* This is for transfering edge of numbers to figure out how far move the axis label off */
    float numEdgeY;
/* ohter temprary data area */
    double xa2points, xaLlog, xaHlog;
    CPDFaxis *anAx2;
};		/* this is CPDFdoc struct -- see above CPDFplotDomain */



/* extern CPDFplotDomain *defaultDomain; */	/* default plot domain */
/* extern CPDFplotDomain *currentDomain; */	/* current plot domain */

/* Public API functions ---------------------------------------------------------- */
char *cpdf_version(void);
char *cpdf_platform(void);

CPDFerrorHandler cpdf_setErrorHandler(CPDFdoc *pdf, CPDFerrorHandler handler);
void cpdf_Error(CPDFdoc *pdf, int level, const char* module, const char* fmt, ...);
CPDFglobalErrorHandler cpdf_setGlobalErrorHandler(CPDFglobalErrorHandler handler);
void cpdf_GlobalError(int level, const char* module, const char* fmt, ...);


void cpdf_setGlobalDocumentLimits(int maxPages, int maxFonts, int maxImages, int maxAnnots, int maxObjects);
void cpdf_setViewerPreferences(CPDFdoc *pdf, CPDFviewerPrefs *vP);

CPDFdoc *cpdf_open(int pspdf, CPDFdocLimits *docLimits);
void cpdf_enableCompression(CPDFdoc *pdf, int cmpON);
void cpdf_useContentMemStream(CPDFdoc *pdf, int flag);
void cpdf_setCompressionFilter(CPDFdoc *pdf, const char *command, const char *filterlist);
void cpdf_setDefaultDomainUnit(CPDFdoc *pdf, float defunit);
void cpdf_init(CPDFdoc *pdf);
int cpdf_pageInit(CPDFdoc *pdf, int pagenum, int rot, const char *mediaboxstr, const char *cropboxstr);
void cpdf_finalizeAll(CPDFdoc *pdf);
int cpdf_savePDFmemoryStreamToFile(CPDFdoc *pdf, const char *file);
char *cpdf_getBufferForPDF(CPDFdoc *pdf, int *length);
void cpdf_close(CPDFdoc *pdf);
int cpdf_launchPreview(CPDFdoc *pdf);
int cpdf_openPDFfileInViewer(CPDFdoc *pdf, const char *pdffilepath);
void cpdf_setCreator(CPDFdoc *pdf, char *pname);
void cpdf_setTitle(CPDFdoc *pdf, char *pname);
void cpdf_setSubject(CPDFdoc *pdf, char *pname);
void cpdf_setKeywords(CPDFdoc *pdf, char *pname);
int  cpdf_comments(CPDFdoc *pdf, char *comments);

/* Page related public functions */
int cpdf_setCurrentPage(CPDFdoc *pdf, int page);
void cpdf_finalizePage(CPDFdoc *pdf, int page);
void cpdf_setPageSize(CPDFdoc *pdf, const char *mboxstr, const char *cboxstr);	/* e.g., "0 0 612 792" */
void cpdf_setBoundingBox(CPDFdoc *pdf, int LLx, int LLy, int URx, int URy);
void cpdf_setMediaBox(CPDFdoc *pdf, int LLx, int LLy, int URx, int URy);
void cpdf_setCropBox(CPDFdoc *pdf, int LLx, int LLy, int URx, int URy);

/* page duration and transition */
void cpdf_setPageDuration(CPDFdoc *pdf, float seconds);
int cpdf_setPageTransition(CPDFdoc *pdf, int type, float duration, float direction, int HV, int IO);


/* Annotation and hyper link functions */
void cpdf_setAnnotation(CPDFdoc *pdf, float xll, float yll, float xur, float yur,
	const char *title, const char *annot, CPDFannotAttrib *attr);
void cpdf_setActionURL(CPDFdoc *pdf, float xll, float yll, float xur, float yur,
	const char *linkspec, CPDFannotAttrib *attr);
void cpdf_setLinkAction(CPDFdoc *pdf, float xll, float yll, float xur, float yur,
	const char *action_dict, CPDFannotAttrib *attr);
void cpdf_setLinkGoToPage(CPDFdoc *pdf, float xll, float yll, float xur, float yur, int page,
	const char *fitmode, CPDFannotAttrib *attr);
void cpdf_rawSetLinkGoToPage(CPDFdoc *pdf, float xll, float yll, float xur, float yur, int page,
	const char *fitmode, CPDFannotAttrib *attr);
void cpdf_rawSetAnnotation(CPDFdoc *pdf, float xll, float yll, float xur, float yur,
	const char *title, const char *annot, CPDFannotAttrib *attr);
void cpdf_rawSetActionURL(CPDFdoc *pdf, float xll, float yll, float xur, float yur,
	const char *linkspec, CPDFannotAttrib *attr);
void cpdf_rawSetLinkAction(CPDFdoc *pdf, float xll, float yll, float xur, float yur,
	const char *action_dict, CPDFannotAttrib *attr);
int cpdf_includeTextFileAsAnnotation(CPDFdoc *pdf, float xll, float yll, float xur, float yur,
	const char *title, const char *filename, CPDFannotAttrib *attr);
int cpdf_rawIncludeTextFileAsAnnotation(CPDFdoc *pdf, float xll, float yll, float xur, float yur,
	const char *title, const char *filename, CPDFannotAttrib *attr);

/* ==== Text and Font functions ================================================ */
void cpdf_beginText(CPDFdoc *pdf, int clipmode);
void cpdf_endText(CPDFdoc *pdf);

/* convenient text functions */
void cpdf_text(CPDFdoc *pdf, float x, float y, float orientation, const char *textstr);
void cpdf_rawText(CPDFdoc *pdf, float x, float y, float orientation, const char *textstr);
void cpdf_textAligned(CPDFdoc *pdf, float x, float y, float orientation, int centmode, const char *textstr);
void cpdf_rawTextAligned(CPDFdoc *pdf, float x, float y, float orientation, int centmode, const char *textstr);
/* NOTE: TextBox char *text is not "const", the text gets modified */
char *cpdf_rawTextBoxY(CPDFdoc *pdf, float xl, float yl, float width, float height, float angle,
			float linespace, float *yend, CPDFtboxAttr *tbattr, char *text);
char *cpdf_textBoxY(CPDFdoc *pdf, float xl, float yl, float width, float height, float angle,
			float linespace, float *yend, CPDFtboxAttr *tbattr, char *text);
char *cpdf_textBox(CPDFdoc *pdf, float xl, float yl, float width, float height, float angle,
			float linespace, CPDFtboxAttr *tbattr, char *text);
char *cpdf_rawTextBox(CPDFdoc *pdf, float xl, float yl, float width, float height, float angle,
			float linespace, CPDFtboxAttr *tbattr, char *text);
float cpdf_textBoxFit(CPDFdoc *pdf, float xl, float yl, float width, float height, float angle,
		float inifontsize, float fsdecrement, float linespace, CPDFtboxAttr *tbattr, char *text);
float cpdf_rawTextBoxFit(CPDFdoc *pdf, float xl, float yl, float width, float height, float angle,
		float inifontsize, float fsdecrement, float linespace, CPDFtboxAttr *tbattr, char *text);

/* primitive PDF text operator functions */
void cpdf_textShow(CPDFdoc *pdf, const char *txtstr);
void cpdf_textCRLFshow(CPDFdoc *pdf, const char *txtstr);
void cpdf_textCRLF(CPDFdoc *pdf);
void cpdf_setNextTextLineOffset(CPDFdoc *pdf, float x, float y);
void cpdf_rawSetNextTextLineOffset(CPDFdoc *pdf, float x, float y);
void cpdf_setTextRise(CPDFdoc *pdf, float rise);
void cpdf_setTextRenderingMode(CPDFdoc *pdf, int mode);
void cpdf_setTextMatrix(CPDFdoc *pdf, float a, float b, float c, float d, float x, float y);
void cpdf_concatTextMatrix(CPDFdoc *pdf, float a, float b, float c, float d, float x, float y);
void cpdf_rotateText(CPDFdoc *pdf, float degrees);
void cpdf_skewText(CPDFdoc *pdf, float alpha, float beta);
void cpdf_setTextPosition(CPDFdoc *pdf, float x, float y);
void cpdf_rawSetTextPosition(CPDFdoc *pdf, float x, float y);
void cpdf_setTextLeading(CPDFdoc *pdf, float leading);
void cpdf_setHorizontalScaling(CPDFdoc *pdf, float scale);
void cpdf_setCharacterSpacing(CPDFdoc *pdf, float spacing);
void cpdf_setWordSpacing(CPDFdoc *pdf, float spacing);
char *cpdf_escapeSpecialChars(const char *instr);
char *cpdf_escapeSpecialCharsBinary(const char *instr, long lengthIn, long *lengthOut);

/* Fonts */
int cpdf_setFont(CPDFdoc *pdf, const char *basefontname, const char *encodename, float size);
void cpdf_setFontDirectories(CPDFdoc *pdf, const char *pfmdir, const char *pfbdir);
void cpdf_setFontMapFile(CPDFdoc *pdf, const char *mapfile);
float cpdf_stringWidth(CPDFdoc *pdf, const unsigned char *str);
float cpdf_capHeight(CPDFdoc *pdf);
float cpdf_ascent(CPDFdoc *pdf);
float cpdf_descent(CPDFdoc *pdf);
void multiplyCTM(CPDFctm *T, const CPDFctm *S);

/* Plot Domain functions */
CPDFplotDomain *cpdf_createPlotDomain(float x, float y, float w, float h,
			float xL, float xH, float yL, float yH,
			int xtype, int ytype, int reserved);
CPDFplotDomain *cpdf_createTimePlotDomain(float x, float y, float w, float h,
			struct tm *xTL, struct tm *xTH, float yL, float yH,
			int xtype, int ytype, int reserved);
void cpdf_freePlotDomain(CPDFplotDomain *aDomain);
CPDFplotDomain *cpdf_setPlotDomain(CPDFdoc *pdf, CPDFplotDomain *aDomain);
void cpdf_clipDomain(CPDFplotDomain *aDomain);
void cpdf_fillDomainWithGray(CPDFplotDomain *aDomain, float gray);
void cpdf_fillDomainWithRGBcolor(CPDFplotDomain *aDomain, float r, float g, float b);
void cpdf_setMeshColor(CPDFplotDomain *aDomain, float meshMajorR, float meshMajorG, float meshMajorB,
		       float meshMinorR, float meshMinorG, float meshMinorB);
void cpdf_drawMeshForDomain(CPDFplotDomain *aDomain);
void cpdf_setLinearMeshParams(CPDFplotDomain *aDomain, int xy, float mesh1ValMajor, float intervalMajor,
					      float mesh1ValMinor, float intervalMinor);
void cpdf_suggestMinMaxForLinearDomain(float vmin, float vmax, float *recmin, float *recmax);
void cpdf_suggestLinearDomainParams(float vmin, float vmax, float *recmin, float *recmax,
		float *tic1ValMajor, float *intervalMajor,
		float *tic1ValMinor, float *intervalMinor);
float x_Domain2Points(CPDFdoc *pdf, float x);
float y_Domain2Points(CPDFdoc *pdf, float y);
float x_Points2Domain(CPDFdoc *pdf, float xpt);
float y_Points2Domain(CPDFdoc *pdf, float ypt);
void cpdf_suggestTimeDomainParams(struct tm *xTL, struct tm *xTH, struct tm *recTL, struct tm *recTH);

/* Axis functions */
CPDFaxis *cpdf_createAxis(float angle, float axislength, int typeflag, float valL, float valH);
CPDFaxis *cpdf_createTimeAxis(float angle, float axislength, int typeflag, struct tm *vTL, struct tm *vTH);
void cpdf_freeAxis(CPDFaxis *anAx);
void cpdf_drawAxis(CPDFaxis *anAx);
void cpdf_attachAxisToDomain(CPDFaxis *anAx, CPDFplotDomain *domain, float x, float y);
void cpdf_setAxisLineParams(CPDFaxis *anAx, float axLineWidth, float ticLenMaj, float ticLenMin,
				float tickWidMaj, float tickWidMin);
void cpdf_setTicNumEnable(CPDFaxis *anAx, int ticEnableMaj, int ticEnableMin, int numEnable);
void cpdf_setAxisTicNumLabelPosition(CPDFaxis *anAx, int ticPos, int numPos, int horizNum, int horizLab);
void cpdf_setAxisNumberFormat(CPDFaxis *anAx, const char *format, const char *fontName, float fontSize);
void cpdf_setAxisNumberFormat2(CPDFaxis *anAx, const char *format,
		const char *fontName, float fontSize, const char *encoding);
void cpdf_setTimeAxisNumberFormat(CPDFaxis *anAx, int useMonName, int use2DigYear,
		const char *fontName, float fontSize);
void cpdf_setTimeAxisNumberFormat2(CPDFaxis *anAx, int useMonName, int use2DigYear,
		const char *fontName, float fontSize, const char *encoding);
void cpdf_setAxisLabel(CPDFaxis *anAx, const char *labelstring, const char *fontName,
	const char *encoding, float fontSize);
void cpdf_setLinearAxisParams(CPDFaxis *anAx, float tic1ValMajor, float intervalMajor,
					      float tic1ValMinor, float intervalMinor);
void cpdf_setLogAxisTickSelector(CPDFaxis *anAx, int ticselect);
void cpdf_setLogAxisNumberSelector(CPDFaxis *anAx, int numselect);

float vAxis2Points(CPDFdoc *pdf, float x);


/* Drawing and path constructions functions */
/* These functions use current (scaled linear or log) coordinate system for (x, y) */
void cpdf_moveto(CPDFdoc *pdf, float x, float y);
void cpdf_lineto(CPDFdoc *pdf, float x, float y);
void cpdf_rmoveto(CPDFdoc *pdf, float x, float y);
void cpdf_rlineto(CPDFdoc *pdf, float x, float y);
void cpdf_curveto(CPDFdoc *pdf, float x1, float y1, float x2, float y2, float x3, float y3);
void cpdf_rect(CPDFdoc *pdf, float x, float y, float w, float h);
void cpdf_rectRotated(CPDFdoc *pdf, float x, float y, float w, float h, float angle);
void cpdf_quickCircle(CPDFdoc *pdf, float xc, float yc, float r);		/* center (x,y) and radius r */
void cpdf_arc(CPDFdoc *pdf, float x, float y, float r, float sangle, float eangle, int moveto0);
void cpdf_circle(CPDFdoc *pdf, float x, float y, float r);

/* multi-threaded */
void cpdf_arc(CPDFdoc *pdf, float x, float y, float r, float sangle, float eangle, int moveto0);
void cpdf_circle(CPDFdoc *pdf, float x, float y, float r);

/* These use raw, point-based coordinate system for (x,y) */
void cpdf_rawMoveto(CPDFdoc *pdf, float x, float y);
void cpdf_rawLineto(CPDFdoc *pdf, float x, float y);
void cpdf_rawRmoveto(CPDFdoc *pdf, float x, float y);
void cpdf_rawRlineto(CPDFdoc *pdf, float x, float y);
void cpdf_rawCurveto(CPDFdoc *pdf, float x1, float y1, float x2, float y2, float x3, float y3);
void cpdf_rawRect(CPDFdoc *pdf, float x, float y, float w, float h);
void cpdf_rawRectRotated(CPDFdoc *pdf, float x, float y, float w, float h, float angle);
void cpdf_rawQuickCircle(CPDFdoc *pdf, float xc, float yc, float r);	/* center (x,y) and radius r */
void cpdf_rawArc(CPDFdoc *pdf, float x, float y, float r, float sangle, float eangle, int moveto0);
void cpdf_rawCircle(CPDFdoc *pdf, float xc, float yc, float r);

void cpdf_rawArc(CPDFdoc *pdf, float x, float y, float r, float sangle, float eangle, int moveto0);
void cpdf_rawCircle(CPDFdoc *pdf, float xc, float yc, float r);

/* Operations on current path */
void cpdf_closepath(CPDFdoc *pdf);
void cpdf_stroke(CPDFdoc *pdf);
void cpdf_fill(CPDFdoc *pdf);
void cpdf_eofill(CPDFdoc *pdf);
void cpdf_fillAndStroke(CPDFdoc *pdf);
void cpdf_eofillAndStroke(CPDFdoc *pdf);
void cpdf_clip(CPDFdoc *pdf);
void cpdf_eoclip(CPDFdoc *pdf);
void cpdf_newpath(CPDFdoc *pdf);


/* Color functions */
void cpdf_setgray(CPDFdoc *pdf, float gray);				/* set both fill and stroke grays */
void cpdf_setrgbcolor(CPDFdoc *pdf, float r, float g, float b);	/* set both fill and stroke colors */
void cpdf_setcmykcolor(CPDFdoc *pdf, float c, float m, float y, float k);
void cpdf_setgrayFill(CPDFdoc *pdf, float gray);
void cpdf_setgrayStroke(CPDFdoc *pdf, float gray);
void cpdf_setrgbcolorFill(CPDFdoc *pdf, float r, float g, float b);
void cpdf_setrgbcolorStroke(CPDFdoc *pdf, float r, float g, float b);
void cpdf_setcmykcolorFill(CPDFdoc *pdf, float c, float m, float y, float k);
void cpdf_setcmykcolorStroke(CPDFdoc *pdf, float c, float m, float y, float k);


/* Graphics state functions */
void cpdf_gsave(CPDFdoc *pdf);
void cpdf_grestore(CPDFdoc *pdf);
void cpdf_setdash(CPDFdoc *pdf, const char *dashspec);
void cpdf_nodash(CPDFdoc *pdf);
void cpdf_concat(CPDFdoc *pdf, float a, float b, float c, float d, float e, float f);
void cpdf_rawConcat(CPDFdoc *pdf, float a, float b, float c, float d, float e, float f);
void cpdf_rotate(CPDFdoc *pdf, float angle);
void cpdf_translate(CPDFdoc *pdf, float xt, float yt);
void cpdf_rawTranslate(CPDFdoc *pdf, float xt, float yt);
void cpdf_scale(CPDFdoc *pdf, float sx, float xy);
void cpdf_setlinewidth(CPDFdoc *pdf, float width);
void cpdf_setflat(CPDFdoc *pdf, int flatness);	/* flatness = 0 .. 100 */
void cpdf_setlinejoin(CPDFdoc *pdf, int linejoin);	/* linejoin = 0(miter), 1(round), 2(bevel) */
void cpdf_setlinecap(CPDFdoc *pdf, int linecap);	/* linecap = 0(butt end), 1(round), 2(projecting square) */
void cpdf_setmiterlimit(CPDFdoc *pdf, float miterlimit);
void cpdf_setstrokeadjust(CPDFdoc *pdf, int flag);		/* PDF-1.2 */

/* Data point Marker funcitons */
void cpdf_marker(CPDFdoc *pdf, float x, float y, int markertype, float size);
void cpdf_pointer(CPDFdoc *pdf, float x, float y, int direction, float size);
void cpdf_errorbar(CPDFdoc *pdf, float x, float y1, float y2, float capsize);
void cpdf_highLowClose(CPDFdoc *pdf, float x, float vhigh, float vlow, float vclose, float ticklen);

/* raw (point-based) versions of marker and other plot symbols above */
void cpdf_rawMarker(CPDFdoc *pdf, float x, float y, int markertype, float size);
void cpdf_rawPointer(CPDFdoc *pdf, float x, float y, int direction, float size);
void cpdf_rawErrorbar(CPDFdoc *pdf, float x, float y1, float y2, float capsize);
void cpdf_rawHighLowClose(CPDFdoc *pdf, float x, float vhigh, float vlow, float vclose, float ticklen);

/* Image related functions */
int cpdf_rawImportImage(CPDFdoc *pdf, const char *imagefile, int type, float x, float y, float angle,
	float *width, float *height, float *xscale, float *yscale, int flags);
int cpdf_importImage(CPDFdoc *pdf, const char *imagefile, int type, float x, float y, float angle,
	float *width, float *height, float *xscale, float *yscale, int flags);
int read_JPEG_header(const char *filename, CPDFimageInfo *jInfo);

int cpdf_placeImageData(CPDFdoc *pdf, const char *uniqueID, const void *imagedata, 
	long length, int nx, int ny, int ncomp_per_pixel, int bits_per_sample,
	float x, float y, float angle, float width, float height, int flags, CPDFimgAttr *imattr);
int cpdf_rawPlaceImageData(CPDFdoc *pdf, const char *uniqueID, const void *imagedata, 
	long length, int nx, int ny, int ncomp_per_pixel, int bits_per_sample,
	float x, float y, float angle, float width, float height, int flags, CPDFimgAttr *imattr);

int cpdf_placeInLineImage(CPDFdoc *pdf, const void *imagedata, int len,
		float x, float y, float angle, float width, float height,
		int pixwidth, int pixheight, int bitspercomp, int CSorMask, int gsave);
int cpdf_rawPlaceInLineImage(CPDFdoc *pdf, const void *imagedata, int len,
		float x, float y, float angle, float width, float height,
		int pixwidth, int pixheight, int bitspercomp, int CSorMask, int gsave);
int cpdf_selectImage(CPDFdoc *pdf, int imgsel);
int cpdf_countImagesInTIFF(const char *imagefile);
int cpdf_readTIFFheader(const char *imagefile, CPDFimageInfo *timgInfo);
int cpdf_readTIFFimageData(char **rasterBuf, CPDFimageInfo *timgInfo, int compress);
int cpdf_readPDFimageHeader(const char *imagefile, CPDFimageInfo *timgInfo);
int cpdf_readPDFimageData(char **rasterBuf, CPDFimageInfo *timgInfo, int ImDataOffset);


/* Memory stream:  most are in cpdfMemBuf.c */
CPDFmemStream *cpdf_setCurrentMemoryStream(CPDFdoc *pdf, CPDFmemStream *memStream); /* cpdfInit.c */

CPDFmemStream *cpdf_openMemoryStream(void);
void cpdf_closeMemoryStream(CPDFmemStream *memStream);
int cpdf_writeMemoryStream(CPDFmemStream *memStream, const char *data, int len);
void cpdf_getMemoryBuffer(CPDFmemStream *memStream, char **streambuf, int *len, int *maxlen);
int cpdf_saveMemoryStreamToFile(CPDFmemStream *stream, const char *name);
int cpdf_memPutc(int ch, CPDFmemStream *memStream);
int cpdf_memPuts(char *str, CPDFmemStream *memStream);
void cpdf_clearMemoryStream(CPDFmemStream *aMstrm);

/* Outline (book mark) functions */
CPDFoutlineEntry *cpdf_addOutlineEntry(CPDFdoc *pdf, CPDFoutlineEntry *afterThis, int sublevel,
	int open, int page, const char *title, int mode, float p1, float p2, float p3, float p4);
CPDFoutlineEntry *cpdf_addOutlineAction(CPDFdoc *pdf, CPDFoutlineEntry *afterThis, int sublevel, int open,
		const char *action_dict, const char *title);

/* Misc functions */
void cpdf_hexStringMode(CPDFdoc *pdf, int flag);
char *cpdf_convertBinaryToHex(const unsigned char *datain, char *hexout, long length, int addFEFF);
unsigned char *cpdf_convertHexToBinary(const char *hexin, unsigned char *binout, long *length);
int cpdf_setDocumentID(CPDFdoc *pdf, int documentID);
struct tm *cpdf_localtime(const time_t *timer, struct tm *p_tm);
time_t cpdf_mktime(struct tm *timeptr);
void  cpdf_setPDFLevel(CPDFdoc *pdf, int major, int minor);
void  cpdf_useStdout(CPDFdoc *pdf, int flag);
char  *cpdf_getOutputFilename(CPDFdoc *pdf);
void cpdf_setOutputFilename(CPDFdoc *pdf, const char *file);
float tm_to_NumDays(struct tm *fromDate, struct tm *toDate);
char  *timestring(int fmt, char *TimeBuf);
int   isLeapYear(int year);
long getFileSize(const char *file);
void  rotate_xyCoordinate(float x, float y, float angle, float *xrot, float *yrot);
float getMantissaExp(float v, int *iexp);
char *cpdf_convertUpathToOS(char *pathbuf, char *upath);
int cpdf_setMonthNames(CPDFdoc *pdf, char *mnArray[]);
int _cpdf_freeMonthNames(CPDFdoc *pdf);



#ifdef MacOS8
/* #include "MacTypes.h" */
void SetFileInfo(const char *fileName, OSType fileType, OSType fileCreator);

#endif


#ifdef MAINDEF
char *cpdf_fontnamelist[] = {
	"Helvetica",
	"Helvetica-Bold",
	"Helvetica-Oblique",
	"Helvetica-BoldOblique",
	"Times-Roman",
	"Times-Bold",
	"Times-Italic",
	"Times-BoldItalic",
	"Courier",
	"Courier-Bold",
	"Courier-Oblique",
	"Courier-BoldOblique",
	"Symbol",
	"ZapfDingbats",

	"AvantGarde-Book",
	"AvantGarde-BookOblique",
	"AvantGarde-Demi",
	"AvantGarde-DemiOblique",
	"Bookman-Demi",
	"Bookman-DemiItalic",
	"Bookman-Light",
	"Bookman-LightItalic",
	"Helvetica-Narrow",
	"Helvetica-Narrow-Oblique",
	"Helvetica-Narrow-Bold",
	"Helvetica-Narrow-BoldOblique",
	"NewCenturySchlbk-Roman",
	"NewCenturySchlbk-Italic",
	"NewCenturySchlbk-Bold",
	"NewCenturySchlbk-BoldItalic",
	"Palatino-Roman",
	"Palatino-Italic",
	"Palatino-Bold",
	"Palatino-BoldItalic",
	"ZapfChancery-MediumItalic",
	
	"Helvetica-Condensed",
	"Helvetica-Condensed-Bold",
	"Helvetica-Condensed-Oblique",
	"Helvetica-Condensed-BoldObl",
	"CPDF-Monospace",
	"CPDF-SmallCap",
/* Japanese */
	"HeiseiKakuGo-W5",
	"HeiseiMin-W3",
/* Korean */
	"HYGoThic-Medium",
	"HYSMyeongJo-Medium",
/* Chinese Traditional */
	"MHei-Medium",
	"MSung-Light",
/* Chinese Simplified */
	"STSong-Light"
};

/* Predefined Encodings and CMap names for Type0 font mappings for Chinese, Japanese, and Korean CID fonts.
   See, page 216, PDF Reference Manual v1.3.
*/
char *cpdf_fontEncodings[] = {
	/* Roman chars */
	"WinAnsiEncoding",
	"MacRomanEncoding",
	"MacExpertEncoding",
	"StandardEncoding",

	/* Generic */
	"Identity-H",
	"Identity-V",

	/* Japanese */
	"83pv-RKSJ-H",
	"90ms-RKSJ-H",		/* MS Shift JIS charset, Win95J */
	"90ms-RKSJ-V",		/* MS Shift JIS, vertical version */
	"90msp-RKSJ-H",		/* MS Shift JIS, with proportional latin chars */
	"90msp-RKSJ-V",		/* MS Shift JIS, vertical version */
	"90pv-RKSJ-H",
	"Add-RKSJ-H",		/* Fujitsu FMR extensions */
	"Add-RKSJ-V",
	"EUC-H",		/* EUC-JP encoding */
	"EUC-V",
	"Ext-RKSJ-H",		/* NEC extenstions */
	"Ext-RKSJ-V",
	"H",			/* ISO-2022-JP encoding */
	"V",
	"UniJIS-UCS2-H",	/* Unicode (UCS-2) encoding for Adobe-Japan1 */
	"UniJIS-UCS2-V",
	"UniJIS-UCS2-HW-H",	/* Latin chars are half-width */
	"UniJIS-UCS2-HW-V",

	/* Chinese */
		/* simplified */
	"GB-EUC-H",
	"GB-EUC-V",
	"GBpc-EUC-H",
	"GBpc-EUC-V",
	"GBK-EUC-H",
	"GBK-EUC-V",
	"UniGB-UCS2-H",		/* Unicode */
	"UniGB-UCS2-V",

		/* traditional */
	"B5pc-H",
	"B5pc-V",
	"ETen-B5-H",
	"ETen-B5-V",
	"ETenms-B5-H",		/* propotional Latin chars */
	"ETenms-B5-V",
	"CNS-EUC-H",
	"CNS-EUC-V",
	"UniCNS-UCS-2-H",	/* Unicode */
	"UniCNS-UCS-2-V",

	/* Korean */
	"KSC-EUC-H",		/* EUC-KR */
	"KSC-EUC-V",
	"KSCms-UHC-H",		/* proportional Latin */
	"KSCms-UHC-V",
	"KSCms-UHC-HW-H",	/* Half width Latin */
	"KSCms-UHC-HW-V",
	"KSCpc-EUC-H",
	"UniKS-UCS2-H",		/* Unicode */
	"UniKS-UCS2-V"
};

#else
	extern char *cpdf_fontnamelist[];
	extern char *cpdf_fontEncodings[];
#endif


#ifdef __cplusplus
}
#endif /* __cplusplus */


#endif  /*  __CLIBPDF_H__  */

