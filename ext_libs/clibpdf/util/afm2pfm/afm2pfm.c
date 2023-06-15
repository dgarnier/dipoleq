/*  afm2pfm.c
 *  afm2pfm - convert Adobe .afm file to MS-Windows binary .pfm file
 *  Copyright (C) 1994  Russell Lang <rjl@eng.monash.edu.au>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

/*
 *  Function:
 *    Convert an AFM (Adobe Font Metric) file to a MS-Windows PFM
 *    (Printer Font Metric) file.  This is to allow Type 1 fonts
 *    with AFM files to be used with MS-Window Adobe Type Manager.
 * 
 *  Syntax:
 *    afm2pfm  [-w] [-?] infile [outfile]
 *
 *  Notes:
 *    1. This program was written to produce a PFM file for a simple
 *       Type 1 font that I had created.  More complex AFM files may 
 *       expose limitations in this program.  Warning messages are
 *       produced for AFM lines or sections that are ignored.
 *       Ligatures are silently ignored.  Track kerning and Composites
 *       are ignored.  Unencoded characters are ignored.
 *    2. This has been tested using the Borland C++ 3.1
 *       compiler and the EMX/GCC 0.8h compiler.
 *       For other compilers, be careful of the packing of the
 *       PFM structures.
 */

/* afm2pfm
 * Version 1.0   1994-01-06  rjl
 * Version 1.1   1995-03-10  rjl
 *   fixes provided by Norman Walsh
 */

/*
 *  The following program was very useful in developing afm2pfm:
 *    Title:  pfm2afm - Convert Windows .pfm files to .afm files
 *    Author: Ken Borgendale   10/9/91  Version 1.0
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "afm2pfm.h"

#define BUFSIZE 4096
/*
 *  Function Prototypes
 */
void  help (void);
void  parseargs (int argc, char ** argv);
void  vcheck(char *key, char *value);
void  lcheck(void);
void  openafm(void);
void  initpfm(void);
void  makeheader(void);
void  makeextent(void);
void  makeaverage(void);
void  makekern(void);
void  makecomposites(void);
void  maketrailer(void);
void  writepfm(void);

/*
 *  Global variables
 */
FILE   * inf;                /* Input file */
FILE   * outf;               /* Output file */
char   infname[272];         /* Input file name */
char   outfname[272];        /* Output file name */

PFMHEADER	pfm;		/* .pfm header */
PFMEXTENSION	pfmext;
char		*devicename = "PostScript";
char		facename[64];
char		fontname[64];
EXTTEXTMETRIC	etm;		/* metrics extension */
WORD		extent[256];	/* extent table */
int		cextent;	/* count of entries in extent table */
KERNPAIR	kernpair[512];	/* kern pair table */
short		ckernpair;	/* count of kern pairs */
KERNTRACK	kerntrack[16];	/* kern track table */
short		ckerntrack;	/* count of kern track entries */

long 		sumwidth;	/* used for calculating average width */
int		maxwidth;
BBOX		fbbox;		/* font bounding box */
BBOX		cbbox[256];	/* bounding boxes for characters */
int		ascender;	/* from AFM Ascender */
int		descender;	/* from AFM Descender */

char   line[256];            /* input line from afm */
int    linecount;	     /* number of lines read */
int    charcount;            /* number of characters */

char   * name[256];          /* pointer to names of characters */
char   * nbuffer;            /* buffer to hold character names */
char   * ntail;		     /* tail of name buffer */

BOOL   warnflag = TRUE;      /* suppress warnings flag */

/*
 *  Do the function
 */
int main(int argc, char **argv) {

    /* Parse arguments */
    parseargs(argc, argv);

    /* Open the afm file */
    openafm();

    /* Create pfm in memory */
    initpfm();
    makeheader();
    makeextent();
    makekern();
    makecomposites();
    maketrailer();

    /* write out pfm from memory to file */
    writepfm();

    if (nbuffer)
	free(nbuffer);
    fclose(inf);
    fclose(outf);
    return 0;
}

/*
 *  Put out normal help
 */
void  help (void) {
    puts("\nafm2pfm - Convert afm to pfm - Version 1.1   1994-01-06\n");
    puts("This utility converts afm files for Adobe type 1 fonts");
    puts("to binary pfm files for use on MS-Windows.\n");
    puts("afm2pfm [-w] [-?] infile  [outfile]");
    puts("    -w suppresses warnings");
    puts("    The extension .afm is added to the infile if it has none.");
    puts("    The outfile is defaulted from the input file name.");
    puts("Russell Lang - rjl@monu1.cc.monash.edu.au\n");
    exit (1);
}


/*
 *  Parse arguments.
 */
void  parseargs (int argc, char ** argv) {
    char  swchar;
    int    argcnt;
    int    filecnt;
    char * argp;

    if (argc == 1)
        help();

    argp = argv[1];
    if (*argp == '-' || *argp == '/') {
	switch(argp[1]) {
	    case '?':
		help();
	    case 'w':
		warnflag = !warnflag;
		break;
	}
	argv++;
	argc--;
    }

    
    if (argc == 3)  {
        strcpy(infname, argv[1]);
        strcpy(outfname, argv[2]);
    }
    else if (argc == 2) {
        strcpy(infname, argv[1]);
    }
    else
	help();
}

/*
 * Initialise the .pfm structures
 */
void initpfm(void) {
    /* storage for names */
    nbuffer = malloc(BUFSIZE);
    memset(nbuffer, 0, BUFSIZE);
    ntail = nbuffer;

    pfm.dfVersion = 256;
    pfm.dfType = 0x81;	/* PF_DEVICE_REALIZED | PF_VECTOR_TYPE */
    pfm.dfPoints = 10; 
    pfm.dfVertRes = 300; 
    pfm.dfHorizRes = 300;

    pfmext.dfSizeFields = sizeof(PFMEXTENSION);

    etm.etmSize = sizeof(EXTTEXTMETRIC);
    /* set some reasonable defaults - taken from OS/2 helv.pfm */
    etm.etmPointSize = 240;	/* units? */
    etm.etmMasterHeight = 1000;	/* or 300? */
    etm.etmMinScale = 3;
    etm.etmMaxScale = 1000;
    etm.etmMasterUnits = 1000;
    etm.etmCapHeight = 729;
    etm.etmXHeight = 524;
    etm.etmLowerCaseAscent = 729;
    etm.etmUpperCaseDescent = 218;
    etm.etmSuperScript = -500;
    etm.etmSubScript = 250;
    etm.etmSuperScriptSize = 500;
    etm.etmSubScriptSize = 500;
    etm.etmUnderlineOffset = 100;
    etm.etmUnderlineWidth = 50;
    etm.etmDoubleUpperUnderlineOffset = 50;
    etm.etmDoubleLowerUnderlineOffset = 100;
    etm.etmDoubleUpperUnderlineWidth = 25;
    etm.etmDoubleLowerUnderlineWidth = 25;
    etm.etmStrikeOutOffset = 339;
    etm.etmStrikeOutWidth = 50;
}
	
void lcheck(void) {
    if (strlen(line) >= sizeof(line)-1) {
        fprintf(stderr,"Error on line %d, line too long\n", linecount);
        exit(2);
    }
}

void vcheck(char *key, char *value) {
    if (value != (char *)NULL) 
	return;
    fprintf(stderr,"Error on line %d, key \042%s\042 requires a value\n",
    	linecount,key);
    exit(2);
}

void makeheader(void) {
    BOOL inhdr = TRUE;
    BOOL isafm = FALSE;
    BOOL gotcopyright = FALSE;
    char *key, *value;
    key = line;
    while ( inhdr && (fgets(line, sizeof(line), inf)!=(char *)NULL) ) {
    	linecount++;
	lcheck();
	value = strchr(line, '\n');
	if (value)
	    *value = '\0';	/* remove trailine \n */
	value = strchr(line, ' ');
	if (value)
	    *value++ = '\0';	/* separate key from value */
	if (strcmp(key, "StartFontMetrics")==0) {
	    isafm = TRUE;
	    vcheck(key, value);
	    continue;
	}
	if (!isafm) {
	    fputs("Error: Not an AFM file\n", stderr);
	    exit(2);
	}
	if (strcmp(key, "FontName")==0) {
	    vcheck(key, value);
	    strcpy(fontname, value);
	    if (strstr(fontname, "Italic")!=(char *)NULL)
		pfm.dfItalic = 1;
	}
	else if (strcmp(key, "FamilyName")==0) {
	    vcheck(key, value);
	    strcpy(facename, value);
	    /* should scan name and guess dfPitchAndFamily field */
	}
	else if (strcmp(key, "Notice")==0) {
	    vcheck(key, value);
	    strncpy(pfm.dfCopyright, value, sizeof(pfm.dfCopyright)-1);
	}
	else if (strcmp(key, "EncodingScheme")==0) {
		vcheck(key, value);
		if (strcmp(value, "AdobeStandardEncoding")!=0)
		    pfm.dfCharSet = 2;   /* use Symbol */
	}
	else if (strcmp(key, "Weight")==0) {
	    vcheck(key, value);
	    if (strcmp(value,"Bold")==0)
	    	pfm.dfWeight = 700;
	    else if (strcmp(value,"Light")==0)
	        pfm.dfWeight = 300;
	    else
	        pfm.dfWeight = 400;
	}
	else if (strcmp(key, "IsFixedPitch")==0) {
	    vcheck(key, value);
	    if (strcmp(value,"false")==0) 
	        pfm.dfPitchAndFamily |= 1;	/* variable width */
	}
	else if (strcmp(key, "FontBBox")==0) {
	    vcheck(key, value);
	    if (sscanf(value,"%d %d %d %d", &fbbox.llx, &fbbox.lly, &fbbox.urx, &fbbox.ury)!=4) {
	    	fprintf(stderr,"Error on line %d, missing value\n",linecount);
		exit(2);
	    }
	}
	else if (strcmp(key, "CapHeight")==0) {
	    vcheck(key, value);
	    etm.etmCapHeight = atoi(value);
	}
	else if (strcmp(key, "XHeight")==0) {
	    vcheck(key, value);
	    etm.etmXHeight = atoi(value);
	}
	else if (strcmp(key, "Ascender")==0) {
	    vcheck(key, value);
	    ascender = atoi(value);
	}
	else if (strcmp(key, "Descender")==0) {
	    vcheck(key, value);
	    descender = atoi(value);
	}
	else if (strcmp(key, "ItalicAngle")==0) {
	    vcheck(key, value);
	    etm.etmSlant = (short)(atof(value)*10);
	    if (etm.etmSlant)
		pfm.dfItalic = 1;
	}
	else if (strcmp(key, "UnderlinePosition")==0) {
	    vcheck(key, value);
	    etm.etmUnderlineOffset = -atoi(value);
	}
	else if (strcmp(key, "UnderlineThickness")==0) {
	    vcheck(key, value);
	    etm.etmUnderlineWidth = atoi(value);
	}
	else if (strcmp(key, "StartCharMetrics")==0) {
	    vcheck(key, value);
	    charcount = atoi(value);
	    inhdr = FALSE;
	}
	else if (strcmp(key, "FullName")==0) {
	    /* nop */
	}
	else if (strcmp(key, "Comment")==0) {
	    /* nop */
	}
	else if (strcmp(key, "Version")==0) {
	    /* nop */
	}
    	else {
	    if (warnflag) {
		fprintf(stderr,"Warning on line %d: Unknown key: %s ",linecount,key);
    	        if (value)
    	           fputs(value,stderr);
	        fputc('\n',stderr);
	    }
    	}
    }
    /* time for StartCharMetric */
}

void makeextent(void) {
    int c, len;
    BOOL inmetric = TRUE;
    char *key, *value, *p;
    pfm.dfFirstChar = 255;
    pfm.dfLastChar = 0;
    if (strncmp(line, "StartCharMetrics",16)!=0) {
        fprintf(stderr,"Error on line %d, expecting StartCharMetrics\n",linecount);
        exit(2);
    }
    while ( inmetric && (fgets(line, sizeof(line), inf)!=(char *)NULL) ) {
    	linecount++;
	lcheck();
    	c = 0;
        key = line;
	value = strchr(line, '\n');
	if (value)
	    *value = '\0';	/* remove trailine \n */
        if (strcmp(key, "EndCharMetrics")==0) {
            pfm.dfAvgWidth = sumwidth / (pfm.dfLastChar - pfm.dfFirstChar + 1); /* wrong - should be weighted */
            pfm.dfMaxWidth = maxwidth;
            inmetric = FALSE;
        }
	key = strtok(line, ";");
	while (inmetric && key) {
	    while (*key == ' ')
		key++;
	    value = strchr(key,' ');
	    if (value)
	        value++;
	    while (value && (*value == ' '))
		value++;
	    if (strncmp(key, "C ", 2)==0) {
	        vcheck(key, value);
	        c = atoi(value);
	        if ((c>0) && (c<pfm.dfFirstChar))
	            pfm.dfFirstChar = c;
	        if ((c>0) && (c>pfm.dfLastChar))
	            pfm.dfLastChar = c;
	        if ( (c<0) && warnflag )
	    	    fprintf(stderr,"Warning on line %d, ignoring %s\n",linecount,key);
	    }
	    else if (strncmp(key, "WX ",3)==0) {
	        vcheck(key, value);
	        if (c>=0) {
	            extent[c] = atoi(value);
		    sumwidth += extent[c];
		    if (extent[c] > maxwidth)
		        maxwidth = extent[c];
	        }
	    }
	    else if (strncmp(key, "N ",2)==0) {
		vcheck(key, value);
	        if (c>=0) {
		    p = strchr(value, ' ');
		    if (p)
		        len = p - value;
	            else
	                len = strlen(value);
	            if ( (len+1) < BUFSIZE-(ntail-nbuffer) ) {
	                name[c] = ntail;
	                strncpy(ntail, value, len);
	                ntail += len;
	                *ntail++ = '\0';
	            }
	            else {
	                fprintf(stderr,"Error on line %d: out of memory for names\n",
	                    linecount);
	                exit(2);
	            }
		}
	    }
	    else if (strncmp(key, "B ",2)==0) {
		vcheck(key, value);
	        if (c>=0) {
	            if (sscanf(value,"%d %d %d %d", &cbbox[c].llx, &cbbox[c].lly, 
				&cbbox[c].urx, &cbbox[c].ury)!=4) {
	    	        fprintf(stderr,"Error on line %d, missing BoundingBox value\n",linecount);
		        exit(2);
	            }
		}
	    }
	    else if (strncmp(key, "L ",2)==0) {
	    	/* ignore it */
	    }
	    else {
		if (warnflag)
	    	    fprintf(stderr,"Warning on line %d, ignoring %s\n",linecount,key);
	    }
	    key = strtok(NULL, ";");
	}
    }
    fgets(line, sizeof(line), inf);
    linecount++;
    lcheck();
    if (pfm.dfFirstChar > pfm.dfLastChar)
	pfm.dfFirstChar = pfm.dfLastChar;
    /* now fix up some dimensions */
    if (pfm.dfMaxWidth == 0)
	pfm.dfMaxWidth = fbbox.urx - fbbox.llx;	/* guess from font bbox */
    pfm.dfAscent = fbbox.ury;
    if ((pfm.dfCharSet==0) && name['d'] && *name['d'] && (strcmp(name['d'], "d")==0) && (cbbox['d'].ury))
        etm.etmLowerCaseAscent = cbbox['d'].ury;
    else {
	if (ascender) {
	    if (warnflag)
	        fprintf(stderr,"Warning: setting lower case ascent from AFM Ascender\n");
            etm.etmLowerCaseAscent = ascender;
	}
	else {
	    if (warnflag)
	        fprintf(stderr,"Warning: setting lower case ascent from font bounding box\n");
            etm.etmLowerCaseAscent = fbbox.ury;
        }
    }
    if ((pfm.dfCharSet==0) && name['p'] && *name['p'] && (strcmp(name['p'], "p")==0) && (cbbox['p'].lly))
        etm.etmUpperCaseDescent = cbbox['p'].lly;
    else {
	if (descender) {
	    if (warnflag)
	        fprintf(stderr,"Warning: setting upper case descent from AFM Descender\n");
            etm.etmUpperCaseDescent = descender;
	}
	else {
	    if (warnflag)
	        fprintf(stderr,"Warning: setting upper case descent from font bounding box\n");
            etm.etmUpperCaseDescent = fbbox.lly;
        }
    }
    if ((etm.etmXHeight==0) && (pfm.dfCharSet==0) && *name['x'] && 
        (strcmp(name['x'], "x")==0) && (cbbox['x'].ury))
        etm.etmXHeight = cbbox['x'].ury;
    if ((etm.etmCapHeight==0) && (pfm.dfCharSet==0) && *name['H'] && 
        (strcmp(name['H'], "H")==0) && (cbbox['H'].ury))
        etm.etmCapHeight = cbbox['H'].ury;
    pfm.dfPixHeight = fbbox.ury - fbbox.lly;
    pfm.dfPixWidth  = fbbox.urx - fbbox.llx;
    makeaverage();
}

void makeaverage(void) {
    /* Formula from OS/2 PM documentation for FONTMETRICS */
    if ((pfm.dfCharSet==0) && (pfm.dfFirstChar<=' ') && (pfm.dfLastChar>='z')) {
	pfm.dfAvgWidth = (
	  extent['a']*64L  + extent['b']*14L + extent['c']*27L + extent['d']*35L +
	  extent['e']*100L + extent['f']*20L + extent['g']*14L + extent['h']*42L +
	  extent['i']*63L  + extent['j']* 3L + extent['k']* 6L + extent['l']*35L +
	  extent['m']*20L  + extent['n']*56L + extent['o']*56L + extent['p']*17L +
	  extent['q']* 4L  + extent['r']*49L + extent['s']*56L + extent['t']*71L +
	  extent['u']*31L  + extent['v']*10L + extent['w']*18L + extent['x']* 3L +
	  extent['y']*18L  + extent['z']* 2L + extent[' ']*166L
	  ) / 1000;
	return;
    }
    if (warnflag)
	fprintf(stderr,"Warning: calculating average width from all characters\n");
}


/* untested */
/* ignores kern pairs for unencoded characters */
/* ignores kern track info */
void makekern(void) {
    BOOL inkern = TRUE;
    BOOL inpair = FALSE;
    BOOL intrack = FALSE;
    char *key, *value, *p, *pair1, *pair2;
    int i, c1, c2, kerncount, kerncount2;
    key = line;
    value = strchr(line, '\n');
    if (value)
        *value = '\0';	/* remove trailine \n */
    if (strncmp(key, "StartKernData",13)!=0)
	return;
    while ( inkern && (fgets(line, sizeof(line), inf)!=(char *)NULL) ) {
    	linecount++;
	lcheck();
	value = strchr(line, '\n');
	if (value)
	    *value = '\0';	/* remove trailine \n */
        if (strcmp(key, "EndKernData")==0) {
            inkern = FALSE;
            if (inpair && warnflag)
		fprintf(stderr,"Warning on line %d, expecting EndKernPairs\n",linecount);
            continue;
        }
        else if (strncmp(key, "StartKernTrack", 14)==0) {
	    intrack = TRUE;
	    if (warnflag)
	        fprintf(stderr,"Warning on line %d, ignoring StartKernTrack\n",linecount);
	}
        else if (strcmp(key, "EndKernTrack")==0) {
	    intrack = FALSE;
	}
        else if (strncmp(key, "StartKernPairs", 14)==0) {
            inpair = TRUE;
	    value = strchr(line, ' ');
	    vcheck(key, value);
	    kerncount = atoi(value);
	    kerncount2 = 0;
	    continue;
        }
        else if (strcmp(key, "EndKernPairs")==0) {
            inpair = FALSE;
            etm.etmKernPairs = ckernpair;
            continue;
        }
        else if (inpair) {
	    key = strtok(line, " ");
	    if (!key)
	        continue;
	    if (strcmp(key, "KPX")==0) {
	        pair1 = strtok(NULL," ");
	        vcheck(key, pair1);
	        pair2 = strtok(NULL," ");
	        vcheck(key, pair2);
	        value = strtok(NULL," ");
	        vcheck(key, value);
	        c1 = c2 = -1;
	        for (i = 0; i < 256; i++)
	            if (name[i] && (strcmp(name[i],pair1)==0)) {
	                c1 = i;
	                break;
	            }
	        for (i = 0; i < 256; i++)
	            if (name[i] && (strcmp(name[i],pair2)==0)) {
	                c2 = i;
	                break;
	            }
		kerncount2++;
		if (kerncount2 > kerncount) {
		    if (warnflag)
	    	        fprintf(stderr,"Warning on line %d, too many kern pairs, ignoring extras\n",linecount);
		}
		else if (ckernpair > sizeof(kernpair)/sizeof(KERNPAIR)) {
		    if (warnflag)
	    	        fprintf(stderr,"Warning on line %d, too many kern pairs, ignoring extras\n",linecount);
		}
		else if (c1==-1) {
		    if (warnflag)
	    	        fprintf(stderr,"Warning on line %d, ignoring kerning for unencoded character: %s\n",linecount, pair1);
		}
		else if (c2==-1) {
		    if (warnflag)
	    	        fprintf(stderr,"Warning on line %d, ignoring kerning for unencoded character: %s\n",linecount, pair2);
		}
		else {
		    kernpair[ckernpair].kpPair.each[0] = c1;
		    kernpair[ckernpair].kpPair.each[1] = c2;
		    kernpair[ckernpair].kpKernAmount = atoi(value);
		    ckernpair++;
		}
	    }
	    else {
		if (warnflag)
	    	    fprintf(stderr,"Warning on line %d, ignoring %s\n",linecount,key);
	    }
	    if (pair1)
	        value++;
	    while (value && (*value == ' '))
		value++;
        }
	else if (intrack) {
	    /* ignore it */
	}
	else {
	    if (warnflag)
    	        fprintf(stderr,"Warning on line %d, ignoring %s\n",linecount,key);
	}
    }
    fgets(line, sizeof(line), inf);
    linecount++;
    lcheck();
}

/* currently just ignores composites */
void makecomposites(void) {
    BOOL incomp = TRUE;
    char *key, *value;
    key = line;
    if (strncmp(key, "StartComposites",15)!=0)
	return;
    if (warnflag)
        fprintf(stderr,"Warning on line %d, ignoring %s\n",linecount,key);
    while ( incomp && (fgets(line, sizeof(line), inf)!=(char *)NULL) ) {
    	linecount++;
	lcheck();
	value = strchr(line, '\n');
	if (value)
	    *value = '\0';	/* remove trailine \n */
        if (strcmp(key, "EndComposites")==0) {
            incomp = FALSE;
            continue;
        }
    }
    fgets(line, sizeof(line), inf);
    linecount++;
    lcheck();
}

/* currently just ignores trailer */
void maketrailer(void) {
    char *key, *value;
    key = line;
    if (strncmp(key, "EndFontMetrics",14)==0)
    	    return;
    if (warnflag)
        fprintf(stderr,"Warning on line %d, expecting EndFontMetrics\n",linecount);
}


/*
 *  Create the .pfm file
 */
void  writepfm(void) {
    char  * cp;
    int len;

    /* Add .pfm if there is none */
    if (!*outfname) {
        strcpy(outfname, infname);
        cp = outfname+strlen(outfname)-1;
        while (cp >= outfname && *cp!='.' && *cp!='\\' && *cp!='/' && *cp!=':')
           cp--;
        if (*cp=='.') *cp=0;
        strcat(outfname, ".pfm");
    }
    /* Open the file */
    outf = fopen(outfname, "wb");
    if (!outf) {
        fputs("Error: Unable to open output file - ", stderr);
        fputs(outfname, stderr);
        fputc('\n', stderr);
        exit(3);
    }
    len = sizeof(pfm) + sizeof(pfmext);
    if (*devicename) {
        pfm.dfDevice = len;
        len += strlen(devicename)+1;
    }
    if (*facename) {
        pfm.dfFace = len;
        len += strlen(facename)+1;
    }
    if (*fontname) {
        pfmext.dfDriverInfo = len;
        len += strlen(fontname)+1;
    }
    pfmext.dfExtMetricsOffset = len;
    len += sizeof(etm);
    if (pfm.dfLastChar != pfm.dfFirstChar) {
        pfmext.dfExtentTable = len;
        len += sizeof(WORD) * (pfm.dfLastChar - pfm.dfFirstChar + 1);
    }
    if (ckernpair) {
        pfmext.dfPairKernTable = len;
	len += sizeof(ckernpair);
        len += ckernpair * sizeof(KERNPAIR);
    }
    if (ckerntrack) {
	/* untested */
        pfmext.dfTrackKernTable = len;
        len += ckerntrack * sizeof(KERNTRACK);
    }
    pfm.dfSize = len;
    fwrite(&pfm, 1, sizeof(pfm), outf);
    fwrite(&pfmext, 1, sizeof(pfmext), outf);
    if (*devicename)
        fwrite(devicename, 1, strlen(devicename)+1, outf);
    if (*facename)
        fwrite(facename, 1, strlen(facename)+1, outf);
    if (*fontname)
        fwrite(fontname, 1, strlen(fontname)+1, outf);
    fwrite(&etm, 1, sizeof(etm), outf);
    if (pfm.dfLastChar != pfm.dfFirstChar)
        fwrite(extent+pfm.dfFirstChar, sizeof(WORD), 
           pfm.dfLastChar - pfm.dfFirstChar + 1, outf);
    if (ckernpair) {
	fwrite(&ckernpair, sizeof(ckernpair), 1, outf);
        fwrite(kernpair, sizeof(KERNPAIR), ckernpair, outf);
    }
    if (ckerntrack)
        fwrite(kerntrack, sizeof(KERNTRACK), ckerntrack, outf);
    fclose(outf);
}


/*
 *  Open the .afm file
 */
void  openafm(void) {
    char  * cp;

    /* Check for a file extension */
    cp = infname+strlen(infname)-1;
    while (cp>=infname && *cp!='.' && *cp!='\\' && *cp!='/' && *cp!=':')
       cp--;
    if (*cp!='.')
        strcat(infname, ".afm");

    /* Open the file */
    inf = fopen(infname, "r");
    if (!inf) {
        fputs("Error: Unable to open input file - ", stderr);
        fputs(infname, stderr);
        fputc('\n', stderr);
        exit(4);
    }
}

