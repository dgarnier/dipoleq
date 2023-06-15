/* weatherj.c -- plots NOAA Weather data - Japanese version
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

cc -Wall -o weatherj -I/usr/local/include weatherj.c -lcpdfpm

1999-09-26 [io] Japanese version.
1999-08-22 [io] for ver 2.00
1998       [io] original version

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
/* #include <time.h> */

#include <cpdflib.h>

#define	CITYNAME	"OAKLAND"
#define MAXDATA		1500

char *tempUnit[] = { "度 F", "度 C", "" };

/* Some plot data and globals for settings */
/* --------------------------------------------------------------------------- */
static int   Ndata = 0;			/* # data points */
static float temp[MAXDATA];		/* temp F */
static float humidity[MAXDATA];
static struct tm xtime[MAXDATA];
static struct tm axStart, axEnd, recStart, recEnd;
static CPDFaxis *tAxis;
static char sky[128];

/* --------------------------------------------------------------------------- */


void plot_temp(CPDFdoc *pdf, int useC);
void plot_rHumidity(CPDFdoc *pdf);
void do_TextStuff(CPDFdoc *pdf, int useC);
int read_datafile(char *gafile, int useC);

void plot_temp(CPDFdoc *pdf, int useC)
{
int i;
float vtmax= -999.0, vtmin = 9999.0;
float recmin = 50.0, recmax = 100.0;
float xt;
CPDFplotDomain *tempDomain, *oldDomain;
CPDFaxis *yAxis;

	/* get data max and min values */
	for(i=0; i<Ndata; i++) {
	    if(temp[i] > vtmax) vtmax = temp[i];
	    if(temp[i] < vtmin) vtmin = temp[i];
	}

	cpdf_suggestMinMaxForLinearDomain(vtmin, vtmax, &recmin, &recmax);
	/* recmin = vtmin; */
	/* recmax = vtmax; */
	/* fprintf(stderr, "(%g, %g) -->rec axis min,max = (%g, %g)\n", vtmax, vtmin, recmin, recmax); */

	memcpy(&axStart, &xtime[0], sizeof(struct tm));
	axStart.tm_sec = 0;
	axStart.tm_min = 0;
	axStart.tm_hour = 0;
    	mktime(&axStart);		/* this rebuilds all the fields */
	memcpy(&axEnd, &xtime[Ndata-1], sizeof(struct tm));
	axEnd.tm_sec = 0;
	axEnd.tm_min = 0;
	axEnd.tm_hour = 0;
	axEnd.tm_mday++;
    	mktime(&axEnd);			/* this rebuilds all the fields */
        cpdf_suggestTimeDomainParams(&axStart, &axEnd, &recStart, &recEnd);
    	tempDomain = cpdf_createTimePlotDomain( 1.7*inch, 5.5*inch, 6.0*inch, 2.5*inch,
			    &recStart, &recEnd, recmin, recmax, TIME, LINEAR, 0);
	oldDomain = cpdf_setPlotDomain(pdf, tempDomain);	/* save oldDomain for later restore */
	if(useC)
	    cpdf_setLinearMeshParams(tempDomain, Y_MESH, recmin, 5.0, recmin, 1.0);
	else
	    cpdf_setLinearMeshParams(tempDomain, Y_MESH, recmin, 10.0, recmin, 5.0);
	cpdf_fillDomainWithRGBcolor(tempDomain, 1.0, 0.9, 0.9);  /* color of domain: light pink */
	cpdf_drawMeshForDomain(tempDomain);

	cpdf_setgray(pdf, 0.0);

	/* X-Axis --------------------------------------------------------------------------------- */
	tAxis = cpdf_createTimeAxis(0.0, 6.0*inch, TIME, &recStart, &recEnd);
	cpdf_attachAxisToDomain(tAxis, tempDomain, 0.0, -0.2*inch );
	cpdf_setTimeAxisNumberFormat(tAxis, MONTH_NAME, YEAR_2DIGIT, "Helvetica", 12.0);	
		/* YEAR_FULL for 1998 or YEAR_2DIGIT for 98 */
	cpdf_setAxisLabel(tAxis, "日付", "HeiseiKakuGo-W5", "90ms-RKSJ-H", 16.0);
	cpdf_drawAxis(tAxis);
	/* cpdf_freeAxis(tAxis); */	/* we use it again for plotting other things */

	/* Y-Axis --------------------------------------------------------------------------------- */
	yAxis = cpdf_createAxis(90.0, 2.5*inch, LINEAR, recmin, recmax);
	cpdf_attachAxisToDomain(yAxis, tempDomain, -0.2*inch, 0.0);	/* give a little X offset left */
	cpdf_setAxisNumberFormat(yAxis, "%.f", "Helvetica", 12.0);
	if(useC)
	    cpdf_setAxisLabel(yAxis, "気温 [度 C]", "HeiseiKakuGo-W5", "90ms-RKSJ-H", 16.0);
	else
	    cpdf_setAxisLabel(yAxis, "気温 [度 F]", "HeiseiKakuGo-W5", "90ms-RKSJ-H", 16.0);

	if(useC)
	    cpdf_setLinearAxisParams(yAxis, recmin, 5.0, recmin, 1.0);
	else	
	    cpdf_setLinearAxisParams(yAxis, recmin, 10.0, recmin, 5.0);
	cpdf_drawAxis(yAxis);
	cpdf_freeAxis(yAxis);

	/* Do the tuning curve first */
	cpdf_comments(pdf, "\n%% plot the temperature curve.\n");
	cpdf_setlinewidth(pdf, 1.5);
	cpdf_setrgbcolorStroke(pdf, 1.0, 0.0, 0.0);	/* red */

	for(i=0; i<Ndata; i++) {
	    xt = tm_to_NumDays(&recStart, &xtime[i]);
	    if(i==0) cpdf_moveto(pdf, xt, temp[i]);
	    else     cpdf_lineto(pdf, xt, temp[i]);
	    /* printf("%d-%d-%d %d:%d - xt=%g, temp= %.1f\n", xtime[i].tm_year, xtime[i].tm_mon+1, xtime[i].tm_mday,
	    			xtime[i].tm_hour, xtime[i].tm_min, xt, temp[i]);
	    */
	}
	cpdf_setlinejoin(pdf, 1);
	cpdf_stroke(pdf);
	cpdf_setlinejoin(pdf, 0);

	cpdf_setPlotDomain(pdf, oldDomain);		/* restore previous plot domain */
	cpdf_freePlotDomain(tempDomain);		/* deallocate the plot domain */
}


void plot_rHumidity(CPDFdoc *pdf)
{
int i;
float vtmax= -999.0, vtmin = 9999.0;
float recmin = 0.0, recmax = 100.0;
float xt;
CPDFplotDomain *rhDomain, *oldDomain;
CPDFaxis *yAxis;

	/* get data max and min values */
	for(i=0; i<Ndata; i++) {
	    if(humidity[i] > vtmax) vtmax = humidity[i];
	    if(humidity[i] < vtmin) vtmin = humidity[i];
	}

	/* recmin = vtmin; */
	/* recmax = vtmax; */
	cpdf_suggestMinMaxForLinearDomain(vtmin, vtmax, &recmin, &recmax);
	/* fprintf(stderr, "rec axis min=%g, max=%g\n", recmin, recmax); */

        cpdf_suggestTimeDomainParams(&axStart, &axEnd, &recStart, &recEnd);
    	rhDomain = cpdf_createTimePlotDomain( 1.7*inch, 1.7*inch, 6.0*inch, 2.5*inch,
			    &recStart, &recEnd, recmin, recmax, TIME, LINEAR, 0);
	oldDomain = cpdf_setPlotDomain(pdf, rhDomain);	/* save oldDomain for later restore */

	cpdf_setLinearMeshParams(rhDomain, Y_MESH, recmin, 10.0, recmin, 5.0);
	cpdf_fillDomainWithRGBcolor(rhDomain, 1.0, 0.9, 0.9);  /* color of domain: light pink */
	cpdf_drawMeshForDomain(rhDomain);

	cpdf_setgray(pdf, 0.0);

	/* X-Axis - we already have the time axis for temp plot ----------------------------------- */
	/* tAxis = cpdf_createTimeAxis(0.0, 6.0*inch, TIME, &recStart, &recEnd); */
	cpdf_attachAxisToDomain(tAxis, rhDomain, 0.0, -0.2*inch );
	/* cpdf_setTimeAxisNumberFormat(tAxis, MONTH_NAME, YEAR_2DIGIT, "Helvetica", 12.0); */
		/* YEAR_FULL for 1998 or YEAR_2DIGIT for 98 */
	/* cpdf_setAxisLabel(tAxis, "日付", "HeiseiKakuGo-W5", "90ms-RKSJ-H", 16.0); */
	cpdf_drawAxis(tAxis);
	cpdf_freeAxis(tAxis);	/* we are done, no more use */

	/* Y-Axis --------------------------------------------------------------------------------- */
	yAxis = cpdf_createAxis(90.0, 2.5*inch, LINEAR, recmin, recmax);
	cpdf_attachAxisToDomain(yAxis, rhDomain, -0.2*inch, 0.0);	/* give a little X offset left */
	cpdf_setAxisNumberFormat(yAxis, "%.f", "Helvetica", 12.0);
	cpdf_setAxisLabel(yAxis, "湿度 [%]", "HeiseiKakuGo-W5", "90ms-RKSJ-H", 16.0);

	cpdf_setLinearAxisParams(yAxis, recmin, 10.0, recmin, 5.0);
	cpdf_drawAxis(yAxis);
	cpdf_freeAxis(yAxis);

	/* Do the tuning curve first */
	cpdf_comments(pdf, "\n%% plot the humidity curve.\n");
	cpdf_setlinewidth(pdf, 1.5);
	cpdf_setrgbcolorStroke(pdf, 0.0, 0.0, 1.0);	/* blue */

	for(i=0; i<Ndata; i++) {
	    xt = tm_to_NumDays(&recStart, &xtime[i]);
	    if(i==0) cpdf_moveto(pdf, xt, humidity[i]);
	    else     cpdf_lineto(pdf, xt, humidity[i]);
	    /* printf("%d-%d-%d %d:%d - xt=%g, temp= %.1f\n", xtime[i].tm_year, xtime[i].tm_mon+1, xtime[i].tm_mday,
	    			xtime[i].tm_hour, xtime[i].tm_min, xt, humidity[i]);
	    */
	}
	cpdf_setlinejoin(pdf, 1);
	cpdf_stroke(pdf);
	cpdf_setlinejoin(pdf, 0);

	cpdf_setPlotDomain(pdf, oldDomain);		/* restore previous plot domain */
	cpdf_freePlotDomain(rhDomain);		/* deallocate the plot domain */
}


/* -------------------------------------------------------------------------------- */
/* Key words for data file entries */

int read_datafile(char *gafile, int useC)
{
float vt, vdew, vrh;
int year, month, day, hour, minute, second;
char linebuf[1024];		/* input line buffer */
char strbuf1[128];
FILE *fp;
	if((fp = fopen(gafile, "r")) == NULL) {
	    fprintf(stderr, "Can't open data file: %s\n", gafile);
	    exit(1);
	}

	Ndata = 0;
	while( ( fgets(linebuf, 1020, fp) ) != NULL ) {
	    if(linebuf[0] == '#' || linebuf[0] == '%' || strlen(linebuf) < 3)
		continue;		/* pretty much blank line */
	    sscanf( linebuf, "%s", strbuf1);	/* first word */
	    if(strstr(strbuf1, "19") || strstr(strbuf1, "20")) {
		/* date field */
		sscanf(linebuf, "%d-%d-%d %d:%d:%d", &year, &month, &day, &hour, &minute, &second);
		year -= 1900;
		xtime[Ndata].tm_year = year;
		xtime[Ndata].tm_mon = month-1;
		xtime[Ndata].tm_mday = day;
		xtime[Ndata].tm_hour = hour;
		xtime[Ndata].tm_min = minute;
		xtime[Ndata].tm_sec = second;
    		mktime(&xtime[Ndata]);		/* this rebuilds all the fields */
	    }
	    else if(strncmp(strbuf1, CITYNAME, 3) == 0) {
		vt=-999.0; vdew=-999.0; vrh = -999.0;	/* init to invalid values */
		sscanf(linebuf+15, "%s %f %f %f", sky, &vt, &vdew, &vrh);
		if(strncmp(sky, "NOT", 3) == 0 || strncmp(sky, "N/A", 3) == 0
			|| strncmp(sky, "MISG", 4) == 0 ||strncmp(sky, "MM", 2) == 0)
				continue;
		if(vt < -990.0 || vdew < -990.0 || vrh < -990.0)
				continue;
		/* OK, we are pretty sure we have valid data */
		if(useC)
		    temp[Ndata] = 0.5555555*(vt - 32.0);  /* convert to centigrade */
		else
		    temp[Ndata] = vt;
		humidity[Ndata] = vrh;
		Ndata++;
	    }
	} /* end while */

	fclose(fp);
	fprintf(stderr, "%d records read.\n", Ndata);
	return(0);
}

void do_TextStuff(CPDFdoc *pdf, int useC)
{
int i = Ndata -1;
char tbuf[1024];
    cpdf_setgray(pdf, 0.0);
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "HeiseiMin-W3", "90ms-RKSJ-H", 22.0);
    cpdf_text(pdf, 1.0, 9.8, 0.0, "オークランド国際空港の天気 (カリフォルニア州）");
    cpdf_endText(pdf);

    cpdf_setgray(pdf, 0.0);
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "HeiseiKakuGo-W5", "90ms-RKSJ-H", 14.0);
    sprintf(tbuf, "年月日: %d-%d-%d     時刻: %2d:00",
	xtime[i].tm_year+1900, xtime[i].tm_mon+1, xtime[i].tm_mday, xtime[i].tm_hour);
    cpdf_text(pdf, 1.7, 9.4, 0.0, tbuf);
    sprintf(tbuf, "気温: %.1f [%s],    相対湿度: %.0f [%%]",
	temp[i], tempUnit[useC], humidity[i]);
    cpdf_text(pdf, 1.7, 9.1, 0.0, tbuf);
    sprintf(tbuf, "空模様: %s", sky);
    cpdf_text(pdf, 1.7, 8.8, 0.0, tbuf);
    cpdf_endText(pdf);

}

/* -------------------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
char *gfile;
int   useCentigrade = 0;		/* if non-zero, plot in centigrade */
CPDFdoc *pdf;

	if(argc < 2) {
	    printf("Usage: %s datafile [C|F]\n", argv[0]);
	    printf("        Arg2 for temp unit: C=centigrade, F=Fahrenheit\n");
	}

	if(argc == 3) {
	    if(argv[2][0] == 'C' || argv[2][0] == 'c')
		useCentigrade =1;
	}

	gfile = argv[1];				/* use arg1 as filename */
	read_datafile(gfile, useCentigrade);

	/* == OK, data file read in. Start plotting */
	pdf = cpdf_open(0, NULL);
	cpdf_enableCompression(pdf, 1);		/* use Flate/Zlib compression */
	cpdf_init(pdf);
	cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */

	/* PDF Info object */
	cpdf_setCreator(pdf, "weather.c");
	cpdf_setTitle(pdf, "Weather of OAK California");
	cpdf_setSubject(pdf, "Weather of Oakland International Airport, California");
	cpdf_setKeywords(pdf, "Weather");

	/* Your drawing code in PDF below */
	/* ----------------------------------------------------------------------------------------- */

	/* plot curves */
	plot_temp(pdf, useCentigrade);
	plot_rHumidity(pdf);
	do_TextStuff(pdf, useCentigrade);

	/* ----------------------------------------------------------------------------------------- */
	/* Your drawing code in PDF above */
	cpdf_finalizeAll(pdf);			/* PDF file is actually written here */
	cpdf_savePDFmemoryStreamToFile(pdf, "weatherj.pdf");

	cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */
	cpdf_close(pdf);			/* PDF file is actually written here */
	return(0);
}

