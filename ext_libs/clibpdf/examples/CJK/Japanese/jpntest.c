/* jpntest.c -- test program for Japanese language support.
 * Copyright (C) 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf
------------------------------------------------------------------------------------
1999-10-01 [io] Reorganized, and changed default unit to centimeter.
1999-09-25 [io] Starting with Japanese, and then Korean.

cc -Wall -o jpntest -I/usr/local/include jpntest.c -lcpdfpm
*/

#include <stdio.h>
#include <stdlib.h>
#include <cpdflib.h>

/* 90ms-RKSJ-H (Shift-JIS) text */
char *sjis_text[] = {
"ClibPDFをあなたのプログラムに組み込むことにより、",
"アドビ・システムズのアクロバット等による手作業",
"無しでＰＤＦファイルを自動的に作成することができます。",
"ウエブ・サーバ上でユーザからの入力、サーバサイドでの",
"データベース・アクセスによって、一定の書式の書類を",
"大量にＰＤＦファイルとして作るような用途に最適です。",
NULL
};

char *sjis_english[] = {
"By linking ClibPDF into your program,",
"It is possible to generate PDF files automatically",
"without manually operating Adobe Systems' Acrobat software",
"and other related plug-ins.  It is suitable for generating",
"a large number of documents of a given format based on user",
"input from Web browsers and server-side data from databases.",
NULL
};

int do_Paragraph(CPDFdoc *pdf, float *x, float *y, float ybump, int alignmode, float fontsize,
			const char *font, const char *encoding, const char *paragraph[]);

int main(int argc, char *argv[])
{
CPDFdoc *pdf;
float fontsize = 12.0;	/* in points */
float x = 2.5;
float y = 27.0;		/* in centimeters */
int alignmode = TEXTPOS_LL;
float ybump = 0.5;	/* in cm */

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);
    cpdf_setDefaultDomainUnit(pdf, POINTSPERCM);	/* change default unit to cm */
    cpdf_pageInit(pdf, 1, PORTRAIT, A4, A4);	/* orientation and size. It's A4 in Japan */

    /* Try 5 different font variations out of 2 base fonts: HeiseiKakuGo-W5 and HeiseiMin-W3. */
    do_Paragraph(pdf, &x, &y, ybump, alignmode, fontsize,
			"HeiseiKakuGo-W5", "90ms-RKSJ-H", sjis_text);

    do_Paragraph(pdf, &x, &y, ybump, alignmode, fontsize,
			"HeiseiKakuGo-W5,Italic", "90ms-RKSJ-H", sjis_text);

    do_Paragraph(pdf, &x, &y, ybump, alignmode, fontsize,
			"HeiseiMin-W3", "90ms-RKSJ-H", sjis_text);

    do_Paragraph(pdf, &x, &y, ybump, alignmode, fontsize,
			"HeiseiMin-W3,Bold", "90ms-RKSJ-H", sjis_text);

    do_Paragraph(pdf, &x, &y, ybump, alignmode, fontsize,
			"HeiseiMin-W3,BoldItalic", "90ms-RKSJ-H", sjis_text);

    /* Now, do English translation. */
    do_Paragraph(pdf, &x, &y, ybump, alignmode, fontsize,
			"Helvetica", "WinAnsiEncoding", sjis_english);

    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "jpntest.pdf");
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */

    /* == Clean up == */
    cpdf_close(pdf);			/* shut down */
    return 0;
}

int do_Paragraph(CPDFdoc *pdf, float *x, float *y, float ybump, int alignmode, float fontsize,
			const char *font, const char *encoding, const char *paragraph[])
{
int i;
char buff[1024];
    *y -= ybump;
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Times-Italic", "WinAnsiEncoding", 10.0);
    sprintf(buff, "[Font: %s      Encoding (CMap): %s ]", font, encoding);
    cpdf_textAligned(pdf, *x, *y, 0.0, alignmode, buff);
    *y -= ybump;

    cpdf_setFont(pdf, font, encoding, fontsize);
    for(i = 0; paragraph[i] != NULL; i++) {
	cpdf_textAligned(pdf, *x, *y, 0.0, alignmode, paragraph[i]);
	*y -= ybump;
    }
    cpdf_endText(pdf);
    return 0;
}

