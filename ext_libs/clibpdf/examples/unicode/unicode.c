/* unicode.c -- test program for Unicode and Japanese language support.
 * Copyright (C) 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf
------------------------------------------------------------------------------------
1999-11-27 [io] Unicode example for outlines and annotations.
1999-10-01 [io] Reorganized, and changed default unit to centimeter.
1999-09-25 [io] Starting with Japanese, and then Korean.

cc -Wall -o unicode -I/usr/local/include unicode.c -lcpdfpm
*/

#include <stdio.h>
#include <stdlib.h>
#include <cpdflib.h>

/* In this example, I am storing Unicode text as HEX strings, because
original Unicode data cannot be represented easily as C strings.
In real use, you will be using cpdf_convertBinaryToHex() function to
convert binary Unicode data to HEX format, and then passing that
to ClibPDF text functions.
*/

/* Unicode (UniJIS-UCS2-H or UniJIS-UCS2-HW-H) Japanese text in HEX */
/* two outline items */
char *pdfref_adobe = "FEFF00500044004630EA30D530A130EC30F330B9FF2030A230C930D3";
char *bookmark2 = "FEFF3057304A308AFF12FF1A30DA30FC30B8306E4E0B65B93078";

/* annotation */
char *annot_title = "FEFF65E5672C8A9E306E6CE891C8306E4F8B";
char *annotUnicode = "FEFF0043006C0069006200500044004630D030FC30B830E730F30032\
002E003000323088308A65E5672C8A9E306E6CE891C830843057304A308A304C30B530DD30FC\
30C83055308C307E30593002000D000D0043006C0069006200500044004630676CE891C80028\
0061006E006E006F0074006100740069006F006E002930843057304A308A0028006F00750074\
006C0069006E0065002F0062006F006F006B006D00610072006B002930924ED83051308B305F\
3081306B306F300100480065007800610064006500630069006D0061006C5F625F0F306E0055\
006E00690063006F00640065306E30C630AD30B930C830924F7F3044307E30593002000D000D\
0063007000640066005F0068006500780053007400720069006E0067004D006F006400650028\
0059004500530029003B002030924F7F3063306630C630AD30B930C800410050004930920048\
0065007830E230FC30C9306B30573066304B30893001901A5E38306E00410050004930924F7F\
306330664E0B305530443002000D000D00480065007800610064006500630069006D0061006C\
5F625F0F3078306E590963DB306E305F3081306E95A2657030010063007000640066005F0063\
006F006E007600650072007400420069006E0061007200790054006F00480065007800280029\
002C00200063007000640066005F0063006F006E00760065007200740048006500780054006F\
00420069006E006100720079002800297B49304C65B0305F306B4ED8305152A030483089308C\
307E3057305F3002000D";

/* This one without "FEFF" - byte-order identifier */
char *hexUnicodeContent = "0043006C00690062005000440046002030D030FC30B830E730F30032002E00300032\
00203088308A65E5672C8A9E306E6CE891C8304C30B530DD30FC30C83055308C307E30593002";

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
float y = 27.5;		/* in centimeters */
int alignmode = TEXTPOS_LL;
float ybump = 0.5;	/* in cm */
float p1=200.0, p2=0.0, p3=0.0, p4=0.0;
CPDFoutlineEntry *currTopOL;
CPDFviewerPrefs vP = { PM_OUTLINES, 0, 0, 0, 0, 0, 0, 0 };
char *act_dict;

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);
    /* Have outline (book marks) pane visible upon document open */
    cpdf_setViewerPreferences(pdf, &vP);
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

    /* Test Unicode for annotation title and text */
    cpdf_hexStringMode(pdf, YES);
    cpdf_rawSetAnnotation(pdf, 0.0, 100.0, 300.0, 350.0, annot_title, annotUnicode, NULL);

    /* Try Unicode in Outline (bookmarks) */
    act_dict = "/Type /Action /S /URI /URI (http://partners.adobe.com/asn/developer/PDFS/TN/PDFSPEC.PDF)";
    currTopOL = cpdf_addOutlineAction(pdf, NULL, OL_SAME, OL_OPEN, act_dict, pdfref_adobe);
    currTopOL = cpdf_addOutlineEntry(pdf, currTopOL, OL_SAME, OL_OPEN, 1,
					bookmark2, DEST_Y, p1, p2, p3, p4);

    /* Try Unicode in page content text */
    y -= ybump;
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "HeiseiMin-W3,Italic", "UniJIS-UCS2-HW-H", fontsize);
    /* cpdf_textAligned() doesn't work yet */
    cpdf_text(pdf, x, y, 0.0, hexUnicodeContent);	/* No need for "FEFF" for page content? */
    cpdf_endText(pdf);
    cpdf_hexStringMode(pdf, NO);	/* turn off Hex mode before you forget */

    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "unicode.pdf");
    /* cpdf_launchPreview(pdf); */		/* launch Acrobat/PDF viewer on the output file */

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

