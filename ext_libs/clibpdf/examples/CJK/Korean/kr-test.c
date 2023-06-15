/* kr-test.c -- test program for Korean language support.
 * Copyright (C) 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf
------------------------------------------------------------------------------------
1999-09-25 [io] Starting with Japanese, and then Korean.
		We will need Chinese help soon.

cc -Wall -o kr-test -I/usr/local/include kr-test.c -lcpdfpm

*/

#include <stdio.h>
#include <stdlib.h>
#include <cpdflib.h>

/* KSC-EUC-H (EUC-KR) test : Curtesy of Dr. Hwang, Jinsoo */
char *euc_kr_text[] = {
"이글은 한글윈도우 95에서 노트패드상에서",
"작성한 글입니다. 이렇게 한글에 대하여도",
"신경을 써주셔서 감사합니다.",
"황 진 수",
NULL
};

/* KSC-EUC-H (EUC-KR) test : Curtesy of Dr. Hwang, Jinsoo */
char *euc_kr_english[] = {
"English translation. :-)",
"This is written under Korean Windows95 with NotePad.",
"I appreciate your efforts to take care of Korean characters.",
"Thank you",
"",
"Hwang, Jinsoo",
NULL
};

int main(int argc, char *argv[])
{
int i;
CPDFdoc *pdf;
float fontsize = 16.0;
float x = 1.0;
float y = 8.5;
int alignmode = TEXTPOS_LL;
float ybump = 0.25;

    /* == Initialization == */
    pdf = cpdf_open(0, NULL);
    cpdf_enableCompression(pdf, YES);		/* use Flate/Zlib compression */
    cpdf_init(pdf);
    cpdf_pageInit(pdf, 1, PORTRAIT, LETTER, LETTER);		/* page orientation */

    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "HYGoThic-Medium", "KSC-EUC-H", fontsize);	/* Gothic type font */
    for(i = 0; euc_kr_text[i] != NULL; i++) {
	cpdf_textAligned(pdf, x, y, 0.0, alignmode, euc_kr_text[i]);
	y -= ybump;
    }
    cpdf_endText(pdf);

    y -= ybump;
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "HYSMyeongJo-Medium", "KSC-EUC-H", fontsize);	/* MyeongJo type font */
    for(i = 0; euc_kr_text[i] != NULL; i++) {
	cpdf_textAligned(pdf, x, y, 0.0, alignmode, euc_kr_text[i]);
	y -= ybump;
    }
    cpdf_endText(pdf);

    y -= ybump;
    cpdf_beginText(pdf, 0);
    cpdf_setFont(pdf, "Helvetica", "MacRomanEncoding", fontsize);
    for(i = 0; euc_kr_english[i] != NULL; i++) {
	cpdf_textAligned(pdf, x, y, 0.0, alignmode, euc_kr_english[i]);
	y -= ybump;
    }
    cpdf_endText(pdf);

    cpdf_finalizeAll(pdf);			/* PDF file/memstream is actually written here */
    cpdf_savePDFmemoryStreamToFile(pdf, "kr-test.pdf");
    cpdf_launchPreview(pdf);		/* launch Acrobat/PDF viewer on the output file */

    /* == Clean up == */
    cpdf_close(pdf);			/* shut down */
    return 0;

}


