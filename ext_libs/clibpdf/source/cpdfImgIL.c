/* cpdfImgIL.c  -- In-line image functions.
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

	Image data are written to the page content stream without using the XObject.
	Only for small images.

1999-08-17 [IO] 
	## MULTI_THREADING_OK ##

1998-09-27 [IO]
*/

/* #define CPDF_DEBUG 1 */

#include "config.h"
#include "version.h"

#include <string.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include "cpdflib.h"		/* This must be included before all other local include files */
#include "cpdf_private.h"
#include "cglobals.h"

int cpdf_placeInLineImage(CPDFdoc *pdf, const void *imagedata, int len,
		float x, float y, float angle, float width, float height,
		int pixwidth, int pixheight, int bitspercomp, int CSorMask, int gsave)
{
    cpdf_rawPlaceInLineImage(pdf, imagedata, len,
		x_Domain2Points(pdf, x), y_Domain2Points(pdf, y), angle, width, height,
		pixwidth, pixheight, bitspercomp, CSorMask, gsave);
    return(0);
}

int cpdf_rawPlaceInLineImage(CPDFdoc *pdf, const void *imagedata, int len,
		float x, float y, float angle, float width, float height,
		int pixwidth, int pixheight, int bitspercomp, int CSorMask, int gsave)
{
int count, bufsize;
char *mbuff;

    if(gsave) cpdf_gsave(pdf);
    cpdf_rawTranslate(pdf, x, y);
    if(fabs(angle) > 0.001)
	cpdf_rotate(pdf, angle);
    cpdf_rawConcat(pdf, width, 0.0, 0.0, height, 0.0, 0.0);	/* concatename CTM */

    cpdf_clearMemoryStream(pdf->scratchMem);
    pdf->inlineImages++;			/* For enabling Image* ProcSet */
    cpdf_memPuts("BI\n", pdf->scratchMem);		/* begin image */

    sprintf(pdf->spbuf, "/W %d\n/H %d\n/BPC %d\n", pixwidth, pixheight, bitspercomp);
    cpdf_writeMemoryStream(pdf->scratchMem, pdf->spbuf, strlen(pdf->spbuf));

    switch(CSorMask & 3) {
	    case IMAGE_MASK:
			cpdf_memPuts("/ImageMask true\n", pdf->scratchMem);
	    		pdf->imageFlagBCI |= 4;		/* /ImageI procset */
			break;
	    case CS_GRAY:
	    		pdf->imageFlagBCI |= 1;		/* /ImageB procset */
			cpdf_memPuts("/ColorSpace /DeviceGray\n", pdf->scratchMem);
			break;
	    case CS_RGB:
			cpdf_memPuts("/ColorSpace /DeviceRGB\n", pdf->scratchMem);
	    		pdf->imageFlagBCI |= 2;		/* /ImageC procset */
			break;
	    case CS_CMYK:
			cpdf_memPuts("/ColorSpace /DeviceCMYK\n", pdf->scratchMem);
	    		pdf->imageFlagBCI |= 2;		/* /ImageC procset */
			break;
    }
    cpdf_memPuts("ID ", pdf->scratchMem);		/* the space char after ID is important */
    /* OK, here we dump the data in binary.  No need to do any compression
	because the page content stream itself is compressed by Flate
	compression anyway (if it is enabled).
    */
    cpdf_writeMemoryStream(pdf->scratchMem, (char *) imagedata, len);

    cpdf_memPuts("\nEI\n", pdf->scratchMem);		/* end image */

    cpdf_getMemoryBuffer(pdf->scratchMem, &mbuff, &count, &bufsize);
    if(pdf->useContentMemStream)
	cpdf_writeMemoryStream(pdf->currentMemStream, mbuff, count);
    else
	fwrite(mbuff, 1, (size_t)count, pdf->fpcontent);

    if(gsave) cpdf_grestore(pdf);

    return(0);
}
