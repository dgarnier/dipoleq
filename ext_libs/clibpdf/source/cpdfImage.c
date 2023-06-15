/* cpdfImage.c
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or at http://www.fastio.com/LICENSE.pdf

1999-11-17 [io]
	/ImageMask support for B/W images.  Useful for watermarks.
1999-08-17 [IO]
	## MULTI_THREADING_OK ##
1998-09-03 [IO]
	JPEG (baseline) image import.
-------------------------------------------------------------------------------------------
Usage Notes on cpdf_importImage(), and cpdf_rawImportImage()

   These two functions are identical except for the domain in which the coordinate
   system for (x, y) is defined.  With  cpdf_importImage(), (x, y) are defined
   in the current domain coordinate system.  Width and height should still be
   defined in points.  With cpdf_rawImportImage(), all of (x, y, width, height)
   are defined in points in the default coordinate system.

   When you import an image file, you may not know image dimensions (in pixels).
   For each X or Y dimension, specify either size (in points) or scaling factor
   and set the unspecified one to zero.
   The unspecified variable will be set upon return (that's why these parameters
   are passed by address).  If both are set, width and height will take
   precedence over scaling factors.  For example if the image file "/tmp/myimage.jpg" has
   the dimension of 150x100 pixels, a call:
   width = 0.0;
   height = 200.0;
   xscale = 0.0;
   yscale = 0.0;
   cpdf_rawImportImage(pdf, "/tmp/myimage.jpg", JPEG_IMG, 100.0, 100.0, &width, &height, &xscale, &yscale, 1);

   will return with width = 300, and xscale = yscale = 2.0, while the aspect ratio is
   preserved.  If the ratio of width/height does not match the aspect ratio of the original
   image, the aspect ratio will not be preserved.

   ** Important Note:
	When the last argument is non-zero, the call brackets image import with
	a cpdf_flags() and cpdf_grestore() pair.  If the last argument is zero,
	the function does not perform flags/grestore internally to allow flexibility
	of having an arbitrary path clip the image or to draw image inside text
	(See Example 14.4 in the  PDF Reference Manual version 1.2).

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

#define M_SOF0  0xC0		/* Baseline JPEG */
#define M_SOF1  0xC1		/* Extended sequential, Huffman */
#define M_SOF2  0xC2		/* Progressive JPEG */

/* Same as cpdf_rawImportImage(), but (x, y) are defined in the current domain
   coordinate system.  Width and height should still be defined in points.
   If 'gsave' is non-zero, it performs gsave/grestore internally.

   #### NOTE on Multi-page TIFF #############################################################
   If info on the non-first page of a multi-page TIFF is needed, cpdf_selectImage() must be
   called first to set the page index (0.. N-1) before calling the functions below.
*/


int cpdf_importImage(CPDFdoc *pdf, const char *imagefile, int type, float x, float y, float angle,
	float *width, float *height, float *xscale, float *yscale, int flags)
{
    return( cpdf_rawImportImage(pdf, imagefile, type, x_Domain2Points(pdf, x), y_Domain2Points(pdf, y), angle,
		width, height, xscale, yscale, flags) );
}


/* same as cpdf_importImage(), but (x,y) are in points.
*/

int cpdf_rawImportImage(CPDFdoc *pdf, const char *imagefile, int type, float x, float y, float angle,
	float *width, float *height, float *xscale, float *yscale, int flags)
{
int foundInPageList = 0;	/* in pageInfos[].imageIdx[] */
int retcode=0;
int imageIndexFound = 0;	/* in imageInfos[] */
int currentImageIndex = 0;
float a, b=0.0, c=0.0, d, e, f;
CPDFimageInfo *newImage;
char imagename[32], *image;

    if(type != TIFF_IMG)
	pdf->imageSelector = 0;	/* only TIFF images can have multiple images in a file */

    /* check if this image has already been included or new image */
    if( _isNewImage(pdf, imagefile, pdf->imageSelector, &imageIndexFound) ) {

	if(pdf->numImages >= pdf->nMaxImages) {
	    cpdf_Error(pdf, 1, "ClibPDF", "Too many images in this PDF: %d. Increase limits by cpdf_open()",
			pdf->numImages);
	    return(1);
	}
	newImage = &pdf->imageInfos[pdf->numImages];
	newImage->data = NULL;
	if(flags & 2)			/* bit-1 of flags is /ImageMask */
	    newImage->imagemask = 1;
	else
	    newImage->imagemask = 0;
	if(flags & 4)			/* bit-2 of flags is invert gray */
	    newImage->invert = 1;
	else
	    newImage->invert = 0;
	newImage->photometric = 1;	/* min is black - tiff */
	newImage->orientation = 1;	/* normal orientation - tiff */
	newImage->imgSelector = pdf->imageSelector;
	newImage->filesize = getFileSize(imagefile);	/* used only for JPEG, overwritten for others */
	currentImageIndex = pdf->numImages;
	/* Process different image file formats here */
	switch(type) {
	    case JPEG_IMG:
		if((retcode = read_JPEG_header(imagefile, newImage))) {
		    if(retcode == -1) {
			return(retcode);	/* cpdf_Error() done already in function */
		    }
		    if(retcode == -2) {
			cpdf_Error(pdf, 1, "ClibPDF", "Not a JPEG file?: %s", imagefile);
			return(retcode);
		    }
		    if( newImage->process != M_SOF0 ) {
			/* Acrobat 4.0x and later can accept progressive JPEG, so don't refuse those
			   just set the PDF level to 1.3 to indicate appropriate requirement.
			*/
			/* cpdf_Error(pdf, 1, "ClibPDF", "JPEG file is not baseline type: %s", imagefile); */
			/* return(-3); */
			cpdf_setPDFLevel(pdf, 1, 3);	/* PDF-1.3 */
		    }
		}
		break;

	    case TIFF_IMG:
		if((retcode = cpdf_readTIFFheader(imagefile, newImage))) {
		    if(retcode)
			return(retcode);
		}
		break;

	    case GIF_IMG:	/* It won't be implemented unless we decide to license LZW patent. */
		break;

	    case CPDF_IMG:	/* PDFIMG - PDF Image XObject dump format */
		if((retcode = cpdf_readPDFimageHeader(imagefile, newImage))) {
		    if(retcode < 0)
			return(retcode);
		}
		break;

	    default:
		/* cpdf_Error(); */
		return(-4);
		break;
	} /* end switch(type) {} */

	if(newImage->ncomponents == 1)
	    pdf->imageFlagBCI |= 1;		/* /ImageB procset */
	else if(newImage->ncomponents > 2 )
	    pdf->imageFlagBCI |= 2;		/* /ImageC procset */
        sprintf(imagename, "IMcpdf%d", pdf->numImages);
	image = imagename;
	/* memcpy(newImage, &tmpImage, sizeof(CPDFimageInfo)); */
	newImage->name = (char *)malloc((size_t)(strlen(imagename) + 1));
	_cpdf_malloc_check((void *)newImage->name);
	newImage->filepath = (char *)malloc((size_t)(strlen(imagefile) + 1));
	_cpdf_malloc_check((void *)newImage->filepath);
	strcpy(newImage->name, imagename);
	strcpy(newImage->filepath, imagefile);
	newImage->type = type;
	pdf->numImages++;
    }
    else {	/* else clause of if( _isNewImage() ) .. */
	/* This image has been used already, and is in imageInfos[imageIndexFound] */
	image = pdf->imageInfos[imageIndexFound].name;
	currentImageIndex = imageIndexFound;
    }

    /*  Store the index of this image in the used image list in pageInfos[].
	Multiple pages can share images.
    */
    {   /* braces for local variable */
	CPDFintList *iLprev, *iL = pdf->pageInfos[pdf->currentPage].imageIdx;
	iLprev = iL;
	foundInPageList = 0;
	while(iL != NULL) {
	    if(iL->value == currentImageIndex) {
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
	    iLprev->value = currentImageIndex;	/* store the index to current font */
	    iL->next = NULL;
	    iL->value = -1;
	    pdf->pageInfos[pdf->currentPage].npImage++;
	}
    }

/*
    fprintf(stderr, "Image file: %s\n", pdf->imageInfos[imageIndexFound].filepath);
    fprintf(stderr, "File size: %ld\n", pdf->imageInfos[imageIndexFound].filesize);
    fprintf(stderr, "Resolution: (%d x %d)\n", pdf->imageInfos[imageIndexFound].width,
				pdf->imageInfos[imageIndexFound].height);
*/

    /* work out positioning and scaling for use with concat operator */
    e = x;
    f = y;
    if(fabs(*width) > 0.00001) {
	/* use the width value passed */
	a = *width;			/* scale unit image to this */
	*xscale = a / pdf->imageInfos[imageIndexFound].width;
    }
    else if( fabs(*xscale) > 0.00001) {
	/* width not specified, but we have xscale */
	a = *xscale * pdf->imageInfos[imageIndexFound].width;
	*width = a;
    }
    /* neither width or xscale is set, see if Y specs are present */
    else if( fabs(*height) > 0.00001) {
 	/* Nothing set for X, but height is set, so use the Y scaling factor for xscale */
	*yscale = *height / pdf->imageInfos[imageIndexFound].height;
	*xscale = *yscale;
	a = *xscale * pdf->imageInfos[imageIndexFound].width;
    }
    else if( fabs(*yscale) > 0.00001) {
	/* Nothing set for X, but y scale is set, so use that for xscale too */
	*xscale = *yscale;
	a = *xscale * pdf->imageInfos[imageIndexFound].width;
	*width = a;
    }
    else {
	/* Nothign is specified. Use image xpixel count as number of points */
	a = (float)pdf->imageInfos[imageIndexFound].width;
	*width = a;
	*xscale = 1.0;
    }

    /* NOW do the same for y scaling and height */
    if(fabs(*height) > 0.00001) {
	/* use the width value passed */
	d = *height;			/* scale unit image to this */
	*yscale = d / pdf->imageInfos[imageIndexFound].height;
    }
    else if( fabs(*yscale) > 0.00001) {
	/* height not specified, but we have yscale */
	d = *yscale * pdf->imageInfos[imageIndexFound].height;
	*height = d;
    }
    /* Nothing is set for Y, see if there is X scaling info we can use */
    else if( fabs(*width) > 0.00001) {
 	/* Nothing set for Y, but width is set, so use the X scaling factor for yscale */
	*xscale = *width / pdf->imageInfos[imageIndexFound].width;
	*yscale = *xscale;
	d = *yscale * pdf->imageInfos[imageIndexFound].height;
    }
    else if( fabs(*xscale) > 0.00001) {
	/* Nothing set for Y, but x scale is set, so use that for yscale too */
	*yscale = *xscale;
	d = *yscale * pdf->imageInfos[imageIndexFound].height;
	*height = d;
    }
    else {
	/* Really nothing is set, so use image xpixel count as number of points */
	d = (float)pdf->imageInfos[imageIndexFound].height;
	*height = d;
	*yscale = 1.0;
    }

    /* Write image use line to current Contents stream */
    if(flags & 1) cpdf_gsave(pdf);		/* bit-0 is gsave flag */
    cpdf_rawTranslate(pdf, x, y);
    if(fabs(angle) > 0.001)
	cpdf_rotate(pdf, angle);

    /* fprintf(stderr, "CTM = [%g %g %g %g %g %g]\n", a, b, c, d, e, f); */
    cpdf_rawConcat(pdf, a, b, c, d, 0.0, 0.0);	/* concatename CTM */

    if(pdf->useContentMemStream) {
	sprintf(pdf->spbuf, "/%s Do\n", image);
	cpdf_writeMemoryStream(pdf->currentMemStream, pdf->spbuf, strlen(pdf->spbuf));
    }
    else
        fprintf(pdf->fpcontent, "/%s Do\n", image);
    if(flags & 1) cpdf_grestore(pdf);
    return(0);
}


int _isNewImage(CPDFdoc *pdf, const char *filepath, int imgSel, int *imageFound)
{
int i, isNew = 1;
    /* look in the list */
    for(i=0; i < pdf->numImages; i++) {
	if( strcmp(filepath, pdf->imageInfos[i].filepath) == 0 && imgSel == pdf->imageInfos[i].imgSelector) {
	    isNew = 0;		/* already defined */
	    *imageFound = i;	/* return found image index */
	    break;
	}
    }
    if(isNew)
	*imageFound = pdf->numImages;		/* index of the new image is this */
    return(isNew);
}

/* Set index (0..N-1) to a current image in multi-page image format such as TIFF */
int cpdf_selectImage(CPDFdoc *pdf, int imgsel)
{
int oldIndex;
    oldIndex = pdf->imageSelector;
    if(imgsel >= 0)
        pdf->imageSelector = imgsel;
    return(oldIndex);
}


int  _cpdf_freeAllImageInfos(CPDFdoc *pdf)
{
CPDFimageInfo *timage;
int i;
    for(i=0; i< pdf->numImages; i++) {
	timage = &pdf->imageInfos[i];
	if(timage->name) {
	    free(timage->name);
	    timage->name = NULL;
	}
	if(timage->filepath) {
	    free(timage->filepath);
	    timage->filepath = NULL;
	}
	if(timage->data) {
	    free(timage->data);
	    timage->data = NULL;
	}
    }
    return(0);
}

#ifndef USE_LIBTIFF
/* Dummy functions just for stubs when USE_LIBTIFF is not defined. */


int cpdf_readTIFFheader(const char *imagefile, CPDFimageInfo *timgInfo)
{
    cpdf_GlobalError(1, "ClibPDF", "LIBTIFF library support is not compiled in.");
    return(-1);
}

int cpdf_readTIFFimageData(char **rasterBuf, CPDFimageInfo *timgInfo, int compress)
{
    /* cpdf_GlobalError(1, "ClibPDF", "LIBTIFF library support is not compiled in."); */
    return(-1);
}

#endif  /* USE_LIBTIFF */

#ifndef USE_PDFIMG
/* Dummy functions just for stubs when USE_PDFIMG is not defined. */
int cpdf_readPDFimageHeader(const char *imagefile, CPDFimageInfo *timgInfo)
{
    cpdf_GlobalError(1, "ClibPDF", "PDFIMG support is not compiled in.");
    return(-1);
}

int _cpdf_copyPDFimageHeader(CPDFdoc *pdf, CPDFimageInfo *timgInfo)
{
    return(-1);
}

int cpdf_readPDFimageData(char **rasterBuf, CPDFimageInfo *timgInfo, int ImDataOffset)
{
    return(-1);
}


#endif  /* USE_PDFIMG */
