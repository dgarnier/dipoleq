/*
 * Copyright (C) 1998, 1999 FastIO Systems, All Rights Reserved.
 * For conditions of use, license, and distribution, see LICENSE.pdf
 * included in the distribution or http://www.fastio.com/LICENSE.pdf

cc -Wall -s -object -o hex hex.c
cc -Wall -o hex hex.c

*/

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>

unsigned char _cpdf_nibbleValue(unsigned char hexchar)
{
unsigned char val;
    if(hexchar >= '0' && hexchar <= '9')
	val = hexchar - '0';
    else if(hexchar >= 'A' && hexchar <= 'F')
	val = hexchar - 'A' + 10;
    else if(hexchar >= 'a' && hexchar <= 'f')
	val = hexchar - 'a' + 10;
    else
	val = 0xff;
    return(val);
}

/* binary to HEX conversion function -- for Unicode support */
static char *_cpdf_HexTab= "0123456789ABCDEF00000000";

char *cpdf_convertBinaryToHex(const unsigned char *datain, char *hexout, long length, int addFEFF)
{
long i;
unsigned char *pin = (unsigned char *)datain;
unsigned short *pcheck = (unsigned short *)datain;
char *pout = hexout;
unsigned char ch;
    if(addFEFF && (*pcheck != 0xfeff && *pcheck != 0xfffe)) {
	*pout++ = 'F';
	*pout++ = 'E';
	*pout++ = 'F';
	*pout++ = 'F';
    }
    for(i=0; i < length; i++) {
	ch = *pin >> 4;
	*pout++ = _cpdf_HexTab[ch];
	ch = *pin++ & 0x0f;
	*pout++ = _cpdf_HexTab[ch];
    }
    *pout = '\0';
    return(hexout);
}

/* After the call, length will contain the byte length of binary data. */
unsigned char *cpdf_convertHexToBinary(const char *hexin, unsigned char *binout, long *length)
{
long bcount=0;
char *pin = (char *)hexin;
unsigned char *pout = binout;
unsigned char ch, chin;
int HighNibble = 1;
    while( (chin = *pin++) != '\0') {
	ch = _cpdf_nibbleValue(chin);
	if(ch > 15)
	    continue;		/* not a HEX character */

	if(HighNibble) {
	    *pout = ch << 4;
    	    HighNibble = 0;
	}
	else {
	    *pout++ |= ch;
    	    HighNibble = 1;
	    bcount++;
	}
    }
    *length = bcount;	/* byte count */
    return(binout);
}

/* Returns the size of file in bytes */

long getFileSize(const char *file)
{
long filesize = 0;
#if defined(_WIN32) || defined(WIN32)
    #ifndef __BORLANDC__
        struct _stat filestat;		/* For Windows VC++ */
    #else
        struct stat filestat;
    #endif
#else
struct stat filestat;
#endif

#if defined(_WIN32) || defined(WIN32)
	if(_stat(file, &filestat))	/* Windows */
#else
	if(stat(file, &filestat))
#endif
	    filesize = 0;		/* stat failed -- no such file */
	else {
	    if( !(filestat.st_mode & S_IFREG) )
		filesize = 0;	/* not a regular file */
#ifdef UNIX
	    else if( !(filestat.st_mode & S_IREAD) )
		filesize = 0;	/* not readable */
#endif
	    else
		filesize = filestat.st_size;
	}
	return(filesize);
}

void main(int argc, char *argv[])
{
unsigned char *data;
char *hex;
long length;
FILE *fp;
	length = getFileSize(argv[1]);
	if(length <= 0) exit(1);
	data = (unsigned char *)malloc(length+2);
	hex = (char *) malloc(2*length+4);
	fp = fopen(argv[1], "r");
	if(fp ==NULL) exit(2);
	fread(data, 1, length, fp);
	fclose(fp);
	cpdf_convertBinaryToHex(data, hex, length, 1);
	printf("%s\n\n", hex);
}

/*
void main(int argc, char *argv[])
{
unsigned char data[256], hex[518], dataout[260];
int i;
long length;

    for(i=0; i<256; i++)
	data[i] = i;
    cpdf_convertBinaryToHex(data, hex, 256, 1);
    printf("%s\n\n", hex);

    cpdf_convertHexToBinary(hex, dataout, &length);
    for(i=0; i<length; i++)
	printf("%d\n", dataout[i]);
}
*/

