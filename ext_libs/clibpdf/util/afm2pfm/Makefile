# Makefile for pfm2afm and afm2pfm
# for EMX/GCC 0.8h
# Russell Lang <rjl@eng.monash.edu.au>

EMXDIR=d:/emx

ZIPFILES=README COPYING Makefile Makefile.bc pfm2afm.c pfm2afm.h pfm2afm.exe afm2pfm.c afm2pfm.h afm2pfm.exe

all: afm2pfm.exe pfm2afm.exe

afm2pfm.exe: afm2pfm
	emxbind $(EMXDIR)/bin/emxl.exe afm2pfm afm2pfm.exe

afm2pfm: afm2pfm.c afm2pfm.h
	cc -g -o afm2pfm afm2pfm.c

pfm2afm.exe: pfm2afm
	emxbind $(EMXDIR)/bin/emxl.exe pfm2afm pfm2afm.exe

pfm2afm: pfm2afm.c pfm2afm.h
	cc -g -o pfm2afm pfm2afm.c

clean:
	-del afm2pfm
	-del afm2pfm.o
	-del afm2pfm.exe
	-del pfm2afm
	-del pfm2afm.o
	-del pfm2afm.exe

zip: $(ZIPFILES)
	-del afm2pfm.zip
	zip afm2pfm.zip $(ZIPFILES)

