# Makefile for hdf2nimrod
# D. Garnier 6/26/01

SHELL = /bin/sh

############# you may need to modify this ##############

# need to force c++ since I have a lame mix of file endings now
# this may need to be g++ on your system
#CC = /usr/bin/g++
CC = /usr/bin/c++

################ below should work without mod ##############

.SUFFIXES:
.SUFFIXES: .cc .c .o

vpath %.h ../includes
vpath %.c ../sources
vpath %.cc ../sources

PDFOBJS = plot_contour.o PDFOutput.o
HDFOBJS = HDFOutput.o

PPOBJS = \
    CDipoleIntStable.o CDipoleStd.o \
    contour.o \
    CPlasmaModel.o \
    fpoly.o \
    GetFluxMoments.o GetFluxParameters.o \
    GetPlasmaParameters.o \
    interpolate.o \
    nrRomberg.o \
    nrSpline.o nrutil.o plasma.o psigrid.o rolldown.o 

#for root only
DEFINES = -DDIPOLE
LIBS += -lmfhdf -ldf -ljpeg -lz 
OBJS = $(PPOBJS)

CFLAGS += -O2 -g -I../includes $(DEFINES)

.c.o .cc.o:
	$(CC) $(CFLAGS) -c $< -o $@

apps: hdf2nimrod contourhdf

hdf2nimrod: $(OBJS) HDF2Nimrod.o HDFInput.o
	$(CC) -o $@ $^ $(LIBS)

contourhdf: $(OBJS) PostProc.o HDFInput.o
	$(CC) -o $@ $^ $(LIBS)

.PHONY : clean
clean:
	-rm *.o




