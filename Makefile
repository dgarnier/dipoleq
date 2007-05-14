# Makefile for dipoleq v 0.9
# D. Garnier 6/13/01

SHELL = /bin/sh

############# you may need to modify this ##############

# need to force c++ since I have a lame mix of file endings now
# this may need to be g++ on your system
#CC = /usr/bin/g++
CC = /usr/bin/c++

# in case you don't have the HDF v4 installed, comment out #
USEHDF=1
# in case you don't have ClibPDF from www.fastio.com comment out #
USEPDF=1

################ below should work without mod ##############

.SUFFIXES:
.SUFFIXES: .cc .c .o

vpath %.h includes
vpath %.c sources
vpath %.cc sources

PDFOBJS = plot_contour.o PDFOutput.o
HDFOBJS = HDFOutput.o
NIMRODOBJS = NimrodOutput.o

ROOTOBJS = AddCoilJ.o AddShellJ.o CDipoleIntStable.o CDipoleStd.o \
    CPlasmaModel.o DelChiSqr.o FileInput.o \
    FileOutput.o FindJ.o FindMeasFit.o \
    Find_ShellCurrent.o Find_dJdy.o GetFluxMoments.o GetFluxParameters.o \
    GetPlasmaParameters.o InitJ.o J_DipoleStd.o J_IsoNoFlow.o \
    J_Std.o LeastSquares.o LoadBndryGreens.o \
    LoadMeasGreens.o LoadShellGreens.o \
    PlasmaBoundary.o PsiBoundary.o Restart.o SVDFit.o bcucof.o \
    bcuint.o coil.o contour.o dUnkn.o fpoly.o green.o interpolate.o \
    limiter.o ludcmp.o meas_J0.o meas_bp.o meas_bpangle.o meas_circle.o \
    meas_coilcur.o meas_flux.o meas_mag_Xxx.o \
    meas_plasmacur.o meas_pnorm.o \
    meas_ppsix.o meas_press.o meas_saddle.o \
    measurement.o multitask.o nrRomberg.o \
    nrSpline.o nrutil.o plasma.o polygon.o psigrid.o rolldown.o \
    separatrix.o shell.o stdio_dmatrix.o tokamak.o tokgreen.o

#for root only
DEFINES = -DDIPOLE
LIBS = 
OBJS = $(ROOTOBJS)

#for HDF
ifdef USEHDF
DEFINES += -DHDFOUTPUT
LIBS    += -lmfhdf -ldf -ljpeg -lz
OBJS	+= $(HDFOBJS)
endif

#for PDF
ifdef USEPDF
DEFINES += -DPDFOUTPUT
LIBS    += -lcpdfm
OBJS	+= $(PDFOBJS)
endif

#for Nimrod output
ifdef USENIMROD
DEFINES += -DNIMROD_OUTPUT
OBJS	+= $(NIMRODOBJS)
endif


CFLAGS += -O2 -g -Iincludes $(DEFINES)

.c.o .cc.o:
	$(CC) $(CFLAGS) -c $< -o $@

dipoleq: $(OBJS) SimDipEq.o
	$(CC) -o $@ $^ $(LIBS)

install: dipoleq
	if [ ! -d $(HOME)/bin ]; then \
		mkdir $(HOME)/bin ; \
	fi; 
	install -m 755 dipoleq $(HOME)/bin

.PHONY : clean
clean:
	-rm *.o




