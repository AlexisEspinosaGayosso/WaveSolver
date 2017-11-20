# -*- indent-tabs-mode:t; -*-
#*****************************
#Makefile for 3dWaveFDTD
#*****************************

#Directories --------------------
IDIR     = ../include
ODIR     = obj
LDIR     = ../lib
BDIR		= ../bin

#Name of the executables----------
_EXEOUT	= 3dWaveFDTD
EXEOUT	= $(patsubst %,$(BDIR)/%,$(_EXEOUT))

_EXETEST	= 3dWaveFDTD_Test
EXETEST	= $(patsubst %,$(BDIR)/%,$(_EXETEST))

#Dependencies--------------------
_DEPS    = allvars.h proto.h
DEPS     = $(patsubst %,$(IDIR)/%,$(_DEPS))

#Object file list ---------------
_OBJ     = domain.o fdtd.o communicate.o setup.o dynamicArrays.o main.o
OBJ      = $(patsubst %,$(ODIR)/%,$(_OBJ))

_OBJTEST	= domain.o fdtd.o communication_test.o setup.o main.o
OBJTEST  = $(patsubst %,$(ODIR)/%,$(_OBJTEST))

#Compiling variables-------------
CC        = mpicc
OPTIMIZE = -O3
WARNING   = -Wall
# OPT       = -DOPENMP
OPT      += -DMPI
OPTIONS   = $(OPTIMIZE) $(OPT)
CFLAGS    = $(OPTIONS) -I $(IDIR)
LIBS      = -fopenmp -lm

#------------------------------------------------------------------
#------------------------------------------------------------------

#Preparation Rules ---------------
$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

#Main Rules ----------------------
3dWaveFDTD: $(OBJ)
	$(CC) -o $(EXEOUT) $^ $(CFLAGS) $(LIBS)

test: $(OBJTEST)
	$(CC) -o $(EXETEST) $^ $(CFLAGS) $(LIBS)


#Cleaning Rules ------------------
.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
