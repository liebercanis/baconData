# gnu makefile for nasdaq/devel 
# builds library for nasdaq and makes executable readit
########################################################
# see http://www.gnu.org/software/make/manual/
# Here is a summary of the automatic variables.
# $@    The file name of the target.
# $%    The target member name, when the target is an archive member.
# $<    The name of the first prerequisite.
# $?    The names of all the prerequisites that are newer than the target
# $^    The names of all the prerequisites, omits duplicates
# $+                                      , retains duplicates
# $* The stem with which an implicit rule matches
########################################################
ARCH          = $(shell root-config --arch)
LOCAL        := $(shell pwd)
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
BOBJ=$(HOME)/bacon/baconData/bobj
# Linux with egcs, gcc 2.9x, gcc 3.x (>= RedHat 5.2)
CXX           = g++ -ggdb 
## -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS         += $(ROOTLIBS) 
CXXFLAGS     += -I$(LOCAL) -I/$/usr/include -I/usr/local/include -I$(BOBJ)

CC=gcc -g -Wall -c $(INC_PATH) 
CPP= $(CXX) $(CXXFLAGS) $(INC_PATH) 


CC+=-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
CPP+=-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64

LIBS+=  -L$(FFTW)/lib -lfftw3 $(ROOTLIBS) -L$($BOBJ) -lBaconData

all: convertRaw

convertRaw: convertRaw.cc
	$(CPP) -o ../bin/$@ $<  $(LIBS)

.PHONY : clean show

clean:
	rm -f *.o ../bin/readit 
show:
	@echo include path = $(INC_PATH)
	@echo CXXFLAGS = $(CXXFLAGS)
	@echo LIBS = $(LIBS)
	@echo cpp = $(CPP)
