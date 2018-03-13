# ***************************************************************************
# *                                                                         *
# *           Aflow STEFANO CURTAROLO - Duke University 2003-2013           *
# *                                                                         *
# ***************************************************************************
# MAKEFILE FOR AFLOW

CPP=g++
# -O3
#CPP=g++-4.7 -std=c++0x
# -std=gnu++0x
#CCFLAGS= -DCYGWIN -Wall -W -fno-inline
CCFLAGS= -Wall -W -fno-inline
OPTS=
#-O1
#-O3
# -march=opteron
VERSNUMBER=30864
DATEMAX=20141231

VERS=-DVERSION=\"$(VERSNUMBER)\"
DATES=-DDATEMAX=$(DATEMAX)

DAY=$(shell date +%y)$(shell date +%m)$(shell date +%d)
#TIME="20$(shell date +%y)-$(shell date +%m)-$(shell date +%d) $(shell date +%H):$(shell date +%M):$(shell date +%S)"
TIME="20$(shell date +%y)-$(shell date +%m)-$(shell date +%d)"
LIBS=-lpthread
# -static
# -lglut -lGLU -lGL -D_HOST_=\"`hostname`\"
#VERSION := $(shell cat aflow.cpp | /home/auro/bin/getvalue VERSION)-$(TIME)
#VERSIONM := $(shell cat Makefile | /home/auro/bin/getvalue MAKEFILE)-$(TIME)
#DEBUG=-DXSTR_DEBUG
DEBUG=


AUROSTD_SRC = aurostd.cpp
AUROSTD_DEPENDENCIES_H = aurostd.h aurostd_argv.h aurostd_xcomplex.h aurostd_xmatrix.h aurostd_xoption.h aurostd_xrandom.h aurostd_xscalar.h aurostd_xtensor.h aurostd_xvector.h aurostd_xcombos.h
AUROSTD_DEPENDENCIES_CPP = aurostd.cpp aurostd_main.cpp aurostd_argv.cpp aurostd_xoption.cpp aurostd_boot.cpp aurostd_crc64.cpp aurostd_xscalar.cpp aurostd_xcomplex.cpp aurostd_xvector.cpp aurostd_xmatrix.cpp aurostd_xtensor.cpp aurostd_xcombos.cpp
AUROSTD_OBJ = aurostd.o

clean:
	rm -rf *~ *.o *.rpo *~ *.o *.rpo #*

