#
# Makefile for GENIE test programs
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# May 02, 2004
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES  := $(LHAPDF_LIBRARIES) $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)
#LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)


TGT =	AnalyseGENIE 		 
 	
all: $(TGT)

AnalyseGENIE: FORCE
	$(CXX) $(CXXFLAGS) -c AnalyseGENIE.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) AnalyseGENIE.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/AnalyseGENIE

#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/AnalyseGENIE

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/AnalyseGENIE

FORCE:

# DO NOT DELETE
