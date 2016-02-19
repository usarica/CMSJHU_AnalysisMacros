
ROOTCFLAGS    = $(shell root-config --cflags) -Lrootlib -lEG
ROOTGLIBS     = $(shell root-config --glibs) -Lrootlib -lEG

CXX           = g++
CXXFLAGS      = -fPIC -ansi -D_GNU_SOURCE -g -O2

ARCH         := $(shell root-config --arch)

ifeq ($(ARCH),macosx)
# MacOS X with cc (GNU cc 2.95.2 and gcc 3.3)
CXX           = c++ -lm
CXXFLAGS      = -O2 -pipe -Wall -W -Woverloaded-virtual -fPIC -Wno-deprecated -O -ansi -D_GNU_SOURCE  -j 4
# The SOFLAGS will be used to create the .dylib,
# the .so will be created separately
DllSuf        = dylib
UNDEFOPT      = suppress
endif

CXXFLAGS      += $(ROOTCFLAGS) -I$(ROOFITSYS)/include/ -I$(CMSSW_BASE)/src

NGLIBS         = $(ROOTGLIBS)
# Hack here, because RooFit is removed from ROOT:
NGLIBS        += -L $(ROOFITSYS)/lib/ -lMinuit -lRooFitCore -lRooFit #-lboost_regex
NGLIBS  += -L $(CMSSW_BASE)/lib/$(SCRAM_ARCH)/ -lZZMatrixElementMELA
GLIBS          = $(filter-out -lNew, $(NGLIBS))


COMPILEPATH      = $(PWD)/
INCLUDEDIR       = $(COMPILEPATH)src/
OUTLIB           = $(COMPILEPATH)lib/
RUNDIR           = $(COMPILEPATH)


$(OUTLIB)PDGHelpers.o: $(INCLUDEDIR)PDGHelpers.cc
	@mkdir -p $(OUTLIB)
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)PDGHelpers.o $(GLIBS) $<

$(OUTLIB)LHEParticleSmear.o: $(INCLUDEDIR)LHEParticleSmear.cc
	@mkdir -p $(OUTLIB)
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)LHEParticleSmear.o $(GLIBS) $<




$(OUTLIB)Particle.o: $(INCLUDEDIR)Particle.cc $(OUTLIB)PDGHelpers.o
	@mkdir -p $(OUTLIB)
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)Particle.o $(GLIBS) $<

$(OUTLIB)ZZCandidate.o: $(INCLUDEDIR)ZZCandidate.cc $(OUTLIB)Particle.o
	@mkdir -p $(OUTLIB)
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)ZZCandidate.o $(GLIBS) $<

$(OUTLIB)ParticleComparators.o: $(INCLUDEDIR)ParticleComparators.cc $(OUTLIB)Particle.o
	@mkdir -p $(OUTLIB)
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)ParticleComparators.o $(GLIBS) $<

$(OUTLIB)Event.o: $(INCLUDEDIR)Event.cc $(OUTLIB)ZZCandidate.o $(OUTLIB)ParticleComparators.o
	@mkdir -p $(OUTLIB)
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)Event.o $(GLIBS) $<

$(OUTLIB)melaHelpers.o: $(INCLUDEDIR)melaHelpers.cc $(OUTLIB)Event.o
	@mkdir -p $(OUTLIB)
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)melaHelpers.o $(GLIBS) $<



$(OUTLIB)HiggsComparators.o: $(INCLUDEDIR)HiggsComparators.cc $(OUTLIB)Event.o
	@mkdir -p $(OUTLIB)
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)HiggsComparators.o $(GLIBS) $<

$(OUTLIB)OptionParser.o: $(INCLUDEDIR)OptionParser.cc $(OUTLIB)ParticleComparators.o $(OUTLIB)HiggsComparators.o $(OUTLIB)melaHelpers.o
	@mkdir -p $(OUTLIB)
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)OptionParser.o $(GLIBS) $<


$(OUTLIB)BaseTree.o: $(INCLUDEDIR)BaseTree.cc
	@mkdir -p $(OUTLIB)
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)BaseTree.o $(GLIBS) $<

$(OUTLIB)HVVTree.o: $(INCLUDEDIR)HVVTree.cc $(OUTLIB)BaseTree.o $(OUTLIB)OptionParser.o
	@mkdir -p $(OUTLIB)
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)HVVTree.o $(GLIBS) $<

$(OUTLIB)converter.o: $(INCLUDEDIR)converter.cc $(OUTLIB)HVVTree.o $(OUTLIB)LHEParticleSmear.o
	@mkdir -p $(OUTLIB)
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)converter.o $(GLIBS) $<

$(OUTLIB)convertLHE.o: $(INCLUDEDIR)convertLHE.cc $(OUTLIB)converter.o
	@mkdir -p $(OUTLIB)
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)convertLHE.o $(GLIBS) $<

$(OUTLIB)convertPythia.o: $(INCLUDEDIR)convertPythia.cc $(OUTLIB)converter.o
	@mkdir -p $(OUTLIB)
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)convertPythia.o $(GLIBS) $<

$(OUTLIB)Reader.o: $(INCLUDEDIR)Reader.cc $(OUTLIB)converter.o
	@mkdir -p $(OUTLIB)
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)Reader.o $(GLIBS) $<

$(OUTLIB)main.o: $(INCLUDEDIR)main.cc $(OUTLIB)Reader.o $(OUTLIB)convertPythia.o $(OUTLIB)convertLHE.o
	@mkdir -p $(OUTLIB)
	$(CXX) $(CXXFLAGS) -c -o $(OUTLIB)main.o $(GLIBS) $<


dolinks:
	ln -sf $(COMPILEPATH)macros/loadLib.cc $(RUNDIR)loadLib.cc
	ln -sf $(COMPILEPATH)macros/trimPythia.cc $(RUNDIR)trimPythia.cc
	ln -sf $(COMPILEPATH)macros/loadLib.py $(RUNDIR)loadLib.py
	ln -sf $(COMPILEPATH)macros/trimPythia.py $(RUNDIR)trimPythia.py

analyzer: $(OUTLIB)main.o
	$(CXX) $(CXXFLAGS) -o $(RUNDIR)analyzer $(OUTLIB)*.o $(GLIBS)

all: analyzer dolinks

clean:
	rm -rf $(OUTLIB)
	rm -rf $(RUNDIR)Pdfdata
	rm -f $(RUNDIR)input.DAT
	rm -f $(RUNDIR)process.DAT
	rm -f $(RUNDIR)ffwarn.dat
	rm -f $(RUNDIR)br.sm*
	rm -f $(RUNDIR)*.cc
	rm -f $(RUNDIR)analyze*
	rm -f $(COMPILEPATH)*.o
