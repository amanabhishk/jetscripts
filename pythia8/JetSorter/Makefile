################################################################
## Makefile for building Pythia-ROOT-Fastjet code (jetsorter) ## 
## Aman Abhishek, 23.03.2015                                 ##  
################################################################

SHELL = /bin/bash
CXXFLAGS := -O2 -ansi -pedantic -W -Wall -Wshadow -Wno-shadow -m64

ROOT=$(shell root-config --cflags --glibs)
PYTHIA=$(shell pythia8-config --cflags --libs)
FASTJET=$(shell fastjet-config --cxxflags --libs)

PYTHIALIB=/usr/local/lib/libpythia8.so
QCDLIB=/home/aabhishe/fastjet-3.1.0/plugins/QCDaware/QCDAware.o

.SECONDEXPANSION:
.PHONY: all clean

all:
	@echo "No make for all objects provided."

# Event generation files for 3 samples.
Zjet_eventgen: $$@.cc $(PYTHIALIB)
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(PYTHIA) $(FASTJET) -ldl -o $@.exe 
gammajet_eventgen: $$@.cc $(PYTHIALIB)
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(PYTHIA) $(FASTJET) -ldl -o $@.exe 
dijet_eventgen: $$@.cc $(PYTHIALIB)
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(PYTHIA) $(FASTJET) -ldl -o $@.exe 

#Physics-definition algorithm
physics_def: $$@.cc $(PYTHIALIB)
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(PYTHIA) $(FASTJET) -ldl -o $@.exe 

#Hadron-based definition
hadronic_def: $$@.cc $(PYTHIALIB)
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(PYTHIA) $(FASTJET) -ldl -o $@.exe 

#QCD-aware definition
QCDaware_def: $$@.cc $(PYTHIALIB) $(QCDLIB)
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(PYTHIA) $(FASTJET) -ldl -o $@.exe 

#Feynman generation
feynman_dijet: $$@.cc $(PYTHIALIB)
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(PYTHIA) $(FASTJET) -ldl -o $@.exe 
feynman_Zgamma: $$@.cc $(PYTHIALIB)
	$(CXX) $^ -w $(CXXFLAGS) $(ROOT) $(PYTHIA) $(FASTJET) -ldl -o $@.exe 


# Clean.
clean:
	rm -f *.log *.exe *.d *.so
