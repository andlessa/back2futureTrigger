# Builds DelphesPythia8Matched against an existing Delphes build (libDelphes.so)
# and the Pythia8 installation used by MadGraph.
#
#   make                       (uses the defaults below)
#   make DELPHES=/path/to/Delphes PYTHIA8=/path/to/pythia8

TOP     := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))..)
DELPHES ?= $(TOP)/DelphesLLP
PYTHIA8 ?= $(TOP)/MG5/HEPTools/pythia8

CXX      ?= g++
CXXFLAGS := -O2 $(shell root-config --cflags) \
            -I$(DELPHES) -I$(DELPHES)/external -I$(PYTHIA8)/include
LDFLAGS  := -L$(DELPHES) -lDelphes -Wl,-rpath,$(DELPHES) \
            -L$(PYTHIA8)/lib -lpythia8 -Wl,-rpath,$(PYTHIA8)/lib \
            $(shell root-config --libs) -lEG -ldl -lz

DelphesPythia8Matched: DelphesPythia8Matched.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)

clean:
	rm -f DelphesPythia8Matched

.PHONY: clean
