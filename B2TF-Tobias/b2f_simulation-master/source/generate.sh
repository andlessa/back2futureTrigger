#!/bin/bash

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh

if [ ! -d "../build" ]
then
       mkdir ../build
fi

cd ../build

#if [ ! -d "CMakeFiles" ]
#then
#       asetup Athena,24.0.72
#       cmake -DATLAS_PACKAGE_FILTER_FILE=../package_filters.txt ../../../athena/Projects/WorkDir
#       make -j
#else
#       asetup Athena,24.0.72
#       make -j
#fi

asetup AthGeneration,23.6.36

#source x86_64-el9-gcc13-opt/setup.sh

cd ../run
if [ -d "gen" ]
then
        rm -r gen
fi
mkdir gen

cd gen
mkdir 100000
cd 100000

cp ../../../source/mc_Pythia8_A14_NNPDF23LO_slowLLP.py .

outfile="../evgen.100000.root"
Gen_tf.py --ecmEnergy=13600. --maxEvents=200  --randomSeed=123456 --outputEVNTFile=$outfile --jobConfig=$PWD

