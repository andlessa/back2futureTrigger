#!/bin/bash

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh

if [ ! -d "../build" ]
then
       mkdir ../build
fi

cd ../build

if [ ! -d "CMakeFiles" ]
then
       asetup Athena,24.0.72
       cmake -DATLAS_PACKAGE_FILTER_FILE=../package_filters.txt ../../../athena/Projects/WorkDir
       make -j
else
       asetup Athena,24.0.72
       make -j
fi

source x86_64-el9-gcc13-opt/setup.sh

cd ../run
if [ -d "derivation" ]
then
        rm -r derivation
fi
mkdir derivation
cd derivation

infile=../reco/100000.AOD.root
outfile="mcval.pool.1000000.root"


Derivation_tf.py \
  --CA \
  --inputAODFile ../reco/100000.AOD.root \
  --outputDAODFile mcval.pool.100000.root \
  --formats PHYS \

