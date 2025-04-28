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

export ATHENA_CORE_NUMBER=8

cd ../run
if [ -d "reco" ]
then
        rm -r reco
fi
mkdir reco
cd reco


Reco_tf.py --inputHITSFile /eos/user/t/toheintz/ultraSlowLLP/b2f_simulation/run/sim/hits.100000.root --outputAODFile 100000.AOD.root --maxEvents 200 --AMIConfig r16210 --AMITag r16210
