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
       asetup Athena,24.0,latest
       cmake -DATLAS_PACKAGE_FILTER_FILE=../package_filters.txt ../../../athena/Projects/WorkDir
       make -j
else
       asetup Athena,24.0,latest
       make -j
fi

source x86_64-el9-gcc13-opt/setup.sh

cd ../run
if [ -d "sim" ]
then
        rm -r sim
fi
mkdir sim

#cp ../source/evgen.100000.root .
cd sim

Sim_tf.py \
  --inputEVNTFile   "../gen/evgen.100000.root" \
  --outputHITSFile  "hits.100000.root" \
  --AMIConfig s4369 \
  --AMITag s4369 \
  --outputEVNT_TRFile "StoppedPartPositions.pool.root" \
  --trackRecordType "stopped" \
  --maxEvents 2 \
  --preExec 'flags.Input.SpecialConfiguration={"postInclude":"SlowLLPs.SlowLLPsConfig.BTFT_Cfg","BTFTChi1Mass":"248700", "BTFTChi1Lifetime":"3000", "BTFTChi0Mass":"173700"}' \
  --preInclude 'Campaigns.MC23eSimulationMultipleIoV,G4DebuggingTools.DebugSlowLLPs'

#--preExec "flags.Exec.DebugMessageComponents=['SlowLLPs']" \
  #--outputEVNT_TRFile "StoppedPartPositions.pool.root" \
  #--trackRecordType "stopped" \
  #--maxEvents 2 \


  #--preInclude "G4DebuggingTools.DebugSlowLLPs,SlowLLPs.SlowLLPsConfig.BTFT_Cfg" 

  #--postInclude "SlowLLPs.SlowLLPsConfig.BTFT_Cfg" \
