#!/bin/bash

cd ../AnalysisTutorial

if [ ! -d "build" ]
then
        mkdir build
fi

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh

cd build
asetup AnalysisBase,25.2.29

if [ ! -d "CMakeFiles" ]
then
        cmake ../source
        make -j
else
        make -j
fi
source x86_64-*/setup.sh

cd ../../run
echo "$PWD"

if [ ! -d "nTuples" ]
then
        mkdir nTuples
fi
cd nTuples

echo "$PWD"

rm -rf submitDir
echo "ln -sf ../derivation/DAOD_PHYS.mcval.pool.100000.root"
ln -sf ../derivation/DAOD_PHYS.mcval.pool.100000.root

eval "python ../../AnalysisTutorial/source/MyAnalysis/share/ATestRun_eljob.py -c ../../AnalysisTutorial/source/MyAnalysis/data/config.yaml -i /eos/user/t/toheintz/ultraSlowLLP/b2f_simulation/run/derivation/DAOD_PHYS.mcval.pool.100000.root"
#eval "python ../../l1-calo-llp-validation/source/MCVal/share/MCVal_eljob.py -c ../../source/config.yaml -i  /eos/user/t/toheintz/ultraSlowLLP/b2f_simulation/run/derivation/DAOD_PHYS.mcval.pool.100000.root"
#eval "python ../../l1-calo-llp-validation/source/MCVal/share/MCVal_eljob.py -c ../../source/config.yaml -i  DAOD_PHYS.mcval.pool.100000.root"

