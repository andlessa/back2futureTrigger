#!/usr/bin/env python3

import os, time, sys
sys.path.append('../')
import numpy as np
import pyhepmc
import random
import logging
import subprocess
import multiprocessing
import tqdm
from helper import getModelDict
from triggers_helper import CalRatioLowET_L1,CalRatioHighET_L1,CalRatio_HLT
delphesDir = os.path.abspath("../DelphesLLP")
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(delphesDir,"external")

import ROOT


ROOT.gSystem.Load(os.path.join(delphesDir,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')

c_light = 2.99792458e8

FORMAT = '%(levelname)s: %(message)s at %(asctime)s'
logging.basicConfig(format=FORMAT,datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger()   


def getEfficiencies(inputFile,ijob=0):


    # Set the seed for each run, so the results independ on
    # how many files are run or their ordering
    random.seed(123)

    # Load hepmc File
    if not os.path.isfile(inputFile):
        logger.error(f"File {inputFile} not found")
        return None
    
    try:
        f = ROOT.TFile(inputFile,'read')
        tree = f.Get("Delphes")
    except:
        logger.error(f"Error reading {inputFile}")
        return None
    
    # Total efficiencies (for each tau value)
    cutFlow = {"Nevents" : 0.0,
               "DV(N-1) > 0": 0.0,
               "MET(N-1) > 30 GeV" : 0.0,
               "dPhi(MET,DV)(N-1) < 1.2" : 0.0
               }
    
    
    nevts = tree.GetEntries()
    for ievt in tqdm.tqdm(range(nevts),position=ijob,
                          desc=os.path.basename(inputFile),
                          leave=True):
        tree.GetEntry(ievt) 
        cutFlow['Nevents'] += 1.0

        llpDec = {}
        llpBoost = {}        
        llps = list(tree.llpParticles)
        if len(llps) != 1:
            logger.error(f"{len(llps)} LLPs found!")
            return 0.0
        for illp,llp in enumerate(llps):
            p = np.array([llp.Px,llp.Py,llp.Pz,llp.E])
            beta = np.linalg.norm(p[:3])/p[3]
            llpBoost[illp] = beta
        
        for d in tree.llpDirectDaughters:
            illp = d.M1
            if illp in llpDec:
                continue
            beta = llpBoost[illp]
            x = np.array([d.X,d.Y,d.Z,d.T,beta])
            l = np.linalg.norm(x[:3])*1e-3
            t_readout = x[3]-l/c_light
            t_readout = t_readout*1e9
            x[3] = t_readout
            llpDec[illp] = x

        # Require decay in N-1 and DV in ~4m-8m
        good_DVs = []
        for illp,dec in llpDec.items():
            llp = llps[illp]
            X,Y,Z,t_readout,beta = dec
            if not (0.0 < t_readout < 15.0):
                continue
            l = np.linalg.norm([X,Y,Z])*1e-3
            lxy = np.linalg.norm([X,Y])*1e-3
            lz = np.abs(Z)*1e-3
            if lz < 5.0: # barrel
                if not (abs(llp.Eta) < 0.7):
                    continue
                if not (3.0 < lxy < 8.0):
                    continue
            elif lz < 15.0: # endcap
                if not (1.3 < abs(llp.Eta) < 2.5):
                    continue
                if not (lxy < 10.0):
                    continue
            else:
                continue

            good_DVs.append(illp)
        
        
        if not good_DVs:
            continue

        cutFlow['DV(N-1) > 0'] += 1.0
        llp = llps[good_DVs[0]]
        

        metOnTime = tree.L1METOnTime.At(0)
        if metOnTime.MET < 30.0:
            continue
        cutFlow['MET(N-1) > 30 GeV'] += 1.0


        dphi = np.abs(llp.Phi-metOnTime.Phi)
        if dphi > np.pi:
            dphi = 2*np.pi - dphi
        if dphi > 1.2:
            continue
        cutFlow['dPhi(MET,DV)(N-1) < 1.2'] += 1.0

    modelDict = getModelDict(inputFile,verbose=False)
    for key,val in modelDict.items():
        print(f"{key} = {val}")

    nevts = cutFlow['Nevents']
    for key,evts in cutFlow.items():
        print(f"{key} : {int(evts):d} ({evts/nevts:1.3e})")


    # Store input file name
    cutFlow['inputFile'] = inputFile
    # Get model information
    
    cutFlow.update(modelDict)
    
    return cutFlow


if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Compute the efficiencies for a given input ROOT file" )
    ap.add_argument('-f', '--inputfile',nargs='+',
            help='path to the ROOT file(s).')
    ap.add_argument('-n', '--ncpus',type=int,default=1,
            help='number of parallel jobs to run.')
    ap.add_argument('-o', '--output_suffix',default='',type=str,
            help='suffix to be added to the output files.')
    ap.add_argument('-v', '--verbose', default='info',
            help='verbose level (debug, info, warning or error). Default is info')
    

    args = ap.parse_args()

    level = args.verbose
    levels = { "debug": logging.DEBUG, "info": logging.INFO,
               "warn": logging.WARNING,
               "warning": logging.WARNING, "error": logging.ERROR }
    if level in levels:       
        logger.setLevel(level = levels[level])



    # First make sure the correct env variables have been set:
    LDPATH = subprocess.check_output('echo $LD_LIBRARY_PATH',shell=True,text=True)
    ROOTINC = subprocess.check_output('echo $ROOT_INCLUDE_PATH',shell=True,text=True)
    pythiaDir = os.path.abspath('../MG5/HEPTools/pythia8/lib')
    delphesDir = os.path.abspath('../DelphesLLP/external')
    if pythiaDir not in LDPATH or delphesDir not in ROOTINC:
        print('\n\nEnviroment variables not properly set. Run source setenv.sh first.\n\n')
        sys.exit()


    t0 = time.time()



    # Split input files by distinct models and get recast data for
    # the set of files from the same model:
    ncpus = args.ncpus
    inputFileList = args.inputfile
    output_suffix = args.output_suffix
    ncpus = min(ncpus,len(inputFileList))
    pool = multiprocessing.Pool(processes=ncpus)
    children = []
    for ijob,inputfile in enumerate(inputFileList):
        p = pool.apply_async(getEfficiencies, args=(inputfile,ijob))
        children.append(p)

    nfiles = len(inputFileList)
    
    
    cutFlowStr=''
    logger.info(f'Analyzing {nfiles} files')
    for p in children: 
        cutFlow = p.get()

        
    logger.info("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
