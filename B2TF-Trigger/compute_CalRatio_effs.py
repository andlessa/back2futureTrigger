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

c = 3e8

FORMAT = '%(levelname)s: %(message)s at %(asctime)s'
logging.basicConfig(format=FORMAT,datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger()   


 
def getEffForCalRatio(tree) -> dict:


    evt_cutFlow = {}
    evt_cutFlow['All'] = 1
    

    tracksDelayed = None
    tracksOnTime = None
    if hasattr(tree,'AllTracks'):
        tracksDelayed = []
        tracksOnTime = []
        for track in tree.AllTracks:
            if track.PT < 2.0:
                continue
            l_track = np.sqrt(track.X**2 + track.Y**2 + track.Z**2)
            t_readout =  track.T - l_track*1e-3/c
            if (25e-9 < t_readout < 35e-9):               
                tracksDelayed.append(track)
            elif (t_readout < 10e-9):
                tracksOnTime.append(track)
    
    
    eventOnTime = {
                    'jetsL1' : tree.L1SmallJetOnTime, 
                    'jetsHLT' : tree.HLTJetOnTime, 
                    'tracks' : tracksOnTime
                  }  
                   
    eventDelayed = {
                    'jetsL1' : tree.L1SmallJetDelayed, 
                    'jetsHLT' : tree.HLTJetDelayed, 
                    'tracks' : tracksDelayed
                  }                   
    

    ### Check if the N-1 or N event records pass the triggers
    cutFlows = []
    for event,label in [(eventOnTime,'N-1'),(eventDelayed,'N')]:
        # Check low ET trigger
        l1_cutflow = CalRatioLowET_L1(event['jetsL1'])
        # If failed, check high ET trigger
        if any(val == 0 for val in l1_cutflow.values()):
             l1_cutflow = CalRatioHighET_L1(event['jetsL1'],pTmin = 60.0)
        
        l1_pass = all(val != 0 for val in l1_cutflow.values())
        
        hlt_cutflow = CalRatio_HLT(event['jetsHLT'],event['tracks'])

        hlt_pass = all(val != 0 for val in hlt_cutflow.values())
        
        cutFlow = {f'({label}) {key}' : val for key,val in l1_cutflow.items()}
        cutFlow.update({f'({label}) {key}' : val*l1_pass 
                        for key, val in hlt_cutflow.items()})
        cutFlows.append(cutFlow)
            
        # If the event passed both triggers, skip next event
        if l1_pass and hlt_pass:
            break
    
    # Sort cutFlow by last cut and total number of cuts passed
    cutFlows = sorted(cutFlows, 
                      key=lambda cDict: (list(cDict.values())[-1],sum(cDict.values())),
                      reverse=True)
    # return the best cutflow
    evt_cutFlow.update(cutFlows[0])

    return evt_cutFlow

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
    effsDict = {"Trigger" : 0.0,
                "Nevents" : 0.0}
    
    evts_cutFlow = {}
    
    
    nevts = tree.GetEntries()
    for ievt in tqdm.tqdm(range(nevts),position=ijob,
                          desc=os.path.basename(inputFile),
                          leave=True):
        tree.GetEntry(ievt) 
        effsDict['Nevents'] += 1.0
        evt_cutFlow = getEffForCalRatio(tree)
        if all(val != 0 for val in evt_cutFlow.values()):
            evt_eff = 1.0
        else:
            evt_eff = 0.0

        # Add event efficiency to total efficiencies 
        # for the given selection (for each tau value)
        effsDict['Trigger'] += evt_eff
        for key,val in evt_cutFlow.items():
            if key not in evts_cutFlow:
                evts_cutFlow[key] = 0
            evts_cutFlow[key] += val
        
    f.Close()

    # Sort keys:
    sorted_keys = sorted(evts_cutFlow.keys(), key=lambda k: ('All' in k, '(N-1)' in k, 'L1' in k, evts_cutFlow[k]),
                         reverse=True)
    evts_cutFlow_sorted = {key : evts_cutFlow[key] for key in sorted_keys}
    evts_cutFlow = evts_cutFlow_sorted
    
    # Divide the total efficiency by the number of events:
    effsDict["TriggerErr"] = np.sqrt(effsDict["Trigger"])/effsDict['Nevents']
    effsDict["Trigger"] = effsDict["Trigger"]/effsDict['Nevents']

    # Store input file name
    effsDict['inputFile'] = inputFile
    # Get model information
    modelDict = getModelDict(inputFile,verbose=False)
    effsDict.update(modelDict)
    evts_cutFlow.update(modelDict)

    return effsDict,evts_cutFlow

def saveOutput(effsDict,outputFile,
               evts_cutFlow=None,cutFile=None):
        
     # Extract info for comments
    nevts = effsDict.pop('Nevents')
    inputFile = effsDict.pop('inputFile')

    
    effsDict['eff'] = effsDict.pop('Trigger')
    effsDict['effError'] = effsDict.pop('TriggerErr')
    
    # Get column labels and data
    columnsLabels = list(effsDict.keys())
    data = np.array([list(effsDict.values())])
    header = ';'.join(columnsLabels)
    np.savetxt(outputFile, data, 
                header=f'Input file: {inputFile}\nNumber of events: {nevts}\n{header}',
                delimiter=';',fmt='%1.6e')

    # Save cutflow (if defined)
    # Get column labels and data
    if evts_cutFlow and cutFile:
        columnsLabels = list(evts_cutFlow.keys())
        data = np.array([list(evts_cutFlow.values())])
        header = ';'.join(columnsLabels)
        np.savetxt(cutFile, data,
                header=f'Input file: {inputFile}\nNumber of events: {nevts}\n{header}',
                delimiter=';',fmt='%1.6e')

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
        effsDict,evts_cutFlow = p.get()
        outFile = effsDict['inputFile'].split('.root')[0].split('.hepmc')[0]
        effFile = outFile + output_suffix +'_calRatio_effs.csv'
        cutFile = outFile + output_suffix +'_calRatio_cutFlow.csv'
        cutFlowStr += f'CutFlow for {effsDict['inputFile']}:\n'
        for key,val in evts_cutFlow.items():
            cutFlowStr += f'  {key} : {val/evts_cutFlow['All']:1.3e}\n'
        saveOutput(effsDict,effFile,evts_cutFlow,cutFile)        

    logger.debug(cutFlowStr+'\n')
            
        
    logger.info("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
