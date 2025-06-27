#!/usr/bin/env python3

import os, time, sys
sys.path.append('../')
import numpy as np
import pyhepmc
import random
import logging
import configparser
import subprocess
import multiprocessing
import itertools
import tqdm
from helper import getModelDict
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


 

def getDataFrom(llpParticles,llpDaughters,
                llp_PDGs=[4000023],invisible_PDGs=[4000022]):
    
    # Consistency checks
    if any(abs(p.PID) not in llp_PDGs for p in llpParticles):
        erroMsg = "LLP PDG from input file does not match definitions"
        logger.error(erroMsg)
        raise ValueError(erroMsg)

    # Build list of LLPs, visible and invisible LLP daughters as pyhepmc GenParticles
    llps = []
    invisibleLists = [[] for i in range(len(llpParticles))]
    visibleLists = [[] for i in range(len(llpParticles))]
    for illp,llp in enumerate(llpParticles):
        llp_momentum = pyhepmc.FourVector(llp.Px,llp.Py,llp.Pz,llp.E)
        llps.append(pyhepmc.GenParticle(momentum = llp_momentum, pid = llp.PID))
        # Select set of daughters from the current llp particle
        llp_daughters = [d for d in llpDaughters if d.M1 == illp]
        visible_particles = []
        invisible_particles = []
        for d in llp_daughters:
            mom = pyhepmc.FourVector(d.Px,d.Py,d.Pz,d.E)            
            p = pyhepmc.GenParticle(momentum = mom, pid = d.PID)
            p.generated_mass = d.Mass
            if abs(p.pid) in invisible_PDGs:
                invisible_particles.append(p)
            else:
                visible_particles.append(p)
        # Only accept single decays to Higgs
        if len(visible_particles) == 1 and visible_particles[0].pid != 25:
            erroMsg = f"Single visible particle from LLP decay with {len(llp_daughters)} daughters (can only handle Higgs)"
            logger.error(erroMsg)
            raise ValueError(erroMsg)
        elif len(visible_particles) > 2:
            errorMsg = f"{len(visible_particles)} visible particles found in LLP decay (can only handle up to 2 particles)"
            logger.error(errorMsg)
            raise ValueError(errorMsg)
       
        # Store particles
        invisibleLists[illp] = invisible_particles
        visibleLists[illp] = visible_particles

        

    eventDict = {}
    for illp,llp in enumerate(llpParticles):
        p_visible = pyhepmc.FourVector(0.,0.,0.,0.)
        p_invisible = pyhepmc.FourVector(0.,0.,0.,0.)
        visible_pdgs = []
        for p in visibleLists[illp]:
            p_visible = p_visible + p.momentum
            visible_pdgs.append(p.pid)
        invisible_pdgs = []
        for p in invisibleLists[illp]:
            p_invisible = p_invisible + p.momentum
            invisible_pdgs.append(p.pid)
    
        visParticle = pyhepmc.GenParticle(momentum = p_visible, pid = illp)
        invisParticle = pyhepmc.GenParticle(momentum = p_invisible, pid = illp)
        llp_momentum = pyhepmc.FourVector(llp.Px,llp.Py,llp.Pz,llp.E)
        llpParticle = pyhepmc.GenParticle(momentum = llp_momentum, pid = llp.PID)
        eventDict[illp] = { 'visible' : visParticle, 
                            'invisible' : invisParticle,
                            'parent' : llpParticle,
                            'visiblePDGs' : visible_pdgs,
                            'invisiblePDGs' : invisible_pdgs,
                          }

    return eventDict

def getJetsFrom(treeJets):
    """
    Convert jets from tree to pyhepmc particles
    """

    # Convert jets to pyhepmc particle
    jets_hepmc = []
    for j in sorted(treeJets, key = lambda x: x.PT, reverse=True):
        j_e = np.sqrt(j.Mass**2 + (j.PT*np.cosh(j.Eta))**2)
        j_px = j.PT*np.cos(j.Phi)
        j_py = j.PT*np.sin(j.Phi)
        j_pz = j.PT*np.sinh(j.Eta)
        p = pyhepmc.FourVector(j_px,j_py,j_pz,j_e)
        jets_hepmc.append(pyhepmc.GenParticle(momentum = p))

    return jets_hepmc

def getJetET(j):

    ET = j.momentum.pt()*(j.momentum.e/j.momentum.p3mod())

    return ET

def defineCutFlow():
    
    cutFlowKeys = [
                    'All',
                    'L1: 40 GeV < MET(N-1) < 100 GeV',
                    'L1: 40 GeV < PT Jet1(N)',
                    'L1: DPhi(Jet(N),MET(N-1)) < 1.0',
                    'HLT: Eta Jet < 2.5 PT Jet > 20 GeV',
                    'HLT: Jet(N) EMF < 0.06',
                    'HLT: DR(Tracks(N),Jet(N)) > 0.2'
                ]
    
    return cutFlowKeys

def getEffFor(tree,llps=[4000023],invisibles=[4000022]):

    evt_eff = 0.0
    cutFlowKeys = defineCutFlow()
    evt_cutFlow = {key : 0.0 for key in cutFlowKeys}

    # Extract necessary data from event
    eventDict = getDataFrom(tree.llpParticles,
                                   tree.llpDirectDaughters,
                                    llp_PDGs=llps,invisible_PDGs=invisibles)
    
    if (len(eventDict) != 2):
        errorMsg = f"{len(eventDict)} LLP found (can only handle 2)"
        logger.error(errorMsg)
        raise ValueError(errorMsg)
        # return evt_eff,evt_cutFlow # return zero effs
    
    # Get MET from the L1 (on-time) calorimenter
    metOnTime = tree.L1METOnTime.At(0)

    # Get jets from the L1 (delayed) calorimeter
    jetsDelayedL1 = list(tree.L1JetDelayed)
        

    evt_cutFlow['All'] += 1

    ### Pre-Selection


    ### L1 Trigger
    if not (40.0 < metOnTime.MET < 100.0):
        return evt_eff,evt_cutFlow

    evt_cutFlow['L1: 40 GeV < MET(N-1) < 100 GeV'] += 1.0

    # Keep only jets passing the pt cut
    jetsDelayedL1 = [j for j in jetsDelayedL1[:] if j.PT > 20.0]
    # Keep only jets passing the eta cut
    jetsDelayedL1 = [j for j in jetsDelayedL1[:] if abs(j.Eta) < 3.2]   
    # Sort jets by highest pT
    jetsDelayedL1 = sorted(jetsDelayedL1, 
                         key = lambda j: j.PT, reverse=True)

    if (not jetsDelayedL1) or  (jetsDelayedL1[0].PT < 40.0):
        return evt_eff,evt_cutFlow
    
    evt_cutFlow['L1: 40 GeV < PT Jet1(N)'] += 1.0

    dphi_min = 2*np.pi
    for j in jetsDelayedL1[:6]:
        dphi = np.abs(j.Phi-metOnTime.Phi)
        if dphi > np.pi:
            dphi = 2*np.pi - dphi
        dphi_min = min(dphi,dphi_min)

    if dphi_min > 1.0:
        return evt_eff,evt_cutFlow
    
    evt_cutFlow['L1: DPhi(Jet(N),MET(N-1)) < 1.0'] += 1.0


    ### HLT Trigger
    # Get jets from the L1 (delayed) calorimeter
    jetsDelayedHLT = list(tree.HLTJetDelayed)
    jets = [j for j in jetsDelayedHLT[:] 
            if (abs(j.Eta) < 2.5 and j.PT > 20.0)]
    if not jets:
        return evt_eff,evt_cutFlow
    
    evt_cutFlow['HLT: Eta Jet < 2.5 PT Jet > 20 GeV'] += 1.0

    # For the remaining jets, find at least one with low
    # energy deposit in the ECAL cell closest to the jet
    jet_cells = []
    for j in jets:
        closest_cell = None
        dRmin = 100.0
        for tower_cell in tree.HLTTowerDelayed:
            dphi = np.abs(j.Phi-tower_cell.Phi)
            deta = j.Eta-tower_cell.Eta
            if dphi > np.pi:
                dphi = 2*np.pi - dphi
            dR = np.sqrt(deta**2 + dphi**2)
            if dR < dRmin:
                dRmin = dR
                closest_cell = tower_cell
        
        jet_cells.append(closest_cell)
    
    # Remove jets with large ECAL deposits
    jets_disp = [j for j,cell in zip(jets,jet_cells) 
                 if cell.Eem/(cell.Eem + cell.Ehad) < 0.06]
    if not jets_disp:
        return evt_eff,evt_cutFlow
    
    evt_cutFlow['HLT: Jet(N) EMF < 0.06'] += 1.0

    # Finally for the remaining jets remove jets with tracks close
    # to it and belonging to the delayed event
    tracksDelayed = []
    for track in tree.AllTracks:
        if not (25e-9 < track.T < 35e-9):
            continue
        if track.PT < 2.0:
            continue
        tracksDelayed.append(track)
    
    jets_clean = []
    for j in jets_disp:
        dRmin = 100.0
        for track in tracksDelayed:
            dphi = np.abs((j.Phi-track.Phi))
            if dphi > np.pi:
                dphi = 2*np.pi - dphi
            deta = (j.Eta-track.Eta)
            dR = np.sqrt(deta**2 + dphi**2)
            dRmin = min(dR,dRmin)
        if dRmin > 0.2:
            jets_clean.append(j)

    if not jets_clean:
        return evt_eff,evt_cutFlow
    
    evt_cutFlow['HLT: DR(Tracks(N),Jet(N)) > 0.2'] += 1.0


    evt_eff = 1.0
    
    return evt_eff,evt_cutFlow

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
    
    cutFlowKeys = defineCutFlow()
    evts_cutFlow = {key : 0.0 for key in cutFlowKeys}
    
    
    nevts = tree.GetEntries()
    for ievt in tqdm.tqdm(range(nevts),position=ijob,
                          desc=os.path.basename(inputFile),
                          leave=True):
        tree.GetEntry(ievt) 
        effsDict['Nevents'] += 1.0
        evt_eff,evt_cutFlow = getEffFor(tree)
        # Add event efficiency to total efficiencies 
        # #for the given selection (for each tau value)
        effsDict['Trigger'] += evt_eff
        for key in evts_cutFlow:
            evts_cutFlow[key] += evt_cutFlow[key]
        
    f.Close()
    
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
        effFile = outFile + output_suffix +'_b2tf_effs.csv'
        cutFile = outFile + output_suffix +'_b2tf_cutFlow.csv'
        cutFlowStr += f'CutFlow for {effsDict['inputFile']}:\n'
        for key,val in evts_cutFlow.items():
            cutFlowStr += f'  {key} : {val/evts_cutFlow['All']:1.3e}\n'
        saveOutput(effsDict,effFile,evts_cutFlow,cutFile)        

    logger.debug(cutFlowStr+'\n')
            
        
    logger.info("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
