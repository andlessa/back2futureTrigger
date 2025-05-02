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
import progressbar as P
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


random.seed(123)
 

def getDataFrom(llpParticles,llpDaughters,llp_PDGs,invisible_PDGs):

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
        for p in visibleLists[illp]:
            p_visible = p_visible + p.momentum
        for p in invisibleLists[illp]:
            p_invisible = p_invisible + p.momentum
    
        visParticle = pyhepmc.GenParticle(momentum = p_visible, pid = illp)
        invisParticle = pyhepmc.GenParticle(momentum = p_invisible, pid = illp)
        llp_momentum = pyhepmc.FourVector(llp.Px,llp.Py,llp.Pz,llp.E)
        llpParticle = pyhepmc.GenParticle(momentum = llp_momentum, pid = llp.PID)
        eventDict[illp] = { 'visible' : visParticle, 
                            'invisible' : invisParticle,
                            'parent' : llpParticle,
                            'visiblePDGs' : [p.pid for p in visible_particles],
                            'invisiblePDGs' : [p.pid for p in invisible_particles],
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



def decayTime(genParticle,tau0=1.0):
    gamma = genParticle.momentum.e/genParticle.momentum.m()
    tau = tau0*gamma / c
    t = random.random()
    tdecay = (-tau*np.log(t))*1e9 # Time in nano seconds

    return tdecay

def getDecayLength(genParticle,tdecay):
    
    beta = (genParticle.momentum/genParticle.momentum.e)
    L_vec = beta*c*(1e-9*tdecay) # convert decay time to seconds

    return L_vec.pt(),abs(L_vec.pz)

def passTimingCut(t_onTime,t_delayed):
    """
    Require that the first and second decays to take place
    in the intervals 0 < t1 < 10ns and 25ns < t2 < 35ns
    """

    if not (0.0 <= t_onTime <= 10.0):
        return False
    if not (25.0 <= t_delayed <= 35.0):
        return False
    
    return True


def passDecayLengthCut(r):
    """
    Require that the decay takes place within R < 4m.
    """

    if (r > 4.0):
        return False
    
    return True

def getJetET(j):

    ET = j.momentum.pt()*(j.momentum.e/j.momentum.p3mod())

    return ET

def passEnergyCut(j_ontime,j_delayed):
    """
    Apply the transverse energy cuts for the on-time visible
    decay (j1) and delayed visible decay (j2)
    """


    # Get HT for onTime event:
    ETonTime = getJetET(j_ontime)
    # ETonTime = np.sqrt(j1.momentum.pt()**2 + j1.momentum.m()**2)
    if not (40.0 < ETonTime < 100.0):
        return False
        # Get HT for delayed event:
    ETdelayed= getJetET(j_delayed)
    # ETdelayed = np.sqrt(j2.momentum.pt()**2 + j2.momentum.m()**2)
    if not (40.0 < ETdelayed):
        return False

    return True

def passEtaCut(j):
    """
    Apply the eta cut between the on-time visible
    decay (j1) and delayed visible decay (j2)
    """

    if abs(j.momentum.eta()) > 3.2:
        return False
    
    return True

def passDelaPhiCut(j1,j2):
    """
    Apply the deltaPhi cut between the on-time visible
    decay (j1) and delayed visible decay (j2)
    """

    dphi = pyhepmc.delta_phi(j1.momentum,j2.momentum)
    if abs(dphi) < np.pi -1.0:
        return False
    
    return True


def getEventRecordsFrom(tree,ievt,llp_PDGs,invisible_PDGs,tau0=1.0):

    tree.GetEntry(ievt)
    # Check how many LLPs were found in event
    nLLPs = len(tree.llpParticles)
    nLLPs_A = len(tree.llpParticlesA)
    nLLPs_B = len(tree.llpParticlesB)
    if nLLPs != nLLPs_A+nLLPs_B:
        errorMsg = f"The total number of LLPs ({nLLPs}) does not correspond to the sum of type A ({nLLPs_A}) and type B ({nLLPs_B}) LLPs."
        logger.error(errorMsg)
        return None
    elif nLLPs > 2:
        errorMsg = f"{nLLPs} LLPs found in event (can only handle with up to 2 LLPs)"
        logger.error(errorMsg)
        return None
    
    # Get event records based on how many LLPs are present in the event
    # (on-time = N-1 event record, delayed = N event record)
    eventRecords = {'on-time' : {}, 'delayed' : {}}

    # If nLLPs = 0, set the current event as N-1 and the next event as the N event record.
    if nLLPs == 0: 
        onTimeER = {'jets' : getJetsFrom(tree.GenJet), 'LLPdata' : {}}
        # Get the next event record
        tree.GetEntry(ievt+1)
        delayedER = {'jets' : getJetsFrom(tree.GenJet), 'LLPdata' : {}}
        # Go back to the current event
        tree.GetEntry(ievt) 
        # If no LLPs are present, set decay times to each event record time
        onTimeER['LLPdata'].update({'tdecay_0': 0.0, 'tau0' : 0.0, 'Lxy_0' : 0.0})
        delayedER['LLPdata'].update({'tdecay_0': 25.0, 'tau0' : 0.0, 'Lxy_0' : 0.0})

    # If nLLPs = 1, then we require its decay to correspond to the N event record and everything else is the N-1 event record
    elif nLLPs == 1:
        if nLLPs_B == 0: # There are no particles of type B, so A must be the delayed particle
            llpDict = getDataFrom(tree.llpParticlesA,
                                   tree.llpDirectDaughtersA,
                                    llp_PDGs=llp_PDGs,invisible_PDGs=invisible_PDGs)
            onTimeER = {'jets' : getJetsFrom(tree.GenJetBOnTime), 'LLPdata' : {}}
            delayedER = {'jets' : getJetsFrom(tree.GenJetBDelayed), 'LLPdata' : llpDict}
        elif nLLPs_A == 0: # There are no particles of type A, so B must be the delayed particle:
            llpDict = getDataFrom(tree.llpParticlesB,
                                   tree.llpDirectDaughtersB,
                                    llp_PDGs=llp_PDGs,invisible_PDGs=invisible_PDGs)
            onTimeER = {'jets' : getJetsFrom(tree.GenJetAOnTime), 'LLPdata' : {}}
            delayedER = {'jets' : getJetsFrom(tree.GenJetADelayed), 'LLPdata' : llpDict}
        tdecay =  decayTime(llpDict[0]['parent'],tau0)
        Lxy_0,_ = getDecayLength(llpDict[0]['parent'],tdecay)
        # Set prompt "decays" for the case of no LLPs in the on-time event
        onTimeER['LLPdata'].update({'tdecay_0': 0.0, 'tau0' : 0.0, 'Lxy_0' : 0.0})
        delayedER['LLPdata'].update({'tdecay_0': tdecay, 'tau0' : tau0, 'Lxy_0' : Lxy_0})
    
    # If nLLPs = 2, then we decide on which decay leads to the N-1 event record and which to the N event record based on their decay times
    elif nLLPs == 2:
            llpDictA = getDataFrom(tree.llpParticlesA,
                                   tree.llpDirectDaughtersA,
                                    llp_PDGs=llp_PDGs,invisible_PDGs=invisible_PDGs)
            llpDictB = getDataFrom(tree.llpParticlesB,
                                   tree.llpDirectDaughtersB,
                                    llp_PDGs=llp_PDGs,invisible_PDGs=invisible_PDGs)
            # (note that only the ordering of the decays matter and not their absolute value, so the value of tau0 is irrelevant)
            tdecayA =  decayTime(llpDictA[0]['parent'],tau0)
            Lxy_A,_ = getDecayLength(llpDictA[0]['parent'],tdecayA)
            tdecayB =  decayTime(llpDictB[0]['parent'],tau0)
            Lxy_B,_ = getDecayLength(llpDictB[0]['parent'],tdecayB)
            if tdecayA < tdecayB: # LLPA is on-time and LLPB delayed
                onTimeER = {'jets' : getJetsFrom(tree.GenJetAOnTime), 'LLPdata' : llpDictA}
                delayedER = {'jets' : getJetsFrom(tree.GenJetADelayed), 'LLPdata' : llpDictB}
                onTimeER['LLPdata'].update({'tdecay_0': tdecayA, 'tau0' : tau0, 'Lxy_0' : Lxy_A})
                delayedER['LLPdata'].update({'tdecay_0': tdecayB, 'tau0' : tau0, 'Lxy_0' : Lxy_B})

            else:
                onTimeER = {'jets' : getJetsFrom(tree.GenJetBOnTime), 'LLPdata' : llpDictB}
                delayedER = {'jets' : getJetsFrom(tree.GenJetBDelayed), 'LLPdata' : llpDictA}
                onTimeER['LLPdata'].update({'tdecay_0': tdecayB, 'tau0' : tau0, 'Lxy_0' : Lxy_B})
                delayedER['LLPdata'].update({'tdecay_0': tdecayA, 'tau0' : tau0, 'Lxy_0' : Lxy_A})
    

    eventRecords['on-time'] = onTimeER
    eventRecords['delayed'] = delayedER

    return eventRecords

def getEffFor(eventRecords,tauList):

    evt_effs = np.zeros(len(tauList))
    evt_cutFlow = {'All' : np.zeros(len(tauList)),
                    'Eta cut' : np.zeros(len(tauList)),
                    'ETmax(N-1) < 100 GeV' : np.zeros(len(tauList)),
                    'Energy theshold' : np.zeros(len(tauList)),
                    'Delta Phi' : np.zeros(len(tauList)),              
                    'Decay Position' : np.zeros(len(tauList)),
                    'Decay Time' : np.zeros(len(tauList)),
                   }
    
    onTimeER = eventRecords['on-time']
    delayedER = eventRecords['delayed']

    for i,tau in enumerate(tauList):

        # Rescale decay times and lengths by tau value
        ratio = 1.0
        if onTimeER['LLPdata']['tau0'] != 0.0:
            ratio = tau/onTimeER['LLPdata']['tau0']
        # The decay time and position scale with tau/tau0
        t_decay_ontime = onTimeER['LLPdata']['tdecay_0']*ratio
        Lxy_ontime = onTimeER['LLPdata']['Lxy_0']*ratio
        
        ratio = 1.0
        if delayedER['LLPdata']['tau0'] != 0.0:
            ratio = tau/delayedER['LLPdata']['tau0']
        # The decay time and position scale with tau/tau0
        t_decay_delayed = delayedER['LLPdata']['tdecay_0']*ratio
        Lxy_delayed = delayedER['LLPdata']['Lxy_0']*ratio

        jets_ontime = onTimeER['jets']
        jets_delayed = delayedER['jets']
   
        evt_cutFlow['All'][i] += 1

        # Keep only jets passing the eta cut
        jets_ontime = [j for j in jets_ontime[:] if passEtaCut(j)]
        jets_delayed = [j for j in jets_delayed[:] if passEtaCut(j)] 

        if len(jets_ontime) == 0:
            continue
        if len(jets_delayed) == 0:
            continue
        # # Use leading jet from each event record
        # j_ontime = jets_ontime[0]
        # j_delayed = jets_delayed[0]

        # if not passEtaCut(j_ontime):
        #     continue
        # if not passEtaCut(j_delayed):
        #     continue
        evt_cutFlow['Eta cut'][i] += 1

        # Apply cut on hardest allowed jet for the on-time event
        if getJetET(jets_ontime[0]) > 100.0:
            continue

        evt_cutFlow['ETmax(N-1) < 100 GeV'][i] += 1

        # Keep only the hardest
        j_ontime = jets_ontime[0]
        j_delayed = jets_delayed[0]

        if not passEnergyCut(j_ontime,j_delayed):
            continue
        evt_cutFlow['Energy theshold'][i] += 1

        if not passDelaPhiCut(j_ontime,j_delayed):
            continue
        evt_cutFlow['Delta Phi'][i] += 1

        if not passDecayLengthCut(Lxy_ontime):
            continue
        if not passDecayLengthCut(Lxy_delayed):
            continue

        evt_cutFlow['Decay Position'][i] += 1

        if not passTimingCut(t_decay_ontime,t_decay_delayed):
            continue
        evt_cutFlow['Decay Time'][i] += 1

        evt_effs[i] += 1.0
        
    return evt_effs,evt_cutFlow

def getEfficiencies(inputFile,tauList,
                    llps=[4000023],invisibles=[4000022]):

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
    effsDict = {"Trigger" : np.zeros(len(tauList)),
                "Nevents" : 0}
    
    evts_cutFlow = {'All' : np.zeros(len(tauList)),
                    'Eta cut' : np.zeros(len(tauList)),
                    'ETmax(N-1) < 100 GeV' : np.zeros(len(tauList)),
                    'Energy theshold' : np.zeros(len(tauList)),
                    'Delta Phi' : np.zeros(len(tauList)),              
                    'Decay Position' : np.zeros(len(tauList)),
                    'Decay Time' : np.zeros(len(tauList)),
                   }
    
    
    nevts = tree.GetEntries()
    for ievt in range(nevts):
        evtRecords = getEventRecordsFrom(tree,ievt,
                                         llp_PDGs=llps,
                                         invisible_PDGs=invisibles,
                                         tau0=tauList[0])
        effsDict['Nevents'] += 1
        evt_effs,evt_cutFlow = getEffFor(evtRecords,tauList)
        # Add event efficiency to total efficiencies 
        # #for the given selection (for each tau value)
        effsDict['Trigger'] += np.array(evt_effs)
        for key in evts_cutFlow:
            evts_cutFlow[key] += evt_cutFlow[key]
        
    f.Close()

    # Divide the total efficiency by the number of events:
    effsDict["TriggerErr"] = np.sqrt(effsDict["Trigger"])/effsDict['Nevents']
    effsDict["Trigger"] = effsDict["Trigger"]/effsDict['Nevents']

    # Store ctau list
    effsDict['ctau'] = tauList[:]
    evts_cutFlow['ctau'] = tauList[:]
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
    columns = []
    columnsLabels = []
    for key,val in effsDict.items():
        columnsLabels.append(key)
        if isinstance(val,(float,int)):
            columns.append(itertools.repeat(val))
        else:
            columns.append(val)
    data = np.array(list(zip(*columns)))
    header = ','.join(columnsLabels)
    
    np.savetxt(outputFile, data, 
                header=f'Input file: {inputFile}\nNumber of events: {nevts}\n{header}',
                delimiter=',',fmt='%1.3e')

    # Save cutflow (if defined)
    # Get column labels and data
    columns = []
    columnsLabels = []
    for key,val in evts_cutFlow.items():
        columnsLabels.append(key)
        if isinstance(val,(float,int)):
            columns.append(itertools.repeat(val))
        else:
            columns.append(val)
    data = np.array(list(zip(*columns)))
    header = ','.join(columnsLabels)
    np.savetxt(cutFile, data, 
                header=f'Input file: {inputFile}\nNumber of events: {nevts}\n{header}',
                delimiter=',',fmt='%1.3e')

if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Compute the efficiencies for a given input HepMC file" )
    ap.add_argument('-f', '--inputfile',nargs='+',
            help='path to the HepMC file(s).')
    ap.add_argument('-p', '--parfile',default='parameters_getEff.ini',
            help='parameters file name, where additional options are defined [parameters_getEff.ini].')
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

    parser = configparser.ConfigParser(inline_comment_prefixes="#")   
    ret = parser.read(args.parfile)
    if ret == []:
        logger.error(f"No such file or directory: {args.parfile}")
        sys.exit()

    
    try:
        tauList = parser.get("options","ctau_values")
        tauList = eval(tauList)
    except:
        logger.error("Error getting list of ctau values")
        tauList = np.geomspace(args.tmin,args.tmax,args.ntau)

    try:
        llpPDGs = parser.get("model","llp_pdgs").split(',')
        llpPDGs = [int(pdg) for pdg in llpPDGs]
    except:
        logger.error("Error getting list of LLP PDG values")
        llpPDGs = [35]

    try:
        invisiblePDGs = parser.get("model","invisible_pdgs").split(',')
        invisiblePDGs = [int(pdg) for pdg in invisiblePDGs]
    except:
        logger.error("Error getting list of LLP PDG values")
        invisiblePDGs = [12,14,16]

    t0 = time.time()


    # Check for output file suffix
    output_suffix = ""
    if parser.has_option("options","output_suffix"):
        output_suffix = parser.get("options","output_suffix")
   

    # Split input files by distinct models and get recast data for
    # the set of files from the same model:
    ncpus = int(parser.get("options","ncpus"))
    inputFileList = args.inputfile
    ncpus = min(ncpus,len(inputFileList))
    if ncpus == 1:
        effsDict,evts_cutFlow = getEfficiencies(inputFileList[0],
                                                tauList,
                                                llpPDGs,
                                                invisiblePDGs)
        outFile = effsDict['inputFile'].split('.root')[0]
        effFile = outFile + output_suffix  +'_b2tf_effs_v3.csv'
        cutFile = outFile + output_suffix +'_b2tf_cutFlow_v3.csv'
        saveOutput(effsDict,effFile,evts_cutFlow,cutFile)
    else:
        pool = multiprocessing.Pool(processes=ncpus)
        children = []
        for inputfile in inputFileList:
            p = pool.apply_async(getEfficiencies, args=(inputfile,tauList,
                            llpPDGs,invisiblePDGs))
            children.append(p)

        nfiles = len(inputFileList)
        progressbar = P.ProgressBar(widgets=[f"Reading {nfiles} files", 
                                    P.Percentage(),P.Bar(marker=P.RotatingMarker()), P.ETA()])
        progressbar.maxval = nfiles
        progressbar.start()
        ndone = 0
        for p in children: 
            effsDict,evts_cutFlow = p.get()
            outFile = effsDict['inputFile'].split('.root')[0].split('.hepmc')[0]
            effFile = outFile + output_suffix +'_b2tf_effs_v3.csv'
            cutFile = outFile + output_suffix +'_b2tf_cutFlow_v3.csv'
            saveOutput(effsDict,effFile,evts_cutFlow,cutFile)
            ndone += 1
            progressbar.update(ndone)
        
    logger.info("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
