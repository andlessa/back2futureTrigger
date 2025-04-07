#!/usr/bin/env python3

import os, time, sys
import numpy as np
import pyhepmc
import random
import logging
import configparser
import multiprocessing
import progressbar as P
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

def getDataFrom(tree,llps,invisibles):

    # Consistency checks
    if any(abs(p.PID) not in llps for p in tree.llpParticles):
        erroMsg = "LLP PDG from input file does not match definitions"
        logger.error(erroMsg)
        raise ValueError(erroMsg)
    
    if any(abs(p.PID) not in invisibles for p in tree.dmParticles):
        erroMsg = "DM PDG from input file does not match definitions"
        logger.error(erroMsg)
        raise ValueError(erroMsg)
    

    llpParticles = tree.llpParticles
    llpDaughters = tree.llpDirectDaughters
    eventDict = {}
    if len(llpParticles) != 2:
        erroMsg = f"{len(llps)} LLP decay vertices found (can only handle 2 decay vertices)"
        logger.error(erroMsg)
        raise ValueError(erroMsg)
    for illp,llp in enumerate(llpParticles):
        # Select set of daughters from the current llp particle
        llp_daughters = [d for d in llpDaughters if d.M1 == illp]
        visible_particles = []
        invisible_particles = []
        for d in llp_daughters:
            mom = pyhepmc.FourVector(d.Px,d.Py,d.Pz,d.E)            
            p = pyhepmc.GenParticle(momentum = mom, pid = d.PID)
            p.generated_mass = d.Mass
            if abs(p.pid) in invisibles:
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
        
        p_visible = pyhepmc.FourVector(0.,0.,0.,0.)
        p_invisible = pyhepmc.FourVector(0.,0.,0.,0.)
        for p in visible_particles:
            p_visible = p_visible + p.momentum
        for p in invisible_particles:
            p_invisible = p_invisible + p.momentum

        visParticle = pyhepmc.GenParticle(momentum = p_visible, pid = illp)
        invisParticle = pyhepmc.GenParticle(momentum = p_invisible, pid = illp)
        llp_momentum = pyhepmc.FourVector(llp.Px,llp.Py,llp.Pz,llp.E)
        llpParticle = pyhepmc.GenParticle(momentum = llp_momentum, pid = llp.PID)
        eventDict[illp] = {'visible' : visParticle, 
                              'invisible' : invisParticle,
                              'parent' : llpParticle,
                              'visiblePDGs' : [p.pid for p in visible_particles],
                              'invisiblePDGs' : [p.pid for p in invisible_particles],
                              }

    return eventDict

def decayTime(genParticle,tau0):
    gamma = genParticle.momentum.e/genParticle.momentum.m()
    tau = tau0*gamma / c
    t = random.random()
    tdecay = (-tau*np.log(t))*1e9 # Time in nano seconds

    return tdecay

def getDecayLength(genParticle,tdecay):
    
    beta = (genParticle.momentum/genParticle.momentum.e)
    L_vec = beta*c*(1e-9*tdecay) # convert decay time to seconds

    return L_vec.pt(),abs(L_vec.pz)

def passTimingCut(t1ns,t2ns):
    """
    Require that the first and second decays to take place
    in the intervals 0 < t1 < 10ns and 25ns < t2 < 35ns
    """

    t_onTime = min(t1ns,t2ns)
    t_delayed = max(t1ns,t2ns)
    if not (0.0 < t_onTime < 10.0):
        return False
    if not (25.0 < t_delayed > 35.0):
        return False
    
    return True

def passDecayLengthCut(r1,r2):
    """
    Require that both decays take place within R < 4m.
    """

    if r1 > 4.0 or r2 > 4.0:
        return False
    
    return True

def passEnergyCut(j1,j2):
    """
    Apply the transverse energy cuts for the on-time visible
    decay (j1) and delayed visible decay (j2)
    """


    # Get HT for onTime event:
    ETonTime = j1.momentum.pt()
    if not (40.0 < ETonTime < 100.0):
        return False
        # Get HT for delayed event:
    ETdelayed= j2.momentum.pt()
    if not (40.0 < ETdelayed):
        return False

    return True

def passAngleCuts(j1,j2):
    """
    Apply the deltaPhi and deltaEta cuts between the on-time visible
    decay (j1) and delayed visible decay (j2)
    """

    dphi = pyhepmc.delta_phi(j1.momentum,j2.momentum)
    if abs(dphi) < np.pi -1.0:
        return False
    deta = pyhepmc.delta_eta(j1.momentum,j2.momentum)
    if abs(deta) > 3.2:
        return False
    
    return True



def getEffFor(tree,tauList,llps,invisibles):


    evt_effs = np.zeros(len(tauList))
    evt_cutFlow = {'Total' : np.zeros(len(tauList)),
                    'TimingCut' : np.zeros(len(tauList)),
                    'DecayLengthCut' : np.zeros(len(tauList)),
                   'EnergyCut' : np.zeros(len(tauList)),
                   'AngleCuts' : np.zeros(len(tauList))}

    # Extract necessary data from event
    eventDict = getDataFrom(tree,llps=llps,invisibles=invisibles)
    if len(eventDict) != 2:
        errorMsg = f"{len(eventDict)} LLP decay vertices found (can only handle 2 decay vertices)"
        logger.error(errorMsg)
        # raise ValueError(errorMsg)
        return evt_effs # return zero effs
    
    # Loop over tau values and compute the decay
    # and time positions:
    for i,tau in enumerate(tauList):

        evt_cutFlow['Total'][i] += 1

        t1_decay = decayTime(eventDict[0]['parent'],tau)
        t2_decay = decayTime(eventDict[1]['parent'],tau)
        if not passTimingCut(t1_decay,t2_decay):
            continue
        evt_cutFlow['TimingCut'][i] += 1
        L1xy,_ = getDecayLength(eventDict[0]['parent'],t1_decay)
        L2xy,_ = getDecayLength(eventDict[1]['parent'],t2_decay)
    
        if not passDecayLengthCut(L1xy,L2xy):
            continue
        evt_cutFlow['DecayLengthCut'][i] += 1
        # Get visible particles for the on-time decay:
        if t1_decay < t2_decay:
            eventOnTimeDict = eventDict[0]
            eventDelayedDict = eventDict[1]
        else:
            eventOnTimeDict = eventDict[1]
            eventDelayedDict = eventDict[0]

        j1 = eventOnTimeDict['visible']
        j2 = eventDelayedDict['visible']


        if not passEnergyCut(j1,j2):
            continue
        evt_cutFlow['EnergyCut'][i] += 1
        if not passAngleCuts(j1,j2):
            continue
        evt_cutFlow['AngleCuts'][i] += 1
        evt_effs[i] += 1.0
        
    return evt_effs,evt_cutFlow


def getEfficiencies(inputFile,tauList,
                    llps=[35],invisibles=[12,14,16]):

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
    
    evts_cutFlow = {'Total' : np.zeros(len(tauList)),
                    'TimingCut' : np.zeros(len(tauList)),
                    'DecayLengthCut' : np.zeros(len(tauList)),
                   'EnergyCut' : np.zeros(len(tauList)),
                   'AngleCuts' : np.zeros(len(tauList))}
    
    
    nevts = tree.GetEntries()
    for ievt in range(nevts):
        tree.GetEntry(ievt)   
        effsDict['Nevents'] += 1
        evt_effs,evt_cutFlow = getEffFor(tree,tauList,llps,invisibles)
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
    # Store input file name
    effsDict['inputFile'] = inputFile

    ntotal = evts_cutFlow['Total']
    print(evts_cutFlow)
    for key in evts_cutFlow:
        print(key,f'{evts_cutFlow[key]/ntotal:1.2e}')

    return effsDict

def saveOutput(effsDict,outputFile):
        
    tauList = effsDict['ctau']
    data = np.array(list(zip(tauList,effsDict['Trigger'],effsDict['TriggerErr'])))

    
    np.savetxt(outputFile, data, 
                header=f'Input file: {effsDict['inputFile']}\nNumber of events: {effsDict['Nevents']}\nctau(m),eff(Trigger),effErr',
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
        effsDict = getEfficiencies(inputFileList[0],tauList,
                                  llpPDGs,invisiblePDGs)
        outFile = effsDict['inputFile'].split('.root')[0]
        outFile = outFile + output_suffix  +'_b2tf_effs.csv'
        saveOutput(effsDict,outFile)
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
            effsDict = p.get()
            outFile = effsDict['inputFile'].split('.root')[0].split('.hepmc')[0]
            outFile = outFile + output_suffix +'_b2tf_effs.csv'
            saveOutput(effsDict,outFile)
            ndone += 1
            progressbar.update(ndone)
        
    logger.info("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
