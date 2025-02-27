#!/usr/bin/env python3

import os, time, sys
import numpy as np
import readMapNew as rmN
import pyhepmc
import random
import logging
import configparser
import multiprocessing
import progressbar as P

c = 3e8

FORMAT = '%(levelname)s: %(message)s at %(asctime)s'
logging.basicConfig(format=FORMAT,datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger()   


random.seed(123)

def getLLPdecays(event,llps=[35]):

    vertexDict = {}
    for ivertex,vertex in enumerate(event.vertices):
        if len(vertex.particles_in) != 1:
            continue
        if len(vertex.particles_out) == 1:
            continue
        part_in = vertex.particles_in[0]
        # Check if incoming particle appear in pdgList
        if abs(part_in.pid) not in llps:
            continue

        part_out_list = vertex.particles_out
        # Skip vertices which could correspond to FSR (initial state appears in final state)
        if any(p.pid == part_in.pid for p in part_out_list):
            continue

        vertexDict[ivertex] = vertex
    
    return vertexDict

def getDataFrom(event,llps=[35],invisibles=[12,14,16]):

    vertexDict = getLLPdecays(event,llps=llps)
    eventDict = {}
    if len(vertexDict) != 2:
        logger.error(f"{len(vertexDict)} LLP decay vertices found (can only deal with 2 decay vertices)")
        raise ValueError(f"{len(vertexDict)} LLP decay vertices found (can only deal with 2 decay vertices)")
    for ivertex,vertex in vertexDict.items():
        visible_particles = [p for p in vertex.particles_out 
                                if abs(p.pid) not in invisibles]
        if len(visible_particles) != 2:
            logger.error(f"{len(visible_particles)} visible particles found in LLP decay (can only deal with 2 particles)")
            raise ValueError(f"{len(visible_particles)} visible particles found in LLP decay (can only deal with 2 particles)")
        p_visible = pyhepmc.FourVector(0.,0.,0.,0.)
        for p in visible_particles:
            p_visible = p_visible + p.momentum

        visParticle = pyhepmc.GenParticle(momentum = p_visible, pid = ivertex)
        eventDict[ivertex] = {'visible' : visParticle, 
                              'parent' : vertex.particles_in[0]}

    return eventDict

def lifetime(avgtau = 4.3):
    avgtau = avgtau / c
    t = random.random()
    return -1.0 * avgtau * np.log(t)

def getDecayLength(genParticle,tau):
    
    lt = lifetime(tau) # set mean lifetime
    gamma = genParticle.momentum.e/genParticle.momentum.m()
    L_vec = (genParticle.momentum/genParticle.momentum.e)*c*lt*gamma

    return L_vec.pt(),abs(L_vec.pz)

def getEfficiencies(hepmcFile,tauList,
                    llps=[35],invisibles=[12,14,16]):

    # Load hepmc File
    if not os.path.isfile(hepmcFile):
        logger.error(f"File {hepmcFile} not found")
        return None
    
    try:
        f = pyhepmc.open(hepmcFile) # Open HEPMC file
    except:
        logger.error(f"Error reading {hepmcFile}")
        return None
    
    # Total efficiencies (for each tau value)
    effsDict = {"high-ET" : np.zeros(len(tauList)),
                "low-ET" : np.zeros(len(tauList)),
                "Nevents" : 0}
    
    # Get next event
    event = f.read()
    while event is not None:
        effsDict["Nevents"] += 1
        # Extract necessary data from event
        eventDict = getDataFrom(event,llps=llps,invisibles=invisibles)
        if len(eventDict) != 2:
            logger.error(f"{len(eventDict)} LLP decay vertices found (can only deal with 2 decay vertices)")
            raise ValueError(f"{len(eventDict)} LLP decay vertices found (can only deal with 2 decay vertices)")
        llpList = [d['parent'] for d in eventDict.values()]
        visList = [d['visible'] for d in eventDict.values()]
        # llp_mass = llpList[0].generated_mass
        # if any(llp.generated_mass != llp_mass for llp in llpList):
        #     raise ValueError(f"Can not deal with distinct masses found for LLPs.")
        
        # Get total momentum of LLP pair (equals momentum of parent)
        pTot = llpList[0].momentum + llpList[1].momentum
        if pTot.m() >= 400:
            sr = "high-ET"
        else:
            sr = "low-ET"

        
        p1 = visList[0].momentum        
        p1_pt = p1.pt()
        p1_eta = p1.eta()
        p1_pdgs = set([abs(c.pid) for c in llpList[0].children 
                   if abs(c.pid) not in invisibles])
        if len(p1_pdgs) != 1:
            raise ValueError(f'Can not handle decays to distinct pair of fermions (e.g. {p1_pdgs})')
        else:
            p1_pdgs = list(p1_pdgs)[0]
        
        p2 = visList[1].momentum
        p2_pt = p2.pt()
        p2_eta = p2.eta()
        p2_pdgs = set([abs(c.pid) for c in llpList[1].children 
                   if abs(c.pid) not in invisibles])
        if len(p2_pdgs) != 1:
            raise ValueError(f'Can not handle decays to distinct pair of fermions (e.g. {p2_pdgs})')
        else:
            p2_pdgs = list(p2_pdgs)[0]
        
        
        evt_effs = []
        for tau in tauList:
            L1xy,L1z = getDecayLength(llpList[0],tau)
            L2xy,L2z = getDecayLength(llpList[1],tau)

            eff = rmN.queryMapFromKinematics(p1_pt,p1_eta,L1xy,L1z,p1_pdgs,
                                             p2_pt,p2_eta,L2xy,L2z,p2_pdgs,
                                             selection = sr)
            evt_effs.append(eff)

        # Add event efficiency to total efficiencies 
        # #for the given selection (for each tau value)
        effsDict[sr] += np.array(evt_effs)
        # Get next event
        event = f.read()

    f.close()
    # Divide the total efficiency by the number of events:
    effsDict["high-ET"] = effsDict["high-ET"]/effsDict['Nevents']
    effsDict["low-ET"] = effsDict["low-ET"]/effsDict['Nevents']
    # Store ctau list
    effsDict['ctau'] = tauList[:]
    # Store input file name
    effsDict['hepmcFile'] = hepmcFile

    return effsDict

def saveOutput(effsDict):
        
    tauList = effsDict['ctau']
    data = np.array(list(zip(tauList,effsDict['low-ET'],
                                effsDict['high-ET'])))

    outputFile = effsDict['hepmcFile'].split('.hepmc')[0]+'_effs.csv'
    np.savetxt(outputFile, data, 
                header=f'Input file: {effsDict['hepmcFile']}\nNumber of events: {effsDict['Nevents']}\nctau(m),eff(low-ET),eff(high-ET)',
                delimiter=',',fmt='%1.3e')
    


if __name__ == "__main__":
    
    import argparse    
    ap = argparse.ArgumentParser( description=
            "Compute the efficiencies for a given input HepMC file" )
    ap.add_argument('-f', '--inputfile',nargs='+',
            help='path to the HepMC file(s).')
    ap.add_argument('-p', '--parfile',default='parameters.ini',
            help='parameters file name, where additional options are defined [parameters.ini].')
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


    # Split input files by distinct models and get recast data for
    # the set of files from the same model:
    ncpus = parser.get("options","ncpus")
    inputFileList = args.inputFile
    ncpus = min(ncpus,len(inputFileList))
    pool = multiprocessing.Pool(processes=ncpus)
    children = []
    for inputfile in inputFileList:
        p = pool.apply_async(getEfficiencies, args=(inputfile,tauList,
                           llpPDGs,invisiblePDGs,),callback=saveOutput)
        children.append(p)

    nfiles = len(inputFileList)
    progressbar = P.ProgressBar(widgets=[f"Reading {nfiles} HepMC files", 
                                P.Percentage(),P.Bar(marker=P.RotatingMarker()), P.ETA()])
    progressbar.maxval = nfiles
    progressbar.start()
    ndone = 0
    for p in children: 
        p.get()
        ndone += 1
        progressbar.update(ndone)
        
    logger.info("\n\nDone in %3.2f min" %((time.time()-t0)/60.))
