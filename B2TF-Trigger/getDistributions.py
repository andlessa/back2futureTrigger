#!/usr/bin/env python3
import re
import numpy as np
import pandas as pd
import os
import logging

c = 3e8

FORMAT = '%(levelname)s: %(message)s at %(asctime)s'
logging.basicConfig(format=FORMAT,datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger()


# ### Get kinematical variables for each file
L1metStr = r'$E_T^{\rm miss}$ (Hardware Trigger, N-1) (GeV)'
L1njStr = r'$n_{j}$ (Hardware Trigger, N)'
L1pTj1Str = r'Leading jet $E_{T}$ (Hardware Trigger, N) (GeV)'
L1dPhi = r'$\Delta \phi^{min} (MET,j)$ (Hardware Trigger)'
L1metPartonStr = r'$E_T^{\rm miss}$ (Parton Level, N-1) (GeV)'

HLTpTj1Str = r'Leading jet $E_{T}$ (Off-line, N) (GeV)'
HLTnjStr = r'$n_{j}$ (Off-line, N)'
HLTemfStr = r'EMF$_{\rm min}$ (Off-line, N)'

bb1dPhi = r'$\Delta \phi (b_1,\bar{b}_1)$ (Parton Level)'
bb1dR = r'$\Delta R (b_1,\bar{b}_1)$ (Parton Level)'
bb2dPhi = r'$\Delta \phi (b_2,\bar{b}_2)$ (Parton Level)'
bb2dR = r'$\Delta R (b_2,\bar{b}_2)$ (Parton Level)'


cols = [L1metStr,L1njStr,L1pTj1Str,L1dPhi,L1metPartonStr,HLTnjStr,HLTpTj1Str,HLTemfStr,
        bb1dPhi,bb1dR,bb2dPhi,bb2dR]


# ### Define CSV filename
csvFiles = {L1metStr : './results/l1_met.csv',
                L1njStr : './results/l1_nj.csv',
                L1pTj1Str : './results/l1_et.csv',
                L1dPhi : './results/l1_del_phi.csv',
                L1metPartonStr : './results/l1_metParton.csv',
                HLTnjStr : './results/offline_nj.csv',
                HLTpTj1Str : './results/offline_pt.csv',
                HLTemfStr : './results/offline_EMF.csv',
                bb1dPhi : './results/bbar1_dPhi.csv',
                bb1dR : './results/bbar1_dR.csv',
                bb2dPhi : './results/bbar2_dPhi.csv',
                bb2dR : './results/bbar2_dR.csv',
                }

if __name__ == "__main__":
    

    import argparse
    import progressbar as P
    from tqdm import tqdm
    ap = argparse.ArgumentParser( description=
            "Read the output of a Delphes ROOT file and save the distributions to CSV files." )
    ap.add_argument('-f', '--inputfile',
            help='path to the ROOT file.',required=True)
    ap.add_argument('-l', '--label',
            help='label to the output file.',required=True)    
    ap.add_argument('-v', '--verbose', default='info',
            help='verbose level (debug, info, warning or error). Default is info')
    

    args = ap.parse_args()
    file = args.inputfile
    if not os.path.isfile(file):
        raise ValueError(f'File {file} does not exist!')
    
    label = args.label
        

    delphesDir = os.path.abspath("../DelphesLLP")
    os.environ['ROOT_INCLUDE_PATH'] = os.path.join(delphesDir,"external")

    import ROOT
    ROOT.gSystem.Load(os.path.join(delphesDir,"libDelphes.so"))

    ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
    ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
    ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')

    dataList = []
    f = ROOT.TFile(file,'read')
    tree = f.Get("Delphes")
    nevts = tree.GetEntries()
    print(f"Reading {file}\n")
        
    for ievt in tqdm(range(nevts)):
        tree.GetEntry(ievt)

        # Get parton level MET and b-bar angular separation
        llps = list(tree.llpParticles)
        invisibles = [p for p in tree.llpDirectDaughters 
                    if abs(p.PID) == 4000022]
        bquarks = [(p.M1,p) for p in tree.llpDirectDaughters 
                    if abs(p.PID) == 5]
        # Group bquarks by parent
        bquark_pairs = {illp : [] for illp,_ in bquarks}
        for illp,b  in bquarks:
            bquark_pairs[illp].append(b)
        invisibles = sorted(invisibles, key = lambda p: p.M1)
        pInvTot = np.zeros(3)
        for illp,llp in enumerate(llps):
            daughter = invisibles[illp]
            decayTime = daughter.T
            if decayTime < 10e-9: # if LLP decays on-time, add its daughter momentum
                pInv = np.array([daughter.Px,daughter.Py,
                                daughter.Pz])
            else: # add the LLP momentum
                pInv = np.array([llp.Px,llp.Py,
                                llp.Pz])
            pInvTot += pInv
        
        metParton = np.linalg.norm(pInvTot[:2])

        bb1dphi = 0.0
        bb1dr = 0.0
        bb2dphi = 0.0
        bb2dr = 0.0
        for illp,bpair in bquark_pairs.items():
            if len(bpair) != 2:
                continue
            b = bpair[0]
            bbar = bpair[1]
            bbdphi = np.abs(b.Phi-bbar.Phi)
            if bbdphi > np.pi:
                bbdphi = 2*np.pi-bbdphi
            bbdr = np.sqrt((b.Eta-bbar.Eta)**2 + bbdphi**2)
            if illp == 0:
                bb1dphi = bbdphi
                bb1dr = bbdr
            else:
                bb2dphi = bbdphi
                bb2dr = bbdr


        # Get detector level variables
                
        metOnTime = tree.L1METOnTime.At(0)
        jetsDelayed = list(tree.L1JetDelayed)
        jetsDelayed = sorted(jetsDelayed, 
                        key = lambda j: j.PT, reverse=True)

        met = metOnTime.MET
        nj = len(jetsDelayed)
        if nj > 0:
            pTj1 = jetsDelayed[0].PT
            dphi_min = 10000.0
            for j in jetsDelayed[:6]:
                dphi = np.abs(j.Phi-metOnTime.Phi)
                if dphi > np.pi:
                    dphi = 2*np.pi-dphi
                dphi_min = min(dphi,dphi_min)
        else:
            pTj1 = 0.0
            dphi_min = 5.0


        jetsDelayedHLT = list(tree.HLTJetDelayed)
        jets = [j for j in jetsDelayedHLT[:] if abs(j.Eta) < 2.5]
        jets = [j for j in jets[:] if j.PT > 20.0]
        jets = sorted(jets, key = lambda j: j.PT, reverse=True)
        njHLT = len(jets)     
        pTHLT = 0.0
        emf_min = -1.0   
        if njHLT > 0:
            pTHLT = jets[0].PT
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
            # Compute minimum EMF:
            emf_min = min([cell.Eem/(cell.Eem + cell.Ehad) for cell in jet_cells])        


        dataList.append([met,nj,pTj1,dphi_min,metParton,njHLT,pTHLT,emf_min,
                         bb1dphi,bb1dr,bb2dphi,bb2dr])
    f.Close()
    
    df = pd.DataFrame(columns=cols,data=dataList)
            

    # ### Get binning from Tobias curves
    tobiasCurves = {L1metStr : '../B2TF-Tobias/results/l1_met_new.csv',
                    L1njStr : '../B2TF-Tobias/results/l1_nj_new.csv',
                    L1pTj1Str : '../B2TF-Tobias/results/l1_et_new.csv',
                    L1dPhi : '../B2TF-Tobias/results/l1_del_phi_new.csv',
                    L1metPartonStr : '../B2TF-Tobias/results/l1_met_new.csv',

                    HLTnjStr : '../B2TF-Tobias/results/offline_nj.csv',
                    HLTpTj1Str : '../B2TF-Tobias/results/offline_pt.csv',
                    }

    # Define default (fallback) binning
    binsDict = {L1metStr : np.arange(0.,300.,20.), 
                L1njStr : np.arange(-0.5,5.,1.),
                L1pTj1Str : np.arange(0.,100.,5.),
                L1dPhi : np.linspace(0.0,5.0,20),
                L1metPartonStr : np.arange(0.,300.,20.), 
                HLTemfStr : np.array([-1.1,-0.9,-0.01,0.01,0.02,0.03,0.04,
                                        0.05,0.06,0.07,0.08,0.09,0.1,0.2,
                                        0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.2]),
                bb1dPhi : np.linspace(0.,np.pi,100),
                bb1dR : np.linspace(0.,10.0,100),
                bb2dPhi : np.linspace(0.,np.pi,100),
                bb2dR : np.linspace(0.,10.0,100)
                }

    # ### Save histograms
    # Make sure the label can be used as a filename:
    label_norm =  re.sub('[^A-Za-z0-9]+', '_', label.strip())        
    for var in cols:
        if not var in csvFiles:
            continue
        
        csvFile = csvFiles[var].replace('.csv','_%s.csv' %label)
        f = open(csvFile,'w')
        f.write(f'# Values for variable: {var}\n')
        if var in tobiasCurves:
                tobiasData = np.genfromtxt(tobiasCurves[var],delimiter=',',names=True)
                bin_centers = tobiasData['bin_center']
                bins = list(tobiasData['bin_left_edge'])
                bins.append(tobiasData['bin_right_edge'][-1])
        else:
            v_max, v_min = df[var].max(),df[var].min()
            if var in binsDict:
                bins = binsDict[var]
            else:
                bins = np.linspace(v_min,v_max,25)
            # print(var,bins)
            bin_centers = 0.5*(bins[1:] + bins[:-1])

        dataDict = {'bin_center' : bin_centers, 
                'bin_left_edge' : bins[:-1],  'bin_right_edge' : bins[1:]}
        dataDF = pd.DataFrame(data=dataDict,index=bin_centers)
        dataDF['bin_count'] = pd.cut(df[var],bins=bins,
                                        right=True,include_lowest=True,
                                        retbins=False,ordered=False,
                                        labels=bin_centers).value_counts()
       
        nTot = float(np.sum(dataDF['bin_count']))
        dataDF['normalized_count'] = dataDF['bin_count']/nTot
        dataDF['error'] = np.sqrt(dataDF['bin_count'])/nTot
        dataDF.drop(columns=['bin_count'],inplace=True)

        f.write(f'# Input file for {label_norm} : {file} ({int(nTot)} MC events)\n')
        f.close()
        print(f'Saving {csvFile}')
        dataDF.to_csv(csvFile,index=False,float_format='%1.4e',mode='a')
