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
N1betaStr = r'$\beta$ (N-1)'
N2betaStr = r'$\beta$ (N)'
N1rhoStr = r'$\rho_{\rm dec}$ (N-1) (m)'
N2rhoStr = r'$\rho_{\rm dec}$ (N) (m)'
N1zStr = r'$z_{\rm dec}$ (N-1) (m)'
N2zStr = r'$z_{\rm dec}$ (N) (m)'
N1toutStr = r'$t_{\rm readout}$ (N-1) (ns)'
N2toutStr = r'$t_{\rm readout}$ (N) (ns)'

cols = [N1betaStr,N2betaStr,N1rhoStr,N2rhoStr,N1zStr,N2zStr,N1toutStr,N2toutStr]


if __name__ == "__main__":
    

    import argparse
    import progressbar as P
    from tqdm import tqdm
    ap = argparse.ArgumentParser( description=
            "Read the output of a Delphes ROOT file and save the distributions to CSV files." )
    ap.add_argument('-f', '--inputfile',
            help='path to the ROOT file.',required=True)
    ap.add_argument('-l', '--label',
            help='label to the output file.',required=False,default="")    
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
        
    c_light = 2.99792458e8

    llpData = []
    for ievt in tqdm(range(nevts)):
        tree.GetEntry(ievt)

        # Get parton level MET and b-bar angular separation
        llps = list(tree.llpParticles)
        llpList = []
        for illp,llp in enumerate(llps):            
            for d in tree.llpDirectDaughters:
                if illp == d.M1:
                    x = np.array([d.X,d.Y,d.Z,d.T])
                    l = np.linalg.norm(x[:3])*1e-3
                    t_readout = x[-1]-l/c_light
                    t_readout = t_readout*1e9
                    break
            
            if t_readout < 10:
                eventRecord = 'N-1'
            elif 25 < t_readout < 35:
                eventRecord = 'N'
            else:
                continue
            beta = np.sqrt(llp.Px**2 + llp.Py**2 + llp.Pz**2)/llp.E
            rho = np.linalg.norm(x[:2])*1e-3
            z = x[2]*1e-3
            llpDict = {'eventRecord' : eventRecord, r'$\beta$' : beta,
                       r'$\rho_{\rm dec}$ (m)' : rho, r'$z$ (m)' : z, 
                       r'$t_{\rm readout}$ (ns)' : t_readout, 'event' : ievt}
        
            llpList.append(llpDict)
        llpList = sorted(llpList, key = lambda llpD: llpD[r'$t_{\rm readout}$ (ns)'])
        for illp,llpD in enumerate(llpList):
            llpD['illp'] = illp
            llpData.append(llpD)
    f.Close()
    
    dataDF = pd.DataFrame(data=llpData).set_index(['event','eventRecord','illp'])
    
    outFile = file.replace('.root','')+label+'.pcl'
    print(f'Saving {outFile}')
    dataDF.to_pickle(outFile)
