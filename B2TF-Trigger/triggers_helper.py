#!/usr/bin/env python3

from typing import List,Union
import numpy as np
from xml.etree import ElementTree as ET



def B2TF_L1(metOnTime, jetsDelayedL1 : List) -> dict:

    l1_cuflow = {
                    'L1: 40 GeV < MET(N-1) < 100 GeV' : 0,
                    'L1: 40 GeV < PT Jet1(N)' : 0,
                    'L1: DPhi(Jet(N),MET(N-1)) < 1.0' : 0,
                }


    if not (40.0 < metOnTime.MET < 100.0):
        return l1_cuflow

    l1_cuflow['L1: 40 GeV < MET(N-1) < 100 GeV'] += 1

    # Keep only jets passing the pt cut
    jetsDelayedL1 = [j for j in jetsDelayedL1 if j.PT > 20.0]
    # Keep only jets passing the eta cut
    jetsDelayedL1 = [j for j in jetsDelayedL1 if abs(j.Eta) < 3.2]   
    # Sort jets by highest pT
    jetsDelayedL1 = sorted(jetsDelayedL1, 
                         key = lambda j: j.PT, reverse=True)

    if (not jetsDelayedL1) or  (jetsDelayedL1[0].PT < 40.0):
        return l1_cuflow
    
    l1_cuflow['L1: 40 GeV < PT Jet1(N)'] += 1

    dphi_min = 2*np.pi
    for j in jetsDelayedL1[:6]:
        dphi = np.abs(j.Phi-metOnTime.Phi)
        if dphi > np.pi:
            dphi = 2*np.pi - dphi
        dphi_min = min(dphi,dphi_min)

    if dphi_min > 1.0:
        return l1_cuflow
    
    l1_cuflow['L1: DPhi(Jet(N),MET(N-1)) < 1.0'] += 1

    return l1_cuflow


def CalRatioLowET_L1(jetsL1 : List) -> dict:

    l1_cuflow = {
                    'L1: Eta Jet < 2.5' : 0,
                    'L1-LowET: jet(HCAL) > 30 GeV and jet(ECAL) < 3 GeV' : 0,
                }

    if not jetsL1:
        return l1_cuflow
    jetsL1 = [j for j in jetsL1 if (abs(j.Eta) < 2.5)]
    jetsL1 = sorted(jetsL1, key=lambda j: j.PT, reverse=True)
    if not jetsL1:
        return l1_cuflow
    
    l1_cuflow['L1: Eta Jet < 2.5'] += 1
            
    # For the remaining jets, check if there is one jet
    # Remove jets with large ECAL deposits
    jets_disp = []
    for j in jetsL1:
        jetE = j.PT*np.cosh(j.Eta)
        Eem = jetE/(1.0+j.EhadOverEem)
        Ehad = jetE-Eem
        if (Ehad > 30.0 and Eem < 3.):
            jets_disp.append(j)
    
    if not jets_disp:
        return l1_cuflow
    
    l1_cuflow['L1-LowET: jet(HCAL) > 30 GeV and jet(ECAL) < 3 GeV'] += 1

    return l1_cuflow

def CalRatioHighET_L1(jetsL1 : List, pTmin : float = 60.0) -> dict:

    l1_cuflow = {
                    'L1: Eta Jet < 2.5' : 0,
                    f'L1-HighET: jet PT > {pTmin:1.0f} GeV' : 0                    
                }

    if not jetsL1:
        return l1_cuflow
    

    jetsL1 = [j for j in jetsL1 if (abs(j.Eta) < 2.5)]
    jetsL1 = sorted(jetsL1, key=lambda j: j.PT, reverse=True)
    if not jetsL1:
        return l1_cuflow
    
    l1_cuflow['L1: Eta Jet < 2.5'] += 1
            
    # For the remaining jets, check if there is one jet satisfying the minimum pT requirement    
    if jetsL1[0].PT < pTmin:
        return l1_cuflow
    
    l1_cuflow[f'L1-HighET: jet PT > {pTmin:1.0f} GeV'] += 1

    return l1_cuflow

def CalRatio_HLT(jetsHLT : List, tracks : List) -> dict:

    hlt_cuflow = {
                    'HLT: Eta Jet < 2.5 PT Jet > 20 GeV' : 0,
                    'HLT: Jet(N) EMF < 0.06' : 0,
                    'HLT: DR(Tracks(N),Jet(N)) > 0.2' : 0
                }

    jetsHLT = [j for j in jetsHLT if (abs(j.Eta) < 2.5 and j.PT > 20.0)]
    if not jetsHLT:
        return hlt_cuflow
    
    hlt_cuflow['HLT: Eta Jet < 2.5 PT Jet > 20 GeV'] += 1
        
    jets_disp = []
    for j in jetsHLT:
        EMF = 1.0/(1.0+j.EhadOverEem)
        if EMF < 0.06:
            jets_disp.append(j)

    if not jets_disp:
        return hlt_cuflow
    
    hlt_cuflow['HLT: Jet(N) EMF < 0.06'] += 1
    

    # Finally for the remaining jets remove jets with tracks close
    # to it and belonging to the delayed event    
    jets_clean = []
    for j in jets_disp:
        dRmin = 100.0
        for track in tracks:
            dphi = np.abs((j.Phi-track.Phi))
            if dphi > np.pi:
                dphi = 2*np.pi - dphi
            deta = (j.Eta-track.Eta)
            dR = np.sqrt(deta**2 + dphi**2)
            dRmin = min(dR,dRmin)
        if dRmin > 0.2:
            jets_clean.append(j)

    if not jets_clean:
        return hlt_cuflow
    
    hlt_cuflow['HLT: DR(Tracks(N),Jet(N)) > 0.2'] += 1
    
    return hlt_cuflow

