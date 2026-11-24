#!/usr/bin/env python3
from typing import List
import numpy as np
from xml.etree import ElementTree as ET



def selectionOnTime(metOnTime, jetsOnTime : List) -> dict:


    l1_cuflow = {
                'L1(N): 200 GeV < MET(N)' : 0,
                'L1(N): 150 GeV < pTj1(N)' : 0,
                'L1(N): Dphi(jet_MET) > 0.4/0.6' : 0
                }


    if not (metOnTime.MET > 200.0):
        return l1_cuflow

    l1_cuflow['L1(N): 200 GeV < MET(N)'] += 1
    
    jets = [j for j in jetsOnTime if (abs(j.Eta) < 2.4) and (j.PT > 25.0)]
    jets = sorted(jetsOnTime, key=lambda j: j.PT, reverse=True)

    if not jets or not (jets[0].PT > 150.0):
        return l1_cuflow
    l1_cuflow['L1(N): 150 GeV < pTj1(N)'] += 1

    
    for j in jets[:4]:
        dphi_min = 2*np.pi
        dphi = np.abs(j.Phi-metOnTime.Phi)
        if dphi > np.pi:
            dphi = 2*np.pi - dphi
        dphi_min = min(dphi,dphi_min)
        if metOnTime.MET <= 250.0 and dphi_min < 0.6:
            return l1_cuflow
        elif metOnTime.MET > 250.0 and dphi_min < 0.4:
            return l1_cuflow
    
    l1_cuflow['L1(N): Dphi(jet_MET) > 0.4/0.6'] += 1

    return l1_cuflow


@np.vectorize
def in_ms(z,rho):

    if 7.0 < abs(z) < 13.0 and rho < 7.5:
        return True
    elif 3.5 < rho < 7.5 and abs(z) < 7.0:
        return True
    else:
        return False

@np.vectorize
def in_N(t):
    
    if t < 10.0:
        return True
    else:
        return False

@np.vectorize
def in_Np1(t):

    if 25.0 < t < 50.0:
        return True
    else:
        return False

@np.vectorize
def in_hcal(z,rho):

    if 4.3 < abs(z) < 5.0 and rho < 3.5:
        return True
    elif 2.0 < rho < 3.5 and abs(z) < 4.3:
        return True
    else:
        return False

def selectionDelayed(llpDecays : List, pdgInvisible=4000022) -> dict:

    c_light = 3e8
    cutflow = {
                '(N+1): 25 ns < T < 50 ns' : 0,
                '(N+1): decay in MS' : 0,
                '(N+1): nRoI >= 3' : 0
                }


    # Select decays in N+1:
    decaysNp1 = []
    for d in llpDecays:
        l_m = np.linalg.norm([d.X,d.Y,d.Z])*1e-3
        t = d.T
        t_readout = (t - l_m/c_light)*1e9 # readout time in ns
        if in_Np1(t_readout):
            decaysNp1.append(d)

    if not decaysNp1:
        return cutflow
    cutflow['(N+1): 25 ns < T < 50 ns'] += 1

    # Select decays within MS fiducial volume:
    decaysMS = []
    for d in decaysNp1:
        rho = np.linalg.norm([d.X,d.Y])*1e-3
        z = d.Z*1e-3
        if in_ms(z,rho):
            decaysMS.append(d)
    
    if not decaysMS:
        return cutflow
    cutflow['(N+1): decay in MS'] += 1

    # Select minimum visible energy in MS:
    nRoI = 0.0
    for d in decaysMS:
        if abs(d.PID) == pdgInvisible:
            continue
        nRoI += (d.E/4.0)

    if nRoI < 3.0:
        return cutflow 
    cutflow['(N+1): nRoI >= 3'] += 1

    return cutflow

