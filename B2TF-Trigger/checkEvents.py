#!/usr/bin/env python3

import os, time, sys
sys.path.append('../')
import numpy as np
import pyhepmc
import logging
import fastjet
import awkward as ak
import vector
vector.register_awkward()

delphesDir = os.path.abspath("../DelphesLLP")
os.environ['ROOT_INCLUDE_PATH'] = os.path.join(delphesDir,"external")

import ROOT


ROOT.gSystem.Load(os.path.join(delphesDir,"libDelphes.so"))

ROOT.gInterpreter.Declare('#include "classes/SortableObject.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')



inputFile = '../pp2chi0chi0_minimalH_scan_jets_fix/Events/run_06/ddmH_mS_2000_m1_979_dm_90_delphes_events.root'
# inputFile = '../DelphesLLP/test.root'

f = ROOT.TFile(inputFile,'read')
tree = f.Get("Delphes")
tree.GetEntry(10)

storeParticles = []

pS = np.array([0.0,0.0,0.0,0.0])
print('############ Case A: chi1 on-time and anti-chi1 delayed ###################')
print('############ ON-TIME:')
print('llp:')
for p in tree.llpParticlesA:
    print(p.PID,p.Status)
    print(p.Px,p.Py,p.Pz,p.E)

# print('mom:')
# for p in tree.llpMothersA:
    # print(p.PID,p.Status)    

print('direct:')
pVis = np.array([0.0,0.0,0.0,0.0])
pInv = np.array([0.0,0.0,0.0,0.0])
for p in tree.llpDirectDaughtersA:
    print(p.PID,p.Status)   
    if abs(p.PID) < 10000:
        pVis += np.array([p.Px,p.Py,p.Pz,p.E])
    else:
        pInv += np.array([p.Px,p.Py,p.Pz,p.E])
pTvis = np.sqrt(pVis[0]**2+pVis[1]**2)
print(f'  (pT visible = {pTvis})')
pTot = pVis+pInv
print(f'  (pTot = {pTot})')
mTot = np.sqrt(pTot[3]**2-np.dot(pTot[:3],pTot[:3]))
print(f'  (mTot  = {mTot})')
pS = pS + pTot

print('final:')
pVis = np.array([0.0,0.0,0.0,0.0])
for p in tree.llpFinalA:
    if abs(p.PID) > 10000 or p.Status != 1:
        print(p.PID,p.Status)
    else:
        pVis += np.array([p.Px,p.Py,p.Pz,p.E])
print('...')
pTvis = np.sqrt(pVis[0]**2+pVis[1]**2)
print(f'  (pT visible = {pTvis})')

print('others:')
for p in tree.llpOthersA:
    if abs(p.PID) > 10000 or p.Status != 1:
        print(p.PID,p.Status)
print('...')

print('############ DELAYED:')
print('llp:')
for p in tree.llpParticlesB:
    print(p.PID,p.Status)
    print(p.Px,p.Py,p.Pz,p.E)

# print('mom:')
# for p in tree.llpMothersB:
#     print(p.PID,p.Status)    

print('direct:')
pVis = np.array([0.0,0.0,0.0,0.0])
pInv = np.array([0.0,0.0,0.0,0.0])
for p in tree.llpDirectDaughtersB:
    print(p.PID,p.Status)   
    if abs(p.PID) < 10000:
        pVis += np.array([p.Px,p.Py,p.Pz,p.E])
    else:
        pInv += np.array([p.Px,p.Py,p.Pz,p.E])
pTvis = np.sqrt(pVis[0]**2+pVis[1]**2)
print(f'  (pT visible = {pTvis})')
pTot = pVis+pInv
print(f'  (pTot = {pTot})')
mTot = np.sqrt(pTot[3]**2-np.dot(pTot[:3],pTot[:3]))
print(f'  (mTot  = {mTot})')
pS = pS + pTot

print('final:')
pVis = np.array([0.0,0.0,0.0,0.0])
for p in tree.llpFinalB:
    if abs(p.PID) > 10000 or p.Status != 1:
        print(p.PID,p.Status)
    else:
        if abs(p.PID) not in [12,14,16,13]:
            storeParticles.append([p.Px,p.Py,p.Pz,p.E])
        pVis += np.array([p.Px,p.Py,p.Pz,p.E])
print('...')
pTvis = np.sqrt(pVis[0]**2+pVis[1]**2)
print(f'  (pT visible = {pTvis})')

print('others:')
for p in tree.llpOthersB:
    if abs(p.PID) > 10000 or p.Status != 1:
        print(p.PID,p.Status)
print('...')

pTS = np.sqrt(pS[0]**2 + pS[1]**2)
print(f'  (pS = {pS}, pTS = {pTS})')
mS = np.sqrt(pS[3]**2-np.dot(pS[:3],pS[:3]))
print(f'  (MS  = {mS})')


print('####### ALL JETS ################')
for j in tree.GenJet:
    print(j.PT,len(j.Constituents))


print('####### ON-TIME JETS ################')
for j in tree.GenJetAOnTime:
    print(j.PT,len(j.Constituents))    



print('####### Delayed JETS ################')
for j in tree.GenJetADelayed:
    print(j.PT,len(j.Constituents))        




print('####### JETS from Delayed decays ################')

jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
array0 = []
for p in storeParticles:
    array0.append({'px' : p[0], 'py' : p[1], 'pz' : p[2], 'E' : p[3]})
array1 = ak.Array(array0)

cluster = fastjet.ClusterSequence(array1, jetdef)
for c in cluster.inclusive_jets(min_pt=20.0):
    print(c,np.sqrt(c['px']**2 + c['py']**2))


storeParticles = []
pS = np.array([0.0,0.0,0.0,0.0])
print('\n\n\n############ Case B: chi1 delayed and anti-chi1 on-time ###################')
print('############ ON-TIME:')
print('llp:')
for p in tree.llpParticlesB:
    print(p.PID,p.Status)
    print(p.Px,p.Py,p.Pz,p.E)

# print('mom:')
# for p in tree.llpMothersB:
    # print(p.PID,p.Status)    

print('direct:')
pVis = np.array([0.0,0.0,0.0,0.0])
pInv = np.array([0.0,0.0,0.0,0.0])
for p in tree.llpDirectDaughtersB:
    print(p.PID,p.Status)   
    if abs(p.PID) < 10000:
        pVis += np.array([p.Px,p.Py,p.Pz,p.E])
    else:
        pInv += np.array([p.Px,p.Py,p.Pz,p.E])
pTvis = np.sqrt(pVis[0]**2+pVis[1]**2)
print(f'  (pT visible = {pTvis})')
pTot = pVis+pInv
print(f'  (pTot = {pTot})')
mTot = np.sqrt(pTot[3]**2-np.dot(pTot[:3],pTot[:3]))
print(f'  (mTot  = {mTot})')
pS = pS + pTot

print('final:')
pVis = np.array([0.0,0.0,0.0,0.0])
for p in tree.llpFinalB:
    if abs(p.PID) > 10000 or p.Status != 1:
        print(p.PID,p.Status)
    else:
        pVis += np.array([p.Px,p.Py,p.Pz,p.E])
print('...')
pTvis = np.sqrt(pVis[0]**2+pVis[1]**2)
print(f'  (pT visible = {pTvis})')

print('others:')
for p in tree.llpOthersB:
    if abs(p.PID) > 10000 or p.Status != 1:
        print(p.PID,p.Status)
print('...')

print('############ DELAYED:')
print('llp:')
for p in tree.llpParticlesA:
    print(p.PID,p.Status)
    print(p.Px,p.Py,p.Pz,p.E)

# print('mom:')
# for p in tree.llpMothersA:
    # print(p.PID,p.Status)    

print('direct:')
pVis = np.array([0.0,0.0,0.0,0.0])
pInv = np.array([0.0,0.0,0.0,0.0])
for p in tree.llpDirectDaughtersA:
    print(p.PID,p.Status)   
    if abs(p.PID) < 10000:
        pVis += np.array([p.Px,p.Py,p.Pz,p.E])
    else:
        pInv += np.array([p.Px,p.Py,p.Pz,p.E])
pTvis = np.sqrt(pVis[0]**2+pVis[1]**2)
print(f'  (pT visible = {pTvis})')
pTot = pVis+pInv
print(f'  (pTot = {pTot})')
mTot = np.sqrt(pTot[3]**2-np.dot(pTot[:3],pTot[:3]))
print(f'  (mTot  = {mTot})')
pS = pS + pTot

print('final:')
pVis = np.array([0.0,0.0,0.0,0.0])
for p in tree.llpFinalA:
    if abs(p.PID) > 10000 or p.Status != 1:
        print(p.PID,p.Status)
    else:
        if abs(p.PID) not in [12,14,16,13]:
            storeParticles.append([p.Px,p.Py,p.Pz,p.E])
        pVis += np.array([p.Px,p.Py,p.Pz,p.E])
print('...')
pTvis = np.sqrt(pVis[0]**2+pVis[1]**2)
print(f'  (pT visible = {pTvis})')

print('others:')
for p in tree.llpOthersA:
    if abs(p.PID) > 10000 or p.Status != 1:
        print(p.PID,p.Status)
print('...')

pTS = np.sqrt(pS[0]**2 + pS[1]**2)
print(f'  (pS = {pS}, pTS = {pTS})')
mS = np.sqrt(pS[3]**2-np.dot(pS[:3],pS[:3]))
print(f'  (MS  = {mS})')


print('####### ALL JETS ################')
for j in tree.GenJet:
    print(j.PT,len(j.Constituents))


print('####### ON-TIME JETS ################')
for j in tree.GenJetBOnTime:
    print(j.PT,len(j.Constituents))    



print('####### Delayed JETS ################')
for j in tree.GenJetBDelayed:
    print(j.PT,len(j.Constituents))        
f.Close()

print('####### JETS from Delayed decays ################')

jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
array0 = []
for p in storeParticles:
    array0.append({'px' : p[0], 'py' : p[1], 'pz' : p[2], 'E' : p[3]})
array1 = ak.Array(array0)

cluster = fastjet.ClusterSequence(array1, jetdef)
for c in cluster.inclusive_jets(min_pt=20.0):
    print(c,np.sqrt(c['px']**2 + c['py']**2))


f.Close()