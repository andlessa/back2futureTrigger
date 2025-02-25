import sys
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import scipy
import time
import readMapNew as rmN
import tqdm
import mplhep as hep
import lhe_parser as lhe
import hepmc_parser as hepmc
import uproot
import random
import os
random.seed(123)
hep.style.use("ATLAS")

# set cst
c = 3e8# Light velocity in m/s

#All plots are made with 10 000 events, if you want to try with other numbers of events, you will have to change the number of events in lines 285-314-551-580 for the calculation of the efficiency.


def getLLPdecays(event,llpPdgList=[35]):

    vertexList = []
    for vertex in event.vertex.values():
        if len(vertex.incoming) != 1:
            continue
        if len(vertex.outgoing) == 1:
            continue
        part_in = vertex.incoming[0]
        # Check if incoming particle appear in pdgList
        if abs(part_in.pdg) not in llpPdgList:
            continue

        part_out_list = list(vertex.outcoming)
        # Skip vertices which could correspond to FSR (initial state appears in final state)
        if any(p.pdg == part_in.pdg for p in part_out_list):
            continue

        vertexList.append(vertex)
    
    return vertexList

#########################################################################################
#Parsing the hepmc file from the hadronization of the MG outputs to recover the data from the process.
#########################################################################################

def parsing_hepmc_new(events,llps=[35],invisibles=[12,14,16]):

    px_TOT = []
    py_TOT = []
    pz_TOT = []
    E_TOT = []
    mass_TOT = []
    pdg_TOT = []

    for event in events:
        vertexList = getLLPdecays(event,llpPdgList=llps)
        px_event = []
        py_event = []
        pz_event = []
        E_event = []
        mass_event = []
        pdg_event = []
        if len(vertexList) != 2:
            raise ValueError(f"{len(vertexList)} LLP decay vertices found (can only deal with 2 decay vertices)")
        for vertex in vertexList:
            visible_particles = [p for p in vertex.outcoming 
                                 if abs(p.pdg) not in invisibles]
            if len(visible_particles) != 2:
                raise ValueError(f"{len(visible_particles)} visible particles found in LLP decay (can only deal with 2 particles)")
            px_visible = sum([p.px for p in visible_particles])
            py_visible = sum([p.py for p in visible_particles])
            pz_visible = sum([p.pz for p in visible_particles])
            E_visible =  sum([p.e for p in visible_particles])
            mass_visible = np.sqrt(E_visible**2 - px_visible**2 - py_visible**2 - pz_visible**2)
            px_event.append(px_visible)
            py_event.append(py_visible)
            pz_event.append(pz_visible)
            E_event.append(E_visible)
            mass_event.append(mass_visible)
            pdg_event.append([p.pdg for p in visible_particles])
        px_TOT.append(px_event)
        py_TOT.append(py_event)
        pz_TOT.append(pz_event)
        E_TOT.append(E_event)
        mass_TOT.append(mass_event)
        pdg_TOT.append(pdg_event)

    return px_TOT, py_TOT, pz_TOT, E_TOT, mass_TOT, pdg_TOT


def parsing_hepmc(events):

    px_TOT = []
    py_TOT = []
    pz_TOT = []
    E_TOT = []
    mass_TOT = []
    pdg_TOT = []

    for ie , event in enumerate(events):
        count=0
        for id, vertex in event.vertex.items():
            if [p.pdg for p in vertex.incoming] == [25] and [p.pdg for p in vertex.outcoming] == [35, 35]: # PDGID 25 = Higgs, PDGID 35 Dark Higgs
                px_TOT.append(list(p.px for p in vertex.outcoming)) # recover the x momenta in GeV
                py_TOT.append(list(p.py for p in vertex.outcoming)) # recover the y momenta in GeV
                pz_TOT.append(list(p.pz for p in vertex.outcoming)) # recover the z momenta in GeV
                E_TOT.append(list(p.E for p in vertex.outcoming)) # recover the Energy in GeV
                mass_TOT.append(list(p.mass for p in vertex.outcoming)) # recover the mass in GeV

            if [p.pdg for p in vertex.incoming] == [35]: # PDGID 35 Dark Higgs
                pdg_TOT.append((list(p.pdg for p in vertex.outcoming))) # recover the PDG ID of the particle produced

                count = count+1
                if count==2: ##
                    break
                    pass

    return px_TOT, py_TOT, pz_TOT, E_TOT, mass_TOT, pdg_TOT

#########################################################################################
#The data recovered are list of list, we need to convert them into one list to be able to separate the contribution of each LLP.
#########################################################################################

def conversion_one_list(px_TOT, py_TOT, pz_TOT, E_TOT, mass_TOT, pdg_TOT):

    px_tot = []
    for i in range(len(px_TOT)):
        for y in range(len(px_TOT[i])):
            px_tot.append(px_TOT[i][y])

    py_tot = []
    for i in range(len(py_TOT)):
        for y in range(len(py_TOT[i])):
            py_tot.append(py_TOT[i][y])

    pz_tot = []
    for i in range(len(pz_TOT)):
        for y in range(len(pz_TOT[i])):
            pz_tot.append(pz_TOT[i][y])

    E_tot = []
    for i in range(len(E_TOT)):
        for y in range(len(E_TOT[i])):
            E_tot.append(E_TOT[i][y])

    mass_tot = []
    for i in range(len(mass_TOT)):
        for y in range(len(mass_TOT[i])):
            mass_tot.append(mass_TOT[i][y])

    pdg_tot = []
    for i in range(len(pdg_TOT)):
        for y in range(len(pdg_TOT[i])):
            pdg_tot.append(pdg_TOT[i][y])

    return px_tot, py_tot, pz_tot, E_tot, mass_tot, pdg_tot

#########################################################################################
# Recovering the data from each LLP (px,py,pz,E,mass,PDG ID).
#########################################################################################

def recover(px_tot, py_tot, pz_tot, E_tot, mass_tot,pdg_tot):

    px_DH1 = []
    px_DH2 = []

    py_DH1 = []
    py_DH2 = []

    pz_DH1 = []
    pz_DH2 = []

    E_DH1 = []
    E_DH2 = []

    mass_DH1 = []
    mass_DH2 = []

    pdg_tot_DH1 = []
    pdg_tot_DH2 = []

    for i in range(0, len(px_tot),2): # in the list, each even value is for DH1
        px_DH1.append(px_tot[i])
        py_DH1.append(py_tot[i])
        pz_DH1.append(pz_tot[i])
        E_DH1.append(E_tot[i])
        mass_DH1.append(mass_tot[i])

    for i in range(1, len(px_tot),2): # in the list, each odd value is for DH2
        px_DH2.append(px_tot[i])
        py_DH2.append(py_tot[i])
        pz_DH2.append(pz_tot[i])
        E_DH2.append(E_tot[i])
        mass_DH2.append(mass_tot[i])

    for i in range(0, len(pdg_tot),4): # recover the PDG ID of the particle produced by decay of DH1
        pdg_tot_DH1.append(pdg_tot[i])

    for i in range(2, len(pdg_tot),4): # recover the PDG ID of the particle produced by decay of DH2
        pdg_tot_DH2.append(pdg_tot[i])

    # Convert all lists into arrays

    px_DH1 = np.array(px_DH1)/c # GeV/c
    py_DH1 = np.array(py_DH1)/c
    pz_DH1 = np.array(pz_DH1)/c
    E_DH1 = np.array(E_DH1)
    mass_DH1 = np.array(mass_DH1)
    pdg_tot_DH1 = np.array(pdg_tot_DH1)

    px_DH2 = np.array(px_DH2)/c
    py_DH2 = np.array(py_DH2)/c
    pz_DH2 = np.array(pz_DH2)/c
    E_DH2 = np.array(E_DH2)
    mass_DH2 = np.array(mass_DH2)
    pdg_tot_DH2 = np.array(pdg_tot_DH2)

    return px_DH1, px_DH2, py_DH1, py_DH2, pz_DH1, pz_DH2, pdg_tot_DH1, pdg_tot_DH2, E_DH1, E_DH2, mass_DH1, mass_DH2

#########################################################################################
# Computation of the kinematics variable for LLP1 (velocities, beta, gamma, pT the transverse momenta, eta the pseudo-rapidity).
#########################################################################################

def kinematics_DH(px_DH1, py_DH1, pz_DH1, E_DH1):

    vx_DH1 = (px_DH1*c**2)/E_DH1 #compute the velocities in each direction
    vy_DH1 = (py_DH1*c**2)/E_DH1
    vz_DH1 = (pz_DH1*c**2)/E_DH1
    beta_DH1 = np.sqrt(vx_DH1**2 + vy_DH1**2 + vz_DH1**2)/c # compute beta
    gamma_DH1 = 1/(np.sqrt(1-beta_DH1**2)) # compute gamma

    pT_DH1 = np.sqrt(px_DH1**2 + py_DH1**2)*c # compute the transverse momenta
    eta_DH1 = np.arctanh(pz_DH1/(np.sqrt(px_DH1**2 + py_DH1**2 + pz_DH1**2))) # compute the pseudorapidity

    return beta_DH1, gamma_DH1, pT_DH1, eta_DH1


#########################################################################################
# lifetime function.
#########################################################################################

def lifetime(avgtau = 4.3):
    import math
    avgtau = avgtau / c
    t = random.random()
    return -1.0 * avgtau * math.log(t)

#########################################################################################
# Decay lenght computation for LLP1.
#########################################################################################

def decaylenghtDH(px_DH1, py_DH1, pz_DH1, E_DH1, gamma_DH1, tauN):

    Lx_tot_DH1 = []
    Ly_tot_DH1 = []
    Lz_tot_DH1 = []
    Lxy_tot_DH1 = []

    for ctau in range(len(tauN)):

        Lx_DH1 = []
        Ly_DH1 = []
        Lz_DH1 = []
        Lxy_DH1 = []

        for i in range(len(gamma_DH1)):
            lt = lifetime(tauN[ctau]) # set mean lifetime
            Lx_DH1.append((px_DH1[i]/E_DH1[i])*c**2 * lt * gamma_DH1[i]) # compute the decay lenght in x,y,z
            Ly_DH1.append((py_DH1[i]/E_DH1[i])*c**2 * lt * gamma_DH1[i])
            Lz_DH1.append((abs(pz_DH1[i])/E_DH1[i])*c**2 * lt  * gamma_DH1[i] )
            Lxy_DH1.append(np.sqrt((Lx_DH1[i])**2 + (Ly_DH1[i])**2)) # compte the transverse decay lenght

        Lx_tot_DH1.append(Lx_DH1)
        Ly_tot_DH1.append(Ly_DH1)
        Lz_tot_DH1.append(Lz_DH1)
        Lxy_tot_DH1.append(Lxy_DH1)
    return Lxy_tot_DH1, Lz_tot_DH1

#########################################################################################
# Computation of the efficiency with the map from the data obtained with MG+Pythia8 for the high-ET samples (mH >= 400GeV).
#########################################################################################

def eff_map_High(pT_DH1, eta_DH1, Lxy_tot_DH1, Lz_tot_DH1, pdg_tot_DH1, pT_DH2, eta_DH2, Lxy_tot_DH2, Lz_tot_DH2, pdg_tot_DH2, tauN, nevent, mass_phi, mass_s):

    eff_highETX = []

    for index in tqdm.tqdm(range(len(tauN))):
        queryMapResult = []
        for iEvent in range(len(pT_DH1)):
            queryMapResult.append(rmN.queryMapFromKinematics(pT_DH1[iEvent],
                                                            eta_DH1[iEvent],
                                                            Lxy_tot_DH1[index][iEvent],
                                                            Lz_tot_DH1[index][iEvent],
                                                            abs(pdg_tot_DH1[iEvent]),
                                                            pT_DH2[iEvent],
                                                            eta_DH2[iEvent],
                                                            Lxy_tot_DH2[index][iEvent],
                                                            Lz_tot_DH2[index][iEvent],
                                                            abs(pdg_tot_DH2[iEvent]),
                                                            selection = "high-ET"))
        eff_highETX.append(sum(queryMapResult))
    queryMapResult = np.array(queryMapResult) #convertion into array
    eff_highETX = np.array(eff_highETX) #convertion into array
    eff_highETX = eff_highETX/nevent #efficiency/(nbr of event)

    return eff_highETX

#########################################################################################
# Computation of the efficiency with the map from the data obtained with MG+Pythia8 for the low-ET samples (mH <= 400GeV).
#########################################################################################

def eff_map_Low(pT_DH1, eta_DH1, Lxy_tot_DH1, Lz_tot_DH1, pdg_tot_DH1, pT_DH2, eta_DH2, Lxy_tot_DH2, Lz_tot_DH2, pdg_tot_DH2, tauN,nevent, mass_phi, mass_s):

    eff_lowETX = []

    for index in tqdm.tqdm(range(len(tauN))):
        queryMapResult = []
        for iEvent in range(len(pT_DH1)):
            queryMapResult.append(rmN.queryMapFromKinematics(pT_DH1[iEvent],
                                                            eta_DH1[iEvent],
                                                            Lxy_tot_DH1[index][iEvent],
                                                            Lz_tot_DH1[index][iEvent],
                                                            abs(pdg_tot_DH1[iEvent]),
                                                            pT_DH2[iEvent],
                                                            eta_DH2[iEvent],
                                                            Lxy_tot_DH2[index][iEvent],
                                                            Lz_tot_DH2[index][iEvent],
                                                            abs(pdg_tot_DH2[iEvent]),
                                                            selection = "low-ET"))
        eff_lowETX.append(sum(queryMapResult))
    queryMapResult = np.array(queryMapResult) #convertion into array
    eff_lowETX = np.array(eff_lowETX) #convertion into array
    eff_lowETX = eff_lowETX/nevent #efficiency/(nbr of event)

    return eff_lowETX


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

#########################################################################################
#Parsing the lhe file from the MG output to recover the data from the process.
#########################################################################################


def parsing_LHE(MG_events):
    px = []
    py = []
    pz = []
    pdg = []
    E = []
    MASS = []
    for event in MG_events:
        for particle in event:
            pdg.append(particle.pdg)
            if particle.pdg == 35: # PDGID 35 Dark Higgs
                px.append(particle.px)
                py.append(particle.py)
                pz.append(particle.pz)
                E.append(particle.E)
                MASS.append(particle.mass)

    px = np.array(px)/c # GeV/c
    py = np.array(py)/c
    pz = np.array(pz)/c

    return px, py, pz, pdg, E, MASS


#########################################################################################
# Recovering the data from MG (LHE) (PDG ID, px,py,pz,E,mass).
#########################################################################################

def recover_MG_DH1(px, py, pz, E, MASS, pdg):

    MG_pdg_DH1_1 = []
    for i in range(5,len(pdg),9):
        MG_pdg_DH1_1.append(pdg[i]) #List with the PDG ID of the particle produced by the decay of the LLP1

    MG_E_DH1 = []
    for i in range(0,len(px),2):
        MG_E_DH1.append(E[i]) #List with the energy of the LLP1

    MG_px_DH1 = []
    for i in range(0,len(px),2):
        MG_px_DH1.append(px[i]) #List with x momenta from LLP1

    MG_py_DH1 = []
    for i in range(0,len(px),2):
        MG_py_DH1.append(py[i]) #List with y momenta from LLP1

    MG_pz_DH1 = []
    for i in range(0,len(px),2):
        MG_pz_DH1.append(pz[i]) #List with z momenta from LLP1

    MG_mass_DH1 = []
    for i in range(0,len(px),2):
        MG_mass_DH1.append(MASS[i]) #List with the mass from LLP1

    MG_px_DH1 = np.array(MG_px_DH1) # convertion into arrays
    MG_py_DH1 = np.array(MG_py_DH1)
    MG_pz_DH1 = np.array(MG_pz_DH1)
    MG_E_DH1 = np.array(MG_E_DH1)
    MG_mass_DH1 = np.array(MG_mass_DH1)

    return MG_px_DH1, MG_py_DH1,MG_pz_DH1,MG_E_DH1,MG_mass_DH1,MG_pdg_DH1_1

#########################################################################################
# Computation of the kinematics variable for LLP1 (velocities, beta, gamma, pT the transverse momenta, eta the pseudo-rapidity).
#########################################################################################

def kinematics_MG_DH(MG_px_DH1,MG_py_DH1,MG_pz_DH1,MG_E_DH1):

    MG_vx_DH1 = (MG_px_DH1*c**2)/MG_E_DH1 #compute the velocities in each direction
    MG_vy_DH1 = (MG_py_DH1*c**2)/MG_E_DH1
    MG_vz_DH1 = (MG_pz_DH1*c**2)/MG_E_DH1
    MG_beta_DH1 = np.sqrt(MG_vx_DH1**2 + MG_vy_DH1**2 + MG_vz_DH1**2)/c # compute beta
    MG_gamma_DH1 = 1/(np.sqrt(1-MG_beta_DH1**2)) # compute gamma
    MG_pT_DH1 = np.sqrt(MG_px_DH1**2 + MG_py_DH1**2)*c # compute the transverse momenta
    MG_eta_DH1 = np.arctanh(MG_pz_DH1/(np.sqrt(MG_px_DH1**2 + MG_py_DH1**2 + MG_pz_DH1**2))) # compute the pseudorapidity

    return MG_pT_DH1,MG_eta_DH1, MG_gamma_DH1

#########################################################################################
# Recovering the data from LLP2 (PDG ID, px,py,pz,E,mass).
#########################################################################################

def recover_MG_DH2(px, py, pz, E, MASS, pdg):

    MG_pdg_DH2_1 = []
    for i in range(7,len(pdg),9):
        MG_pdg_DH2_1.append(pdg[i]) #List with the PDG ID of the particle produced by the decay of the LLP2

    MG_E_DH2 = []
    for i in range(1,len(px),2):
        MG_E_DH2.append(E[i]) #List with the energy of the LLP2

    MG_px_DH2 = []
    for i in range(1,len(px),2):
        MG_px_DH2.append(px[i]) #List with x momenta from LLP2

    MG_py_DH2 = []
    for i in range(1,len(px),2):
        MG_py_DH2.append(py[i]) #List with y momenta from LLP2

    MG_pz_DH2 = []
    for i in range(1,len(px),2):
        MG_pz_DH2.append(pz[i]) #List with z momenta from LLP2

    MG_mass_DH2 = []
    for i in range(1,len(px),2):
        MG_mass_DH2.append(MASS[i]) #List with the mass from LLP2

    MG_px_DH2 = np.array(MG_px_DH2) # convertion into arrays
    MG_py_DH2 = np.array(MG_py_DH2)
    MG_pz_DH2 = np.array(MG_pz_DH2)
    MG_E_DH2 = np.array(MG_E_DH2)
    MG_mass_DH2 = np.array(MG_mass_DH2)

    return MG_px_DH2, MG_py_DH2,MG_pz_DH2,MG_E_DH2,MG_mass_DH2,MG_pdg_DH2_1


#########################################################################################
# Decay lenght computation for LLP1.
#########################################################################################

def decaylenght_MG_DH(MG_px_DH1, MG_py_DH1, MG_pz_DH1, E_DH1, MG_gamma_DH1, tauN):

    MG_Lx_tot_DH1 = []
    MG_Ly_tot_DH1 = []
    MG_Lz_tot_DH1 = []
    MG_Lxy_tot_DH1 = []

    for ctau in range(len(tauN)):

        MG_Lx_DH1 = []
        MG_Ly_DH1 = []
        MG_Lz_DH1 = []
        MG_Lxy_DH1 = []

        for i in range(len(MG_gamma_DH1)):
            MG_lt = lifetime(tauN[ctau]) # set the mean lifetime
            MG_Lx_DH1.append((MG_px_DH1[i]/E_DH1[i])*c**2 * MG_lt * MG_gamma_DH1[i]) # compute the decay lenght in x,y,z
            MG_Ly_DH1.append((MG_py_DH1[i]/E_DH1[i])*c**2 * MG_lt * MG_gamma_DH1[i])
            MG_Lz_DH1.append((abs(MG_pz_DH1[i])/E_DH1[i])*c**2 * MG_lt  * MG_gamma_DH1[i] )
            MG_Lxy_DH1.append(np.sqrt((MG_Lx_DH1[i])**2 + (MG_Ly_DH1[i])**2)) # compute the transverse decay lenght

        MG_Lx_tot_DH1.append(MG_Lx_DH1) # convertion into arrays
        MG_Ly_tot_DH1.append(MG_Ly_DH1)
        MG_Lz_tot_DH1.append(MG_Lz_DH1)
        MG_Lxy_tot_DH1.append(MG_Lxy_DH1)

    return MG_Lxy_tot_DH1, MG_Lz_tot_DH1


#########################################################################################
# Computation of the efficiency with the map from the data obtained with MG for the high-ET samples (mH <= 400GeV).
#########################################################################################

def eff_map_MG_high(MG_pT_DH1, MG_eta_DH1,MG_Lxy_tot_DH1, MG_Lz_tot_DH1, MG_pdg_DH1_1, MG_pT_DH2, MG_eta_DH2, MG_Lxy_tot_DH2, MG_Lz_tot_DH2, MG_pdg_DH2_1, tauN, nevent, mass_phi, mass_s):

    MG_eff_highETX = []

    for index in tqdm.tqdm(range(len(tauN))):
        MG_queryMapResult = []
        for iEvent in range(len(MG_pT_DH1)):
            MG_queryMapResult.append(rmN.queryMapFromKinematics(MG_pT_DH1[iEvent],
                                                            MG_eta_DH1[iEvent],
                                                            MG_Lxy_tot_DH1[index][iEvent],
                                                            MG_Lz_tot_DH1[index][iEvent],
                                                            abs(MG_pdg_DH1_1[iEvent]),
                                                            MG_pT_DH2[iEvent],
                                                            MG_eta_DH2[iEvent],
                                                            MG_Lxy_tot_DH2[index][iEvent],
                                                            MG_Lz_tot_DH2[index][iEvent],
                                                            abs(MG_pdg_DH2_1[iEvent]),
                                                            selection = "high-ET"))
        MG_eff_highETX.append(sum(MG_queryMapResult))
    MG_queryMapResult = np.array(MG_queryMapResult) # convertion into arrays
    MG_eff_highETX = np.array(MG_eff_highETX) # convertion into arrays
    MG_eff_highETX = MG_eff_highETX/nevent #eff/nbrevent

    
    return MG_eff_highETX

#########################################################################################
# Computation of the efficiency with the map from the data obtained with MG for the low-ET samples (mH <= 400GeV).
#########################################################################################

def eff_map_MG_low(MG_pT_DH1, MG_eta_DH1,MG_Lxy_tot_DH1, MG_Lz_tot_DH1, MG_pdg_DH1_1, MG_pT_DH2, MG_eta_DH2, MG_Lxy_tot_DH2, MG_Lz_tot_DH2, MG_pdg_DH2_1, tauN, nevent, mass_phi, mass_s):

    MG_eff_lowETX = []

    for index in tqdm.tqdm(range(len(tauN))):
        MG_queryMapResult = []
        for iEvent in range(len(MG_pT_DH1)):
            MG_queryMapResult.append(rmN.queryMapFromKinematics(MG_pT_DH1[iEvent],
                                                            MG_eta_DH1[iEvent],
                                                            MG_Lxy_tot_DH1[index][iEvent],
                                                            MG_Lz_tot_DH1[index][iEvent],
                                                            abs(MG_pdg_DH1_1[iEvent]),
                                                            MG_pT_DH2[iEvent],
                                                            MG_eta_DH2[iEvent],
                                                            MG_Lxy_tot_DH2[index][iEvent],
                                                            MG_Lz_tot_DH2[index][iEvent],
                                                            abs(MG_pdg_DH2_1[iEvent]),
                                                            selection = "low-ET"))
        MG_eff_lowETX.append(sum(MG_queryMapResult))
    MG_queryMapResult = np.array(MG_queryMapResult) # convertion into arrays
    MG_eff_lowETX = np.array(MG_eff_lowETX) # convertion into arrays
    MG_eff_lowETX = MG_eff_lowETX/nevent #eff/nbrevent
   
    return MG_eff_lowETX


#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################

#########################################################################################
# Computing HEP data
#########################################################################################

def elem_list(HEP, Branch_HEP, File_HEP_limit, Branch_HEP_limit):

    file_HEP = uproot.open(HEP) # open the file from HEP data for the efficiency
    data_HEP = file_HEP[Branch_HEP] # open the branch

    file_HEP_limit = uproot.open(File_HEP_limit) # open the file from HEP data for the limits
    branch_HEP_limit = file_HEP_limit[Branch_HEP_limit] # open the branch

    return data_HEP, branch_HEP_limit

#########################################################################################
# Plots to compare the results of efficiency obtained with MG, MG+Pythia8 (High-ET).
#########################################################################################

def plt_eff_high(MG_eff_highETX, eff_highETX,tauN, data_HEP,  mass_phi , mass_s, plotDir):

    ################## PLOT EFFICIENCY ##################
    fig, ax = plt.subplots()

    ################## Plot efficiency from MG ##################
    plt.plot(tauN,MG_eff_highETX, 'k--', linewidth=2, label = 'MG')

    ################## Plot efficiency from MG+Pythia8 ##################
    plt.plot(tauN,eff_highETX, 'r', linewidth=2, label = 'MG + Pythia')

    ################## Plot efficiency from HEP data ##################
    plt.plot(data_HEP.values(axis='both')[0],data_HEP.values(axis='both')[1], 'b')

    ################ Uncertainties from HEP ##################
    plt.fill_between(data_HEP.values(axis='both')[0], data_HEP.values(axis='both')[1] +  data_HEP.errors('high')[1] , data_HEP.values(axis='both')[1] - data_HEP.errors('high')[1] , color = 'blue', label = r'ATLAS, with $\pm$ 1 $\sigma$ error bands',alpha=.7)

    ################## Uncertainties from Map ##################
    plt.fill_between(tauN, np.array(eff_highETX) + 0.25* np.array(eff_highETX), np.array(eff_highETX) - 0.25 * np.array(eff_highETX), label='MG+Pythia8, with error bands ', alpha=.7)

    ################## Limits of validity ##################
    ax.hlines(y=(0.25*(max(eff_highETX))), xmin=0, xmax=1e2, linewidth=2, color='g', label = 'Limits of validity' )

    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.05, 0.95, f" $ m_Φ $ = {mass_phi} GeV, $m_S$ = {mass_s} GeV" , transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    x = np.linspace(0,100)
    ax.fill_between(x, 0.25*(max(eff_highETX)), color='black', alpha=.2, hatch="/", edgecolor="black", linewidth=1.0) # adding hatch
    plt.ylim(0) # start at 0

    plt.xscale('log')
    plt.xlabel(r'c$\tau$ [m]', fontsize=20)
    plt.ylabel('Efficiency', fontsize=20 )
    plt.legend(fontsize = 11, loc=1) # set the legend in the upper right corner
    plt.savefig(os.path.join(plotDir,f"Efficiency_comparison_mH{mass_phi}_mS{mass_s}.png"))
    plt.close()


#########################################################################################
# Plots to compared the reasults of efficiency obtained with MG, MG+Pythia8 (Low-ET).
#########################################################################################

def plt_eff_low(MG_eff_lowETX, eff_lowETX,tauN, data_HEP,  mass_phi , mass_s, plotDir):

    ################## PLOT EFFICIENCY ##################
    fig, ax = plt.subplots()

    ################## Plot efficiency from MG ##################
    plt.plot(tauN,MG_eff_lowETX, 'k--',linewidth=2, label = 'MG')

    ################## Plot efficiency from MG+Pythia8 ##################
    plt.plot(tauN,eff_lowETX, 'r', linewidth=2 ,label = 'MG + Pythia')

    ################## Plot efficiency from HEP data ##################
    plt.plot(data_HEP.values(axis='both')[0],data_HEP.values(axis='both')[1], 'b')

    ################ Uncertainties from HEP ##################
    plt.fill_between(data_HEP.values(axis='both')[0], data_HEP.values(axis='both')[1] +  data_HEP.errors('high')[1] , data_HEP.values(axis='both')[1] - data_HEP.errors('high')[1] , color = 'blue', label = r'ATLAS, with $\pm$ 1 $\sigma$ error bands',alpha=.7)

    ################## Uncertainties from Map ##################
    plt.fill_between(tauN, np.array(eff_lowETX) + 0.25* np.array(eff_lowETX), np.array(eff_lowETX) - 0.25*np.array(eff_lowETX), label='MG+Pythia8, with error bands ',alpha=.7)

    ################## Limits of validity ##################
    ax.hlines(y=(0.33*(max(eff_lowETX))), xmin=0, xmax=1e2, linewidth=2, color='g', label = 'Limits of validity' )

    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.05, 0.95, f" $ m_Φ $ = {mass_phi} GeV, $m_S$ = {mass_s} GeV" , transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    x = np.linspace(0,100)
    ax.fill_between(x, 0.33*(max(eff_lowETX)), color='black', alpha=.2, hatch="/", edgecolor="black", linewidth=1.0) # adding hatch
    plt.ylim(0) # start at 0

    plt.xscale('log')
    plt.xlabel(r'c$\tau$ [m]', fontsize=20)
    plt.ylabel('Efficiency', fontsize=20 )
    plt.legend( fontsize = 10, loc=1) # set the legend in the upper right corner
    plt.savefig(os.path.join(plotDir,f"Efficiency_comparison_mH{mass_phi}_mS{mass_s}.png"))
    plt.close()

#########################################################################################
# Plot limits obtained with the map, to compare with those obtain by ATLAS (High-ET).
#########################################################################################

def plt_cross_High(eff_highETX, tauN, mass_phi, mass_s, branch_HEP_limit, factor, plotDir):

    fig, ax = plt.subplots()

    Nsobs = 0.5630 * 26 * factor # nbr of observed events = 26 ( factor )

    Crr_Sec_obs = (Nsobs)/((np.array(eff_highETX)) * 139e3 ) # Luminosity = 139e3 fb**(-1)

    plt.plot(tauN, Crr_Sec_obs, 'r', label ='Map results', linewidth = 2)
    plt.plot(np.array(branch_HEP_limit.values(axis='both')[0]), np.array(branch_HEP_limit.values(axis='both')[1]), 
             'b', label ='Observed', linewidth = 2)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'c$\tau$ [m]')
    plt.ylabel(r'95% CL limit on $\sigma \times B$ [pb]')

    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.05, 0.95, f" $ m_Φ $ = {mass_phi} GeV, $m_S$ = {mass_s} GeV" , transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    plt.legend( fontsize = 10, loc=3)
    plt.savefig(os.path.join(plotDir,f"Cross_section_mH{mass_phi}_mS{mass_s}.png"))
    plt.close()

#########################################################################################
# Plot limits obtained with the map, to those obtain by ATLAS (Low-ET).
#########################################################################################

def plt_cross_Low(eff_lowETX , tauN, mass_phi, mass_s, branch_HEP_limit, factor, plotDir):

    fig, ax = plt.subplots()

    Nsobs = 0.6592 * 26 * factor # nbr of observed events = 26

    Crr_Sec_obs = (Nsobs)/((np.array(eff_lowETX)) * 139e3 )

    plt.plot(tauN, Crr_Sec_obs, 'r', label ='Map results', linewidth = 2)
    plt.plot(tauN, np.array(branch_HEP_limit.values(axis='both')[1]), 'b', label ='Observed', linewidth = 2)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'c$\tau$ [m]')
    plt.ylabel(r'95% CL limit on $\sigma \times B$ [pb]')

    # place a text box in upper left in axes coords
    props = dict(boxstyle='round', facecolor='white', alpha=0.5)
    ax.text(0.05, 0.95, f" $ m_Φ $ = {mass_phi} GeV, $m_S$ = {mass_s} GeV" , transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    plt.legend( fontsize = 10, loc=3)
    plt.savefig(os.path.join(plotDir,f"Cross_section_mH{mass_phi}_mS{mass_s}.png"))
    plt.close()

