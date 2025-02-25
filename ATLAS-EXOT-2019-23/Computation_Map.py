#!/usr/bin/env python3

import sys
import os
sys.path.append('./MG5')

import numpy as np
import mplhep as hep
import lhe_parser as lhe
import hepmc_parser as hepmc
import Computation_Functions as cmfp
import random
import logging

FORMAT = '%(levelname)s: %(message)s at %(asctime)s'
logging.basicConfig(format=FORMAT,datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger()
logger.setLevel(level = logging.DEBUG)    


random.seed(123)
hep.style.use("ATLAS") # Define a style for the plots



fileDict = {
            # (1000,275) : "./Script_mH1000_mS275/Events/run_01/tag_1_pythia8_events.hepmc.gz",
            (1000,275) : "./Script_mH1000_mS275/Events/run_02/highStats_pythia8_events.hepmc.gz",
            # (1000,275) : "./Script_mH1000_mS275_lhapdf/Events/run_01/tag_1_pythia8_events.hepmc.gz",
            # (1000,475) : "./Script_mH1000_mS475/Events/run_01/tag_1_pythia8_events.hepmc.gz",
            }

mg5fileDict = {
            # (1000,275) : './Script_mH1000_mS275/Events/run_01/unweighted_events.lhe.gz',
            (1000,275) : "./Script_mH1000_mS275/Events/run_02/unweighted_events.lhe.gz",
            # (1000,275) : "./Script_mH1000_mS275_lhapdf/Events/run_01/unweighted_events.lhe.gz",
            # (1000,475) : './Script_mH1000_mS475/Events/run_01/unweighted_events.lhe.gz',
            }

fileHEPDict = {(1000,275) : "./ATLAS_data/HEPData-ins2043503-v3-Figure_2e_of_Aux._Mat._1000_275.root",
                (1000,475) : "./ATLAS_data/HEPData-ins2043503-v3-Figure_2e_of_Aux._Mat._1000_475.root",
                } #HEP data files

branchHEPDict = {(1000,275) : "Figure 2e of Aux. Mat. 1000_275/Graph1D_y1;1",
                 (1000,475) : "Figure 2e of Aux. Mat. 1000_475/Graph1D_y1;1",
                } # Branch from HEP data

File_HEP_limit = {(1000,275) : "./ATLAS_data/HEP_Limits/HEPData-ins2043503-v3-Figure_6f_of_Aux._Mat._1000_275.root",
                 (1000,475) : "./ATLAS_data/HEP_Limits/HEPData-ins2043503-v3-Figure_6e_of_Aux._Mat._1000_475.root",
                }

Branch_HEP_limit = {(1000,275) : "Figure 6f of Aux. Mat./Graph1D_y2;1",
                    (1000,475) : "Figure 6e of Aux. Mat./Graph1D_y2;1",
                }


#Constant
c = 3e8# Light velocity in m/s


massPairs = [(1000,275)]
HEP_Lifetime = 95
factor = 1

for mPhi,mS in massPairs:
    tauN=np.geomspace(0.1,1e2,HEP_Lifetime) # New lifetime range
    logger.info(f'Computing results for mPhi = {mPhi} and mS = {mS}')
    if (mPhi,mS) in fileDict:
        #Pythia
        logger.debug('Reading HEPMC')
        events = hepmc.HEPMC_EventFile(fileDict[(mPhi,mS)]) # Open HEPMC file
        nevts_pythia = len(events)
        logger.debug('Parsing HEPMC')
        px_TOT, py_TOT, pz_TOT, E_TOT, mass_TOT,pdg_TOT = cmfp.parsing_hepmc(events) # Parsing the HEPMC file
        px_tot, py_tot, pz_tot, E_tot, mass_tot, pdg_tot = cmfp.conversion_one_list(px_TOT, py_TOT, pz_TOT, E_TOT, mass_TOT, pdg_TOT) # Obtaining data in one list
        px_DH1, px_DH2, py_DH1, py_DH2, pz_DH1, pz_DH2, pdg_tot_DH1, pdg_tot_DH2, E_DH1, E_DH2, mass_DH1, mass_DH2 = cmfp.recover(px_tot, py_tot, pz_tot, E_tot, mass_tot, pdg_tot) # Separate data from DH1 and DH2
        beta_DH1, gamma_DH1, pT_DH1, eta_DH1 = cmfp.kinematics_DH(px_DH1, py_DH1, pz_DH1, E_DH1) # Computing kinematics for DH1
        beta_DH2, gamma_DH2, pT_DH2, eta_DH2 = cmfp.kinematics_DH(px_DH2, py_DH2, pz_DH2, E_DH2) # Computing kinematics for DH2
        Lxy_tot_DH1, Lz_tot_DH1 = cmfp.decaylenghtDH(px_DH1, py_DH1, pz_DH1, E_DH1, gamma_DH1, tauN) # Computing the decay lenght for DH1
        Lxy_tot_DH2, Lz_tot_DH2 = cmfp.decaylenghtDH(px_DH2, py_DH2, pz_DH2, E_DH2, gamma_DH2, tauN) # Computing the decay lenght for DH2

    if (mPhi,mS) in mg5fileDict:
        #MG
        logger.debug('Reading LHE')
        MG_events = lhe.EventFile(mg5fileDict[(mPhi,mS)]) # Open LHE file
        nevts_mg5 = len(MG_events)
        px, py, pz, pdg, E, MASS = cmfp.parsing_LHE(MG_events) #Parsing the LHE file
        MG_px_DH1, MG_py_DH1,MG_pz_DH1,MG_E_DH1,MG_mass_DH1,MG_pdg_DH1_1 = cmfp.recover_MG_DH1(px, py, pz, E, MASS, pdg) # Separate data from DH1 and DH2
        MG_pT_DH1,MG_eta_DH1, MG_gamma_DH1 = cmfp.kinematics_MG_DH1(MG_px_DH1,MG_py_DH1,MG_pz_DH1,MG_E_DH1) # Computing kinematics for DH1
        MG_px_DH2, MG_py_DH2,MG_pz_DH2,MG_E_DH2,MG_mass_DH2,MG_pdg_DH2_1 = cmfp.recover_MG_DH2(px, py, pz, E, MASS, pdg) # Separate data from DH1 and DH2
        MG_pT_DH2,MG_eta_DH2, MG_gamma_DH2 = cmfp.kinemamtics_MG_DH2(MG_px_DH2,MG_py_DH2,MG_pz_DH2,MG_E_DH2) # Computing kinematics for DH2
        MG_Lxy_tot_DH1, MG_Lz_tot_DH1 = cmfp.decaylenght_MG_DH1(MG_px_DH1, MG_py_DH1, MG_pz_DH1, E_DH1, MG_gamma_DH1, tauN) # Computing decay lenght for DH1
        MG_Lxy_tot_DH2, MG_Lz_tot_DH2 = cmfp.decaylenght_MG_DH2(MG_px_DH2, MG_py_DH2, MG_pz_DH2, E_DH2, MG_gamma_DH2, tauN) # Computing decay lenght for DH2

    #HEP data
    data_HEP, branch_HEP_limit = cmfp.elem_list(fileHEPDict[(mPhi,mS)], branchHEPDict[(mPhi,mS)], File_HEP_limit[(mPhi,mS)], Branch_HEP_limit[(mPhi,mS)]) # Recover public data from ATLAS to compare the results

############################################################################################################################################################################
############################################################################################################################################################################
############################################################################################################################################################################

###########################################################Computing the efficiencies and ploting the results###########################################################

    logger.info(f'Computing efficiencies')
    if mPhi >= 400: # Condition if the sample is 'High-ET' or ' Low-ET'
        MG_eff_highETX = cmfp.eff_map_MG_high(MG_pT_DH1, MG_eta_DH1,MG_Lxy_tot_DH1, MG_Lz_tot_DH1, MG_pdg_DH1_1, MG_pT_DH2, MG_eta_DH2, MG_Lxy_tot_DH2, MG_Lz_tot_DH2, MG_pdg_DH2_1, tauN, nevts_mg5,  mPhi, mS) # Compute the efficiency from MG
        MG_Data_Eff_High = np.column_stack(MG_eff_highETX)
        np.savetxt(os.path.join(os.path.dirname(mg5fileDict[(mPhi,mS)]),f'Efficiencies_Text_{mPhi}_{mS}.txt'), MG_Data_Eff_High)

        eff_highETX = cmfp.eff_map_High(pT_DH1, eta_DH1, Lxy_tot_DH1, Lz_tot_DH1, abs(pdg_tot_DH1), pT_DH2, eta_DH2, Lxy_tot_DH2, Lz_tot_DH2, abs(pdg_tot_DH2), tauN, nevts_pythia,  mPhi, mS) # Compute the efficiency from Pythia
        Data_Eff_High = np.column_stack(eff_highETX)
        plotDir = os.path.dirname(fileDict[(mPhi,mS)])
        np.savetxt(os.path.join(plotDir,f'Efficiencies_Text_{mPhi}_{mS}.txt'), Data_Eff_High)

        cmfp.plt_eff_high(MG_eff_highETX, eff_highETX, tauN, data_HEP, mPhi, mS, plotDir ) # Ploting and saving a comparison of all the results of efficiencies
        cmfp.plt_cross_High(eff_highETX, tauN, mPhi, mS, branch_HEP_limit, factor, plotDir)# Ploting and saving a comparison of the limits obtained with the map and by ATLAS.

    else:
        MG_eff_lowETX = cmfp.eff_map_MG_low(MG_pT_DH1, MG_eta_DH1,MG_Lxy_tot_DH1, MG_Lz_tot_DH1, MG_pdg_DH1_1, MG_pT_DH2, MG_eta_DH2, MG_Lxy_tot_DH2, MG_Lz_tot_DH2, MG_pdg_DH2_1, tauN, nevts_mg5, mPhi, mS)
        MG_Data_Eff_Low = np.column_stack(MG_eff_highETX)
        np.savetxt(os.path.join(os.path.dirname(mg5fileDict[(mPhi,mS)]),f'Efficiencies_Text_{mPhi}_{mS}.txt'), MG_Data_Eff_Low)
        eff_lowETX = cmfp.eff_map_Low(pT_DH1, eta_DH1, Lxy_tot_DH1, Lz_tot_DH1, abs(pdg_tot_DH1), pT_DH2, eta_DH2, Lxy_tot_DH2, Lz_tot_DH2, abs(pdg_tot_DH2), tauN, nevts_pythia, mPhi, mS)
        Data_Eff_Low = np.column_stack(eff_lowETX)
        plotDir = os.path.dirname(fileDict[(mPhi,mS)])        
        np.savetxt(os.path.join(plotDir,f'Efficiencies_Text_{mPhi}_{mS}.txt'), Data_Eff_Low)

        cmfp.plt_eff_low(MG_eff_lowETX, eff_lowETX, tauN, data_HEP, mPhi, mS, plotDir)
        cmfp.plt_cross_Low(eff_lowETX, tauN, mPhi, mS, branch_HEP_limit, factor, plotDir)
    
