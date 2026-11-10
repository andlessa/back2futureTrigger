#For the pheno studies we use the following pdf
#FIXME The following lines should be removed for official ATLAS production
import MadGraphControl.MadGraphUtils
MadGraphControl.MadGraphUtils.MADGRAPH_PDFSETTING={
    'central_pdf':315000, # the lhapf id of the central pdf, see https://lhapdf.hepforge.org/pdfsets
    'pdf_variations':[315000], # list of pdfs ids for which all variations (error sets) will be included as weights
    'alternative_pdfs':None, # list of pdfs ids for which only the central set will be included as weights
    'scale_variations':[0.5,1,2], # variations of muR and muF wrt the central scale, all combinations of muF and muR will be evaluated
} 
#FIXME instead, we should un-comment the following line
#import MadGraphControl.MadGraph_NNPDF30NLO_Base_Fragment #Importing ATLAS recommended PDF settings

from MadGraphControl.MadGraphUtils import *
from MadGraphControl.MadGraphUtilsHelpers import get_physics_short
tokens = get_physics_short().split('_')
# Helper for resetting process number
from MadGraphControl.MadGraphUtils import check_reset_proc_number
import os

#---------------------------------------------------------------------------
# Process-specific parameter settings in MadGraph
#---------------------------------------------------------------------------

mS = 450
m2 = 197.97475910501439
m1 = 192.97475910501439
m0 = 132.97475910501439
wS = 'auto'
w2 = 'auto'
w1 = 'auto'
w0 = 'auto'
Chxx10 = 0.1
lambdaUV = 9069683.515309807
cTau = 20
sinAlpha = 0.2

#----------------------------------------------------------------------------
# Process definition
#----------------------------------------------------------------------------

# Make sure PDG IDs are known to TestHepMC:
pdgfile = open("pdgid_extras.txt", "w+")
pdgfile.write("""
4000022
4000023
4000024
-4000022
-4000023
-4000024
55
""")
pdgfile.close()

process="""
import model /eos/user/j/jagargan/cambridge/b2f/splitFrameCTP/mc_gen/monoJet/edmgeneratorstudies/MadGraphModels/DM/DipoleDM2_minimalH_UFO

define p = g u c d s u~ c~ d~ s~
define j = g u c d s u~ c~ d~ s~
define l+ = e+ mu+
define l- = e- mu-
define vl = ve vm vt
define vl~ = ve~ vm~ vt~
define c1 = ~chi1 ~chi1bar
define c2 = ~chi2 ~chi2bar
define c0 = ~chi0 ~chi0bar

generate p p > c1 c2 j / h, (c1 > b b~ c0), (c2 > b b~ c0)

output -f
"""

#----------------------------------------------------------------------------
# Beam energy
#----------------------------------------------------------------------------
beamEnergy=-999
if hasattr(runArgs,'ecmEnergy'):
  beamEnergy = runArgs.ecmEnergy / 2.
else:
  raise RuntimeError('No center of mass energy found.')

#---------------------------------------------------------------------------
# Number of events to generate
#---------------------------------------------------------------------------

nevents=evgenConfig.nEventsPerJob
multiplier=3
nevents = runArgs.maxEvents*multiplier if runArgs.maxEvents>0 else nevents*multiplier
nevents = int(nevents)

skip_events=0
if hasattr(runArgs,'skipEvents'): skip_events=runArgs.skipEvents

#---------------------------------------------------------------------------
# MG5 Run Card
#---------------------------------------------------------------------------

extras_runcard = {
    'nevents'           : int(nevents),
    'lhe_version'       : '3.0',
    'ickkw'             : 0, #MLM
    'maxjetflavor'      : '4',
    'xqcut'		: '0',
    'auto_ptj_mjj'      : 'False',
    'ptj'		: '100',
    'etaj'		: '5.0',
    'time_of_flight'    : '0.0', # threshold (in mm) below which the invariant livetime is not written (-1 means not written), default value is -1
    'sde_strategy'      : '2', #new Phase-Space integration optimization strategy
}

# Build run_card
process_dir = new_process(process)
modify_run_card(process_dir=process_dir,runArgs=runArgs,settings=extras_runcard)

#---------------------------------------------------------------------------
# MG5 parameter Card
#---------------------------------------------------------------------------

# Build param_card.dat

extras_paramcard = {
        'DECAY'  :{
        'WS': wS,
        'w0': w1,
        'w1': w0,
	'w2': w2,
    },
    'MASS'   :{
       'MSd' : mS,
       'M2'  : m2,
       'M1' : m1,
       'M0' : m0,
    },
    'npinputs' :{
        'sina'  : sinAlpha,     #mixing S and h (default 0.2)
        'ychi20' : 0.1,            #only higgs channel
        'ychi21' : 1.5,            #only higgs channel
        'Chxx10' : Chxx10,      #EFT parameter chi1->chi0+h* (default 0.1) #linked to lifetime
        'LambdaUV' : lambdaUV ,
    }
}

modify_param_card(process_dir=process_dir,params=extras_paramcard)

#---------------------------------------------------------------------------
# Generate the events
#---------------------------------------------------------------------------

generate(process_dir=process_dir,runArgs=runArgs)
add_lifetimes(process_dir=process_dir,threshold=0) #see https://gitlab.cern.ch/atlas/athena/-/tree/main/Generators/MadGraphControl#adding-lifetimes
arrange_output(process_dir=process_dir,runArgs=runArgs,lhe_version=3,saveProcDir=True)

#---------------------------------------------------------------------------
# Metadata
#---------------------------------------------------------------------------

evgenConfig.process = "p p > ~chi1 ~chi1bar / h"
evgenConfig.description = "Slow LLP production test"
evgenConfig.keywords = ["Higgs"]
evgenConfig.tune = "A14 NNPDF23LO"

#---------------------------------------------------------------------------
# Shower
#---------------------------------------------------------------------------

#Unlike Madgraph, Pythia cannot run on multi-core mode, go back to single core before executing showering
check_reset_proc_number(opts)

#Pythia Settings
#PYTHIA8_nJetMax=1
#PYTHIA8_qCut=1.5*mS/10 #Pythia's Qcut should be somewhat larger than xqcut; a typical factor would be Qcut = 1.5 * xqcut where xqcut ~ 1/10 of the hard scale
#PYTHIA8_coneRadius=1.00000
#PYTHIA8_nQmatch=4

#include tunes
include("Pythia8_i/Pythia8_A14_NNPDF23LO_EvtGen_Common.py")
include("Pythia8_i/Pythia8_MadGraph.py")
#The following settings are following the Pythia8_MLM_ClassicalMatch.py include
#include("Pythia8_i/Pythia8_MLM_ClassicalMatch.py")

genSeq.EvtInclusiveDecay.whiteList += [4000024,-4000024,4000023,-4000023,4000022,-4000022]
testSeq.TestHepMC.MaxTransVtxDisp = 100000000 #in mm
testSeq.TestHepMC.MaxVtxDisp = 100000000 #in mm
testSeq.TestHepMC.MaxNonG4Energy = 100000000 #in MeV

# decrease relative precision
genSeq.Pythia8.Commands += ["Check:epTolErr               = 1e-2"]

#if "UserHooks" in genSeq.Pythia8.__slots__.keys():
#    genSeq.Pythia8.UserHooks += ['JetMatchingMadgraph']
#else:
#    genSeq.Pythia8.UserHook = 'JetMatchingMadgraph'

genSeq.Pythia8.Commands += ['Tune:pp = 21', #FIXME maybe change/ remove this setting for official ATLAS production
                            'Init:showChangedSettings = on',
                            'Init:showChangedParticleData = off',
                            'Init:showProcesses = on',
                            'Init:showMultipartonInteractions = off',
                            'Next:numberCount = 15000',
                            'Next:numberShowInfo = 1',
                            'Next:numberShowProcess = 1',
                            'Next:numberShowEvent = 1',
                            'Next:numberShowLHA  = 0',
                            'Stat:showProcessLevel = on',
                            'Stat:showErrors = on',
                            'Main:timesAllowErrors = 100',
                            #FIXME we might want to remove the following settings for an official ATLAS production
                            'StringZ:useOldAExtra = off',
                            'StandardModel:sin2thetaW = 0.23120',
                            'StandardModel:sin2thetaWbar = 0.23150',
                            'ParticleDecays:limitTau0 = off'
                            ]

if "computeEfficiency" in genSeq.Pythia8.__slots__.keys():
    genSeq.Pythia8.computeEfficiency = True


#---------------------------------------------------------------------------
# MC Generation Filter
#---------------------------------------------------------------------------

#Add filters below when needed 
#Remember to increase the safety factor accordingly