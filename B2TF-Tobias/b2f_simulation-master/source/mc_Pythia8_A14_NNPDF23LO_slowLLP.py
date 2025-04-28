
evgenConfig.description = "Slow LLP production test"
evgenConfig.keywords = ["Higgs"]
evgenConfig.specialConfig = "DECAYS=True;BTFTChi1Mass=248700;BTFTChi1Lifetime=3000;BTFTChi0Mass=173700"
#evgenConfig.specialConfig = "preInclude=DelayedParticles.DelayedParticlesConfig.DelayedParticlesWritePreInclude; postInclude=DelayedParticles.DelayedParticlesConfig.DelayedParticlesWritePostInclude"
#evgenConfig.specialConfig = "preInclude=DelayedParticles.DelayedParticlesConfig.DelayedParticlesWritePreInclude; postInclude=DelayedParticles.DelayedParticlesConfig.DelayedParticlesWritePostInclude; CosmicFilterVolumeNames=Calo; ParticleID=1000023"
#evgenConfig.specialConfig = "preInclude=DelayedParticles.DelayedParticlesConfig.DelayedParticlesReadPreInclude"

#include("MC15JobOptions/Pythia8_A14_NNPDF23LO_EvtGen_Common.py") what's in the gitlab, but isn't found?
include("Pythia8_i/Pythia8_A14_NNPDF23LO_EvtGen_Common.py")

genSeq.Pythia8.Commands += ['Higgs:useBSM = on']
genSeq.Pythia8.Commands += ['HiggsBSM:gg2H2 = on']
genSeq.Pythia8.Commands += ['35:m0 = 500']
genSeq.Pythia8.Commands += ['35:mwidth = 0.1448849']
genSeq.Pythia8.Commands += ['1000023:m0 = 248.746859276655']
genSeq.Pythia8.Commands += ['1000022:m0 = 173.746859276655']
genSeq.Pythia8.Commands += ['1000022:mayDecay = off']
genSeq.Pythia8.Commands += ['1000023:mayDecay = on']
genSeq.Pythia8.Commands += ['1000023:tauCalc = false']
genSeq.Pythia8.Commands += ['35:oneChannel 1 0.00000000001 100 21 21']
genSeq.Pythia8.Commands += ['35:addChannel 1 0.99999999999 100 1000023 1000023']
genSeq.Pythia8.Commands += ['1000023:oneChannel 1 1.000 100 1000022 22']
genSeq.Pythia8.Commands += ['1000023:tau0 = 3000']

testSeq.TestHepMC.MaxTransVtxDisp = 100000000 #in mm
testSeq.TestHepMC.MaxVtxDisp = 100000000 #in mm
testSeq.TestHepMC.MaxNonG4Energy = 100000000 #in MeV
