/*
 * DelphesPythia8Matched
 *
 * Runs Pythia8 (with MLM jet matching a la MadGraph, i.e. JetMatchingMadgraph)
 * on an LHE file and passes the showered/hadronized events directly to Delphes,
 * without writing an intermediate HepMC file.
 *
 * Differences with respect to the standard DelphesPythia8 reader:
 *  - Check:abortIfVeto is forced on, so that every call to Pythia::next()
 *    consumes exactly one LHE event. This allows us to
 *      (a) count matching vetoes separately from genuine errors (the vanilla
 *          reader counts nothing and silently loops over vetoed events), and
 *      (b) process a well-defined slice of the LHE file
 *          [Beams:nSkipLHEFatInit, Beams:nSkipLHEFatInit + Main:numberOfEvents),
 *          so several instances can run in parallel on the same LHE file.
 *  - Main:numberOfEvents is interpreted as the number of *LHE events* to read
 *    (-1 = until end of file), not the number of accepted events.
 *  - The DelphesLHEFReader side-channel is not used (it cannot read .lhe.gz
 *    and desynchronizes from Pythia when events are vetoed).
 *  - A summary (counts and cross sections) is written to an optional text file.
 *
 * Usage: DelphesPythia8Matched config_file pythia_card output_file [summary_file]
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include <signal.h>

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/JetMatching.h"

#include "TApplication.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TParticlePDG.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "modules/Delphes.h"

#include "ExRootAnalysis/ExRootConfReader.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"

using namespace std;

//---------------------------------------------------------------------------
// JetMatchingMadgraph hook that records whether the current event was vetoed.

class JetMatchingMadgraphTracked : public Pythia8::JetMatchingMadgraph
{
public:
  bool vetoed = false;

  bool doVetoProcessLevel(Pythia8::Event &process) override
  {
    bool veto = Pythia8::JetMatchingMadgraph::doVetoProcessLevel(process);
    vetoed = vetoed || veto;
    return veto;
  }

  bool doVetoPartonLevelEarly(const Pythia8::Event &event) override
  {
    bool veto = Pythia8::JetMatchingMadgraph::doVetoPartonLevelEarly(event);
    vetoed = vetoed || veto;
    return veto;
  }
};

//---------------------------------------------------------------------------
// Same conversion as in readers/DelphesPythia8.cpp

void ConvertInput(Long64_t eventNumber, Pythia8::Pythia *pythia,
  ExRootTreeBranch *branch, DelphesFactory *factory,
  TObjArray *allParticleOutputArray, TObjArray *stableParticleOutputArray, TObjArray *partonOutputArray,
  TStopwatch *readStopWatch, TStopwatch *procStopWatch)
{
  HepMCEvent *element = static_cast<HepMCEvent *>(branch->NewEntry());

  element->Number = eventNumber;
  element->ProcessID = pythia->info.code();
  element->MPI = 1;
  element->Weight = pythia->info.weight();
  element->Scale = pythia->info.QRen();
  element->AlphaQED = pythia->info.alphaEM();
  element->AlphaQCD = pythia->info.alphaS();
  element->ID1 = pythia->info.id1();
  element->ID2 = pythia->info.id2();
  element->X1 = pythia->info.x1();
  element->X2 = pythia->info.x2();
  element->ScalePDF = pythia->info.QFac();
  element->PDF1 = pythia->info.pdf1();
  element->PDF2 = pythia->info.pdf2();
  element->ReadTime = readStopWatch->RealTime();
  element->ProcTime = procStopWatch->RealTime();

  TDatabasePDG *pdg = TDatabasePDG::Instance();

  for(int i = 1; i < pythia->event.size(); ++i)
  {
    Pythia8::Particle &particle = pythia->event[i];

    Candidate *candidate = factory->NewCandidate();

    candidate->PID = particle.id();
    Int_t pdgCode = TMath::Abs(candidate->PID);
    Int_t status = particle.statusHepMC();
    candidate->Status = status;

    candidate->M1 = particle.mother1() - 1;
    candidate->M2 = particle.mother2() - 1;
    candidate->D1 = particle.daughter1() - 1;
    candidate->D2 = particle.daughter2() - 1;

    TParticlePDG *pdgParticle = pdg->GetParticle(particle.id());
    candidate->Charge = pdgParticle ? Int_t(pdgParticle->Charge() / 3.0) : -999;
    candidate->Mass = particle.m();

    candidate->Momentum.SetPxPyPzE(particle.px(), particle.py(), particle.pz(), particle.e());
    candidate->Position.SetXYZT(particle.xProd(), particle.yProd(), particle.zProd(), particle.tProd());

    allParticleOutputArray->Add(candidate);

    if(!pdgParticle) continue;

    if(status == 1)
      stableParticleOutputArray->Add(candidate);
    else if(pdgCode <= 5 || pdgCode == 21 || pdgCode == 15)
      partonOutputArray->Add(candidate);
  }
}

//---------------------------------------------------------------------------

static bool interrupted = false;

void SignalHandler(int sig)
{
  interrupted = true;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  char appName[] = "DelphesPythia8Matched";

  if(argc != 4 && argc != 5)
  {
    cout << " Usage: " << appName << " config_file pythia_card output_file [summary_file]" << endl;
    cout << " config_file  - Delphes configuration file in Tcl format," << endl;
    cout << " pythia_card  - Pythia8 configuration file (Main:numberOfEvents = number of LHE events to read)," << endl;
    cout << " output_file  - output file in ROOT format," << endl;
    cout << " summary_file - (optional) text file with event counts and cross sections." << endl;
    return 1;
  }

  signal(SIGINT, SignalHandler);
  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);

  TFile *outputFile = 0;
  ExRootTreeWriter *treeWriter = 0;
  ExRootConfReader *confReader = 0;
  Delphes *modularDelphes = 0;
  Pythia8::Pythia *pythia = 0;

  try
  {
    stringstream message;

    // ---- Pythia setup ----
    pythia = new Pythia8::Pythia;

    if(!pythia->readFile(argv[2]))
    {
      message << "can't read Pythia8 configuration file " << argv[2];
      throw runtime_error(message.str());
    }
    // Each next() must correspond to exactly one LHE event
    pythia->readString("Check:abortIfVeto = on");

    if(pythia->mode("Beams:frameType") != 4 && pythia->mode("Beams:frameType") != 5)
      throw runtime_error("this reader requires LHE input (Beams:frameType = 4 or 5)");

    const bool doMatching = pythia->flag("JetMatching:merge");
    if(doMatching && pythia->mode("JetMatching:scheme") != 1)
      throw runtime_error("only JetMatching:scheme = 1 (MadGraph MLM) is supported");

    shared_ptr<JetMatchingMadgraphTracked> matchingHook;
    if(doMatching)
    {
      matchingHook = make_shared<JetMatchingMadgraphTracked>();
      pythia->setUserHooksPtr(matchingHook);
    }

    const Long64_t maxLHE = pythia->mode("Main:numberOfEvents");
    const Long64_t nSkip = pythia->mode("Beams:nSkipLHEFatInit");
    const Long64_t timesAllowErrors = pythia->mode("Main:timesAllowErrors");

    // ---- Delphes setup ----
    outputFile = TFile::Open(argv[3], "CREATE");
    if(outputFile == NULL)
    {
      message << "can't create output file " << argv[3];
      throw runtime_error(message.str());
    }

    treeWriter = new ExRootTreeWriter(outputFile, "Delphes");
    ExRootTreeBranch *branchEvent = treeWriter->NewBranch("Event", HepMCEvent::Class());
    ExRootTreeBranch *branchWeight = treeWriter->NewBranch("Weight", Weight::Class());

    confReader = new ExRootConfReader;
    confReader->ReadFile(argv[1]);

    modularDelphes = new Delphes("Delphes");
    modularDelphes->SetConfReader(confReader);
    modularDelphes->SetTreeWriter(treeWriter);

    DelphesFactory *factory = modularDelphes->GetFactory();
    TObjArray *allParticleOutputArray = modularDelphes->ExportArray("allParticles");
    TObjArray *stableParticleOutputArray = modularDelphes->ExportArray("stableParticles");
    TObjArray *partonOutputArray = modularDelphes->ExportArray("partons");

    modularDelphes->InitTask();

    if(!pythia->init())
      throw runtime_error("Pythia initialization failed");

    // ---- Event loop ----
    Long64_t nLHE = 0, nAccepted = 0, nVetoed = 0, nErrors = 0;
    bool reachedEOF = false, tooManyErrors = false;
    TStopwatch readStopWatch, procStopWatch;

    treeWriter->Clear();
    modularDelphes->Clear();
    readStopWatch.Start();

    while(!interrupted && (maxLHE < 0 || nLHE < maxLHE))
    {
      if(matchingHook) matchingHook->vetoed = false;

      bool ok = pythia->next();

      if(!ok && pythia->info.atEndOfFile())
      {
        reachedEOF = true;
        break;
      }
      ++nLHE;

      if(!ok)
      {
        if(matchingHook && matchingHook->vetoed)
          ++nVetoed;
        else if(++nErrors > timesAllowErrors)
        {
          cerr << "** ERROR: too many Pythia errors (" << nErrors << "), aborting." << endl;
          tooManyErrors = true;
          break;
        }
        continue;
      }

      readStopWatch.Stop();
      procStopWatch.Start();
      ConvertInput(nSkip + nLHE - 1, pythia, branchEvent, factory,
        allParticleOutputArray, stableParticleOutputArray, partonOutputArray,
        &readStopWatch, &procStopWatch);
      modularDelphes->ProcessTask();
      procStopWatch.Stop();

      for(size_t iWeight = 0; iWeight < pythia->info.weightValueVector().size(); ++iWeight)
      {
        Weight *weight = static_cast<Weight *>(branchWeight->NewEntry());
        weight->Weight = pythia->info.weightValueVector()[iWeight];
      }

      treeWriter->Fill();
      treeWriter->Clear();
      modularDelphes->Clear();
      ++nAccepted;

      readStopWatch.Start();
      if(nLHE % 1000 == 0)
        cout << "** Processed " << nLHE << " LHE events, accepted " << nAccepted << endl;
    }

    pythia->stat();

    modularDelphes->FinishTask();
    treeWriter->Write();

    // ---- Cross sections ----
    // sigma_LHE: sum of XSECUP from the LHE <init> block (pb).
    // Matched cross section assuming unweighted (equal-weight) LHE events:
    //   sigma_match = sigma_LHE * nAccepted / nLHE
    double sigmaLHE = 0.0;
    for(int i = 0; i < pythia->info.nProcessesLHEF(); ++i)
      sigmaLHE += pythia->info.sigmaLHEF(i);
    double eff = nLHE > 0 ? double(nAccepted) / double(nLHE) : 0.0;
    double effErr = nLHE > 0 ? sqrt(eff * (1.0 - eff) / double(nLHE)) : 0.0;

    ostringstream summary;
    summary.precision(8);
    summary << "nSkip = " << nSkip << "\n"
            << "nLHE = " << nLHE << "\n"
            << "nAccepted = " << nAccepted << "\n"
            << "nVetoed = " << nVetoed << "\n"
            << "nErrors = " << nErrors << "\n"
            << "reachedEOF = " << reachedEOF << "\n"
            << "aborted = " << (tooManyErrors || interrupted) << "\n"
            << "qCut = " << pythia->parm("JetMatching:qCut") << "\n"
            << "sigmaLHE_pb = " << sigmaLHE << "\n"
            << "sigmaMatched_pb = " << sigmaLHE * eff << "\n"
            << "sigmaMatchedErr_pb = " << sigmaLHE * effErr << "\n"
            << "sigmaGenPythia_pb = " << pythia->info.sigmaGen() * 1e9 << "\n"
            << "sigmaErrPythia_pb = " << pythia->info.sigmaErr() * 1e9 << "\n";

    cout << "** Summary:\n" << summary.str();
    if(argc == 5)
    {
      ofstream sumFile(argv[4]);
      sumFile << summary.str();
    }

    cout << "** Exiting..." << endl;

    delete pythia;
    delete modularDelphes;
    delete confReader;
    delete treeWriter;
    delete outputFile;

    return (tooManyErrors || interrupted) ? 2 : 0;
  }
  catch(runtime_error &e)
  {
    if(treeWriter) delete treeWriter;
    if(outputFile) delete outputFile;
    cerr << "** ERROR: " << e.what() << endl;
    return 1;
  }
}
