#include "MCVal/TruthUtils.h"
#ifndef __DC14__
#endif
#include <vector>

#include "xAODBase/IParticleHelpers.h"
#include "TruthUtils/HepMCHelpers.h"
// for documentaiton, see the header file

namespace Truth {
  
  /// print details about the truth particle to the screen
  void printTruthPtcl(const xAOD::TruthParticle *ptcl, TString comment, 
		      int childDepth, int parentDepth, int currentDepth) {
    // indentation in case we print decay chain. Three spaces per level
    TString indent(Form(Form("%%%ds",3*currentDepth),""));
    if (ptcl==NULL) { printf("%sNULL\n",indent.Data()); return; }
    printf("%sTruth part. ID:%5d, status: %2d, pT: %5f, eta: %5f, phi: %5f, mass: %5f %s\n",
           indent.Data(),ptcl->pdgId(),ptcl->status(),ptcl->pt(),ptcl->eta(),ptcl->phi(),ptcl->m(),comment.Data());
    if (childDepth>0||parentDepth>0) {
      int npar=ptcl->nParents(), nchild=ptcl->nChildren();
      printf("%s-> %d parent and %d children\n",indent.Data(),npar,nchild);
      if (parentDepth>0) 
        for (int ip=0;ip<npar;++ip) printTruthPtcl(ptcl->parent(ip),Form("parent %d of ",ip+1)+comment,
                                                   childDepth-1,parentDepth-1,currentDepth+1);
      if (childDepth>0)
        for (int ic=0;ic<nchild;++ic) printTruthPtcl(ptcl->child(ic),Form("child %d of ",ic+1)+comment,
                                                     childDepth-1,parentDepth-1,currentDepth+1);
    }
  }
  
  bool isStable(const xAOD::TruthParticle *ptcl) {
    return ptcl->status() == 1 && ptcl->barcode() < 200000;
  }
  
  bool isDalitz(const xAOD::TruthParticleContainer *truthPtcls) {
    for (auto ptcl:*truthPtcls) // if H -> y*
      if ( abs(ptcl->pdgId())==25 && ptcl->status()==62 && ptcl->nChildren()>=2 &&
          ((ptcl->child(0)&&ptcl->child(0)->pdgId()==22&&ptcl->child(0)->status()!=1)||
           (ptcl->child(1)&&ptcl->child(1)->pdgId()==22&&ptcl->child(1)->status()!=1)) )
        return true;
    return false;
  }
  
  // Return true if not from hadron
  bool notFromHadron(const xAOD::TruthParticle *ptcl) {
    int ID = ptcl->pdgId();
    
    // if the particle is a hadron, return false
    if (MC::isHadron(ID)) return false;
    
    // if there are no parents, not from hadron
    if (ptcl->nParents()==0) return true;
    
    const xAOD::TruthParticle *parent = ptcl->parent(0);
    int parentID = parent->pdgId();
    if (MC::isHadron(parentID)) return false; // from hadron!
    if (parentID==15||parentID==ID) return notFromHadron(parent);
    
    // if we get here, all is good
    return true;
  }
 
  TLorentzVector getTLV(const xAOD::TruthParticle *ptcl) {
    TLorentzVector tempTLV;
    tempTLV.SetPtEtaPhiM(ptcl->pt(),ptcl->eta(),ptcl->phi(),ptcl->m());
    return tempTLV;
  }

  // adds up 4-vectors of all stable particles:
  //  should always give E=m=sqrt(s) \vec{p}=\vec{0} !
  // unless partciels are missing from the file
  TLorentzVector getStableParticle4VectorSum(const xAOD::TruthParticleContainer *truthPtcls) {
    TLorentzVector sum;
    for ( const xAOD::TruthParticle *ptcl : *truthPtcls)
      if (isStable(ptcl)) sum += ptcl->p4();
    return sum;
  }

  /*
  double getTruthIsolation(const xAOD::TruthParticle *ptcl,
                           const xAOD::TruthParticleContainer *truthPtcls,
                           double dr,
                           bool chargeOnly,
                           std::vector<int> ignorePdgIds)
  {
    // Pointer to be used in this function
    const xAOD::TruthParticle *_ptcl = ptcl;

    // Check if this points back to an original particle
    static SG::AuxElement::ConstAccessor<ElementLink<xAOD::IParticleContainer> > acc("originalObjectLink");
    if (acc.isAvailable(*ptcl)) {
      _ptcl = dynamic_cast<const xAOD::TruthParticle*>(xAOD::getOriginalObject(*ptcl));
    }

    // Calculate isolation
    TLorentzVector iso(0,0,0,0);
    for (auto p: *truthPtcls) {
      // Don't count the particles own energy
      if (p->barcode() == _ptcl->barcode())
        continue;

      // Only consider stable particles
      if (not isStable(p))
        continue;

      // Must be withing the dR cone
      if (HWW::DR(p, _ptcl) >= dr)
        continue;

      // Only include charge particles?
      if (chargeOnly && p->threeCharge() == 0)
        continue;

      // Don't consider muons or neutrinos
      if (std::find(ignorePdgIds.begin(), ignorePdgIds.end(), abs(p->pdgId())) != ignorePdgIds.end())
        continue;

      iso += p->p4();
    }

    if (iso.Px() == 0 && iso.Py() == 0)
      return 0.0;

    return iso.Et();
  }
  */

  const xAOD::TruthParticle* getBosonTopAncestor(const xAOD::TruthParticle *ptcl) {
    if(ptcl==nullptr) return ptcl;
    const xAOD::TruthParticle *parent = ptcl->parent();
    if(parent==nullptr) return ptcl;
    if(ptcl->pdgId()==parent->pdgId()) return getBosonTopAncestor(parent);
    if(MC::isW(parent->pdgId())) return parent;
    if(MC::isZ(parent->pdgId())) return parent;
    if(MC::isHiggs(parent->pdgId())) return parent;
    if(MC::isTop(parent->pdgId())) return parent;
    return getBosonTopAncestor(parent);
  }

  const xAOD::TruthParticle* getFirstAncestor(const xAOD::TruthParticle *ptcl) {
    if(ptcl==nullptr) return ptcl;
    const xAOD::TruthParticle *parent = ptcl->parent();
    if(parent==nullptr) return ptcl;
    if(parent->pdgId() != ptcl->pdgId()) return parent;
    return getFirstAncestor(parent);
  }

  bool isFinalHiggs(const xAOD::TruthParticle *part)
  {
    if (!MC::isHiggs(part->pdgId())) return false;
    if (part->child() == nullptr) return false;
    for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
      if (MC::isHiggs(part->child(iChild)->pdgId())) return false;
    }
    return true;
  }

  bool isHiggsToWW(const xAOD::TruthParticle *part)
  {
    if(!isFinalHiggs(part)) return false;
    int nW = 0;
    for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
      if (part->child(iChild) == nullptr) continue;
      if (MC::isW(part->child(iChild)->pdgId())) nW++;
    }
    if(nW) return true;
    else return false;
  }

  bool isHiggsToZZ(const xAOD::TruthParticle *part)
  {
    if(!isFinalHiggs(part)) return false;
    int nZ = 0;
    for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
      if (part->child(iChild) == nullptr) continue;
      if (MC::isZ(part->child(iChild)->pdgId())) nZ++;
    }
    if(nZ) return true;
    else return false;
  }

  bool isHiggsToBB(const xAOD::TruthParticle *part)
  {
    if(!isFinalHiggs(part)) return false;
    int nB = 0;
    for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
      if (part->child(iChild) == nullptr) continue;
      if (abs(part->child(iChild)->pdgId())==5) nB++;
    }
    if(nB) return true;
    else return false;
  }

  TruthPtcls getFinalHiggsBosons(const xAOD::TruthParticleContainer *truthParticles)
  {
    TruthPtcls higgs(SG::VIEW_ELEMENTS);
    for (auto part: *truthParticles) {
      if (isFinalHiggs(part))
        higgs.push_back(part);
    }
    return higgs;
  }
  
  bool isFromHiggs(const xAOD::TruthParticle *ptcl) {
    if (MC::isHiggs(ptcl->pdgId())) return true;
    if (ptcl->parent()==nullptr) return false;
    return isFromHiggs(ptcl->parent());
  }

  bool isFinalNeutrino(const xAOD::TruthParticle *part)
  {
    if (!MC::isNeutrino(part->pdgId())) return false;
    if (part->child() == nullptr) return true;
    for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
      if (MC::isNeutrino(part->child(iChild)->pdgId())) return false;
    }
    return true;
  }

  bool isFinalElectron(const xAOD::TruthParticle *part)
  {
    if (!MC::isElectron(part->pdgId())) return false;
    if (part->child() == nullptr) return true;
    for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
      if (MC::isElectron(part->child(iChild)->pdgId())) return false;
    }
    return true;
  }

  bool isPairProducedElectron(const xAOD::TruthParticle *part)
  {
    if (!MC::isElectron(part->pdgId())) return false;
    if (part->parent() == nullptr) return false;
    if (MC::isElectron(part->parent()->pdgId())) return isPairProducedElectron(part->parent());
    if (MC::isPhoton(part->parent()->pdgId())) return true;
    return false;
  }

  bool isFinalMuon(const xAOD::TruthParticle *part)
  {
    if (!MC::isMuon(part->pdgId())) return false;
    if (part->child() == nullptr) return true;
    for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
      if (MC::isMuon(part->child(iChild)->pdgId())) return false;
    }
    return true;
  }

  bool isFinalTau(const xAOD::TruthParticle *part)
  {
    if (!MC::isTau(part->pdgId())) return false;
    if (part->child() == nullptr) return false;
    for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
      if (MC::isTau(part->child(iChild)->pdgId())) return false;
    }
    return true;
  }

  bool isFinalZ(const xAOD::TruthParticle *part)
  {
    if (!MC::isZ(part->pdgId())) return false;
    if (part->child() == nullptr) return false;
    for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
      if(!part->child(iChild)) continue;
      if (MC::isZ(part->child(iChild)->pdgId())) return false;
    }
    return true;
  }

  bool isInitialZ(const xAOD::TruthParticle *part)
  {
    if (!MC::isZ(part->pdgId())) return false;
    if (part->parent() == nullptr) return true;
    if (MC::isZ(part->parent()->pdgId())) return false;
    return true;
  }

  bool isFinalW(const xAOD::TruthParticle *part)
  {
    if (!MC::isW(part->pdgId())) return false;
    if (part->child() == nullptr) return false;
    for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
      if(!part->child(iChild)) continue;
      if (MC::isW(part->child(iChild)->pdgId())) return false;
    }
    return true;
  }

  bool isInitialW(const xAOD::TruthParticle *part)
  {
    if (!MC::isW(part->pdgId())) return false;
    if (part->parent() == nullptr) return true;
    if (MC::isW(part->parent()->pdgId())) return false;
    return true;
  }

  bool isWLep(const xAOD::TruthParticle *part)
  {
    if (!MC::isW(part->pdgId())) return false;
    if(!isFinalW(part)) {
      for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
        if (MC::isW(part->child(iChild)->pdgId())) return isWLep(part->child(iChild));
      }
    }
    int nLep = 0;
    for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
      if (part->child(iChild) == nullptr) continue;
      if (MC::isLepton(part->child(iChild)->pdgId())) nLep++;
    }
    if(nLep) return true;
    else return false;
  }

  bool isWHad(const xAOD::TruthParticle *part)
  {
    if (!MC::isW(part->pdgId())) return false;
    if(!isFinalW(part)) {
      for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
        if (MC::isW(part->child(iChild)->pdgId())) return isWHad(part->child(iChild));
      }
    }
    int nQ = 0;
    for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
      if (part->child(iChild) == nullptr) continue;
      if (MC::isQuark(part->child(iChild)->pdgId())) nQ++;
    }
    if(nQ) return true;
    else return false;
  }

  bool isZLep(const xAOD::TruthParticle *part)
  {
    if (!MC::isZ(part->pdgId())) return false;
    if(!isFinalZ(part)) {
      for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
        if (MC::isZ(part->child(iChild)->pdgId())) return isZLep(part->child(iChild));
      }
    }
    int nLep = 0;
    for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
      if (part->child(iChild) == nullptr) continue;
      if (MC::isLepton(part->child(iChild)->pdgId())) nLep++;
    }
    if(nLep) return true;
    else return false;
  }

  bool isZHad(const xAOD::TruthParticle *part)
  {
    if (!MC::isZ(part->pdgId())) return false;
    if(!isFinalZ(part)) {
      for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
        if (MC::isZ(part->child(iChild)->pdgId())) return isZHad(part->child(iChild));
      }
    }
    int nQ = 0;
    for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
      if (part->child(iChild) == nullptr) continue;
      if (MC::isQuark(part->child(iChild)->pdgId())) nQ++;
    }
    if(nQ) return true;
    else return false;
  }

  TruthPtcls getFinalWBosons(const xAOD::TruthParticleContainer *truthParticles)
  {
    TruthPtcls ws(SG::VIEW_ELEMENTS);
    for (auto part: *truthParticles) {
      if (isFinalW(part))
        ws.push_back(part);
    }
    return ws;
  }

  bool isFromW(const xAOD::TruthParticle *ptcl) {
    if (MC::isW(ptcl->pdgId())) return true;
    if (ptcl->parent()==nullptr) return false;
    return isFromW(ptcl->parent());
  }

  bool isFromZ(const xAOD::TruthParticle *ptcl) {
    if (MC::isZ(ptcl->pdgId())) return true;
    if (ptcl->parent()==nullptr) return false;
    return isFromZ(ptcl->parent());
  }

  bool isFromTop(const xAOD::TruthParticle *ptcl) {
    if (MC::isTop(ptcl->pdgId())) return true;
    if (ptcl->parent()==nullptr) return false;
    return isFromTop(ptcl->parent());
  }
  
  bool isFromBosonTop(const xAOD::TruthParticle *ptcl) {
    if(isFromHiggs(ptcl)) return true;
    if(isFromW(ptcl)) return true;
    if(isFromZ(ptcl)) return true;
    if(isFromTop(ptcl)) return true;
    return false;
  }

  bool isFromTau(const xAOD::TruthParticle *ptcl) {
    if (MC::isTau(ptcl->pdgId())) return true;
    if (ptcl->parent() == nullptr) return false;
    return isFromTau(ptcl->parent());
  }

  bool isFromGluon(const xAOD::TruthParticle *ptcl) {
    if (MC::isGluon(ptcl->pdgId())) return true;
    if (ptcl->parent() == nullptr) return false;
    if (MC::isGluon(ptcl->parent()->pdgId())) return true;
    if (ptcl->parent()->pdgId() != ptcl->pdgId()) return false;
    return isFromGluon(ptcl->parent());
  }

  bool isFromPhoton(const xAOD::TruthParticle *ptcl) {
    if (MC::isPhoton(ptcl->pdgId())) return true;
    if (ptcl->parent() == nullptr) return false;
    if (MC::isPhoton(ptcl->parent()->pdgId())) return true;
    if (ptcl->parent()->pdgId() != ptcl->pdgId()) return false;
    return isFromPhoton(ptcl->parent());
  }

  bool isFromBhadron(const xAOD::TruthParticle *ptcl) {
    if (MC::isBottomHadron(ptcl->pdgId())) return true;
    if (ptcl->parent()==nullptr) return false;
    return isFromBhadron(ptcl->parent());
  }

  bool isFromParticle(const xAOD::TruthParticle *ptcl, int pdgId) {
    if (std::abs(ptcl->pdgId())==pdgId) return true;
    if (ptcl->parent()==nullptr) return false;
    return isFromParticle(ptcl->parent(),pdgId);
  }

  bool isFinalQuark(const xAOD::TruthParticle *part)
  {
    if (!MC::isQuark(part->pdgId())) return false;
    /*
    if (part->nChildren() == 1) {
      if(part->child(0) != nullptr) {
        if (part->child(0)->pdgId() == part->pdgId()) return false;
        if (part->child(0)->pdgId() > 22) return false;
      }
    }
    */
    for(unsigned int iChild = 0; iChild < part->nChildren(); iChild++) {
      if(part->child(iChild) == nullptr) continue;
      if (part->child(iChild)->pdgId() == part->pdgId()) return false;
    }
    return true;
  }

  bool isInitialQuark(const xAOD::TruthParticle *part)
  {
    if (!MC::isQuark(part->pdgId())) return false;
    if (MC::isQuark(part->parent()->pdgId())) return false;
    if (MC::isGluon(part->parent()->pdgId())) return false;
    return true;
  }

  bool isInitialCorrQuark(const xAOD::TruthParticle *part)
  {
    if (!MC::isQuark(part->pdgId())) return false;
    if (MC::isGluon(part->parent()->pdgId())) return false;
    if (part->parent()->nChildren()==1 && part->nChildren()>1) return true;
    return false;
  }

  bool isDescendant(const xAOD::TruthParticle *descendant, const xAOD::TruthParticle *ancestor) {
    if(descendant->pt()==ancestor->pt() && descendant->eta()==ancestor->eta() && descendant->phi()==ancestor->phi()) return true;
    if (descendant->parent()==nullptr) return false;
    return isDescendant(descendant->parent(),ancestor);
  }

  bool isGoodTruthPhoton(const xAOD::TruthParticle *ptcl) {
    return isStable(ptcl) && MC::isPhoton(ptcl->pdgId()) && notFromHadron(ptcl);
  }

  bool isGoodTruthElectron(const xAOD::TruthParticle *ptcl) {
    return isStable(ptcl) && MC::isElectron(ptcl->pdgId()) && notFromHadron(ptcl);
  }

  bool isGoodTruthMuon(const xAOD::TruthParticle *ptcl) {
    return isStable(ptcl) && MC::isMuon(ptcl->pdgId()) && notFromHadron(ptcl);
  }

  bool isZdecayLepton(const xAOD::TruthParticle *ptcl) {
    return isStable(ptcl) && (MC::isElectron(ptcl->pdgId())||MC::isMuon(ptcl->pdgId())) && isFromZ(ptcl);
  }

  std::vector<const xAOD::TruthParticle*> getGoodTruthPhotonsOld(const xAOD::TruthParticleContainer *truthPtcls) {
    ::std::vector<const xAOD::TruthParticle*> truthPhotons;
    for (const xAOD::TruthParticle *ptcl : *truthPtcls)
      if (isGoodTruthPhoton(ptcl)) truthPhotons.push_back(ptcl);
    return truthPhotons;
  }

  //! /brief returns all stable photons that do not originate from hadrons
  TruthPtcls getGoodTruthPhotons( const xAOD::TruthParticleContainer *truthPtcls ) {
    TruthPtcls ys(SG::VIEW_ELEMENTS);
    for (auto ptcl : *truthPtcls) if ( isGoodTruthPhoton(ptcl) ) ys.push_back(ptcl);
    return ys;
  }

  //! /brief returns all stable electrons that do not originate from hadrons
  TruthPtcls getGoodTruthElectrons( const xAOD::TruthParticleContainer * truthPtcls ) {
    TruthPtcls es(SG::VIEW_ELEMENTS);
    for (auto ptcl : *truthPtcls) if (isGoodTruthElectron(ptcl)) es.push_back(ptcl);
    return es;
  }

  //! /brief returns all stable muons that do not originate from hadrons
  TruthPtcls getGoodTruthMuons( const xAOD::TruthParticleContainer * truthPtcls ) {
    TruthPtcls mus(SG::VIEW_ELEMENTS);
    for (auto ptcl : *truthPtcls) if (isGoodTruthMuon(ptcl)) mus.push_back(ptcl);
    return mus;
  }

  //! /brief returns all stable electrons or muons originating from a Z boson
  TruthPtcls getZdecayLeptons( const xAOD::TruthParticleContainer * truthPtcls ) {
    TruthPtcls ls(SG::VIEW_ELEMENTS);
    for (auto ptcl : *truthPtcls) if (isZdecayLepton(ptcl)) ls.push_back(ptcl);
    return ls;
  }

  TruthPtcls getHadronsAndTheirDecay( const xAOD::TruthParticleContainer * truthPtcls ) {
    TruthPtcls hadrons(SG::VIEW_ELEMENTS);
    for (auto ptcl : *truthPtcls) {
      if (!isStable(ptcl)) continue;
      if (isGoodTruthPhoton(ptcl)) continue;
      if (isGoodTruthElectron(ptcl)) continue;
      if (isGoodTruthMuon(ptcl)) continue;
      hadrons.push_back(ptcl);
    }
    return hadrons;
  }

  TruthPtcls getStableDecayProducts( const xAOD::TruthParticle *ptcl ) {
    TruthPtcls decay(SG::VIEW_ELEMENTS);
    if (isStable(ptcl)) { decay.push_back(ptcl); return decay; }
    for (size_t ichild=0;ichild<ptcl->nChildren();++ichild)
      if (ptcl->child(ichild)) for (auto p:getStableDecayProducts(ptcl->child(ichild))) decay.push_back(p);
    return decay;
  }


  //! /brief returns all stable electrons that do not originate from hadrons
  TruthPtcls getBHadrons( const xAOD::TruthParticleContainer * truthPtcls, double pTcut ) {
    TruthPtcls Bs(SG::VIEW_ELEMENTS);
    for (auto ptcl : *truthPtcls)
      if ( MC::isBottomHadron(ptcl->pdgId()) && (pTcut<0||ptcl->pt()>pTcut) )
        Bs.push_back(ptcl);
    return Bs;
  }

  TruthPtcls getDHadrons( const xAOD::TruthParticleContainer * truthPtcls, double pTcut ) {
    TruthPtcls Ds(SG::VIEW_ELEMENTS);
    for (auto ptcl : *truthPtcls)
      if ( MC::isCharmHadron(ptcl->pdgId()) && (pTcut<0||ptcl->pt()>pTcut) )
        Ds.push_back(ptcl);
    return Ds;
  }

  TruthPtcls getPhotonsFromHiggs( const xAOD::TruthParticleContainer *truthPtcls ) {
    TruthPtcls ys(SG::VIEW_ELEMENTS);
    for (auto ptcl : *truthPtcls)
      if ( isGoodTruthPhoton(ptcl) && isFromHiggs(ptcl) ) ys.push_back(ptcl);
    return ys;
  }

  TruthPtcls getHiggsDecayProducts( const xAOD::TruthParticleContainer *truthPtcls ) {
    TruthPtcls decay(SG::VIEW_ELEMENTS);
    for (auto ptcl : *truthPtcls)
      if ( isStable(ptcl) && isFromHiggs(ptcl) ) decay.push_back(ptcl);
    return decay;
  }


  TruthPtcls getMuonsFromBs( const xAOD::TruthParticleContainer *truthPtcls ) {
    TruthPtcls mus(SG::VIEW_ELEMENTS);
    for (auto ptcl : *truthPtcls)
      if ( isStable(ptcl) && MC::isMuon(ptcl->pdgId()) && isFromBhadron(ptcl) )
        mus.push_back(ptcl);
    return mus;
  }

  /*
  TruthParticleStruct identifyTruthParticles(xAOD::TEvent *event,
      double jet_pTcut) {

    const xAOD::TruthParticleContainer *truthParticles = nullptr;
    if (event->retrieve(truthParticles, "TruthParticle" ).isFailure())
      HWW::fatal("Cannot access TruthParticle");

    const xAOD::JetContainer *truthJets = nullptr;
    if (event->retrieve(truthJets,"AntiKt4TruthJets").isFailure())
      HWW::fatal("Cannot access AntiKt4TruthJets");

    return HWW::identifyTruthParticles(truthParticles, truthJets, jet_pTcut);
  }

  void removeTruthOverlap(DataVector<xAOD::IParticle> &photons,
      DataVector<xAOD::IParticle> &electrons,
      //DataVector<xAOD::IParticle> &muons,
      DataVector<xAOD::IParticle> &jets,
      double jet_pTcut)
  {
    for (auto tj = jets.begin(); tj != jets.end();) {
      // apply a pT cut, if requested
      if (jet_pTcut > 0 && (*tj)->pt() < jet_pTcut) {
        tj = jets.erase(tj);
        continue;
      }

      // ignore jets overlapping with good photons, electrons or muons
      if (HWW::minDRrap(*tj, photons) < 0.4) {
        tj = jets.erase(tj);
        continue;
      }
      if (HWW::minDRrap(*tj, electrons) < 0.4) {
        tj = jets.erase(tj);
        continue; // <<== WZ jets should not do this
      }
      // if (HWW::minDRrap(tj,    muonss)<0.4) continue; // ??
    }
  }

  TruthParticleStruct identifyTruthParticles(const xAOD::TruthParticleContainer *truthPtcls,
      const xAOD::JetContainer *truthJets,
      double jet_pTcut) {
    TruthParticleStruct tp;
    TruthPtcls ys  = getGoodTruthPhotons(truthPtcls);
    TruthPtcls es  = getGoodTruthElectrons(truthPtcls);
    TruthPtcls mus = getGoodTruthMuons(truthPtcls);
    TruthPtcls hads =  getHadronsAndTheirDecay(truthPtcls);

    // TO-DO
    // Dressing should happen here !
    tp.electrons = es; tp.muons = mus;
    // some ys should probably go to hads here
    tp.photons = ys; tp.hadrons = hads;

    tp.photonsFromHiggs = getPhotonsFromHiggs(truthPtcls);
    // this one might be slow ... ?
    tp.HiggsDecay = getHiggsDecayProducts(truthPtcls);

    tp.Bhadrons    = getBHadrons(truthPtcls);
    tp.Dhadrons    = getDHadrons(truthPtcls);
    tp.muonsFromBs = getMuonsFromBs(truthPtcls);

    // TruthJets jets(truthJets->begin(), truthJets->end(), SG::VIEW_ELEMENTS);
    TruthJets jets(SG::VIEW_ELEMENTS);
    TruthJets bjets(SG::VIEW_ELEMENTS);
    TruthJets cjets(SG::VIEW_ELEMENTS);
    TruthJets lightJets(SG::VIEW_ELEMENTS);

    // Here applying a 5 GeV cut for the jet labelling
    TruthPtcls Bs    = getBHadrons(truthPtcls,5.0*GeV);
    TruthPtcls Ds    = getDHadrons(truthPtcls,5.0*GeV);

    // removeTruthOverlap(tp.photons, tp.electrons, tp.muons, jets, jet_pTcut);

    for (const xAOD::Jet *tjet : *truthJets) {
      // apply a pT cut, if requested
      if ( jet_pTcut>0 && tjet->pt()<jet_pTcut ) continue;

      // ignore jets overlapping with good photons, electrons or muons
      if (HWW::minDRrap(tjet,ys)<0.4) continue;
      if (HWW::minDRrap(tjet,es)<0.4) continue; // <<== WZ jets should not do this
      // if (HWW::minDRrap(tjet,mus)<0.4) continue; ??
      jets.push_back(tjet);

      // classify all jets into b, c or light
      if      (HWW::minDRrap(tjet,Bs)<0.4) bjets.push_back(tjet);
      else if (HWW::minDRrap(tjet,Ds)<0.4) cjets.push_back(tjet);
      else lightJets.push_back(tjet);
    }

    // later: further split light jets into: LQ, gluon, unmatched
    tp.jets  = jets;
    tp.bJets = bjets;
    tp.cJets = cjets;
    tp.lightJets = lightJets;
    return tp;
  }
  */

  void printTruthParticles( const TruthParticleStruct &tp ) {
    printf("Identified truth particles:\n");
    printf("  %lu photons\n",tp.photons.size());
    printf("  %lu electrons, %lu muons\n",tp.electrons.size(),tp.muons.size());
    printf("  %lu photons from Higgs\n",tp.photonsFromHiggs.size());
    printf("  %lu B- and %lu D-hadrons\n",tp.Bhadrons.size(),tp.Dhadrons.size());
    printf("  %lu muons from B-hadrons\n",tp.muonsFromBs.size());
    printf("  %lu jets, of which\n",tp.jets.size());
    printf("  %lu b-, %lu c- and %lu light jets\n",
        tp.bJets.size(),tp.cJets.size(),tp.lightJets.size());
  }

} // namespace HWW
