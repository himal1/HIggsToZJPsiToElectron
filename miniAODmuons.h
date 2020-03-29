#ifndef _miniAODmuons_h
#define _miniAODmuons_h

// system include files
#include <memory>
#include <map>
// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "TFile.h"
#include "TTree.h"


//
// class decleration
//

class miniAODmuons : public edm::EDAnalyzer {
public:
  explicit miniAODmuons(const edm::ParameterSet&);
  ~miniAODmuons();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;

    // ----------member data ---------------------------
  
  edm::EDGetTokenT<edm::View<pat::Muon>> dimuon_Label;
  edm::EDGetTokenT<edm::View<pat::Electron>> dielectron_Label;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trakCollection_label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_Label;
  
  //trigger------------
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  //edm::EDGetTokenT<reco::GenParticleCollection> GenCollection_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;

  bool isMC_;

  TTree*      tree_;
  
  std::vector<float>       *Run, *LumiBlock, *Event;
  std::vector<float>       *FourL_mass;
  std::vector<float>       *FourL_px, *FourL_py, *FourL_pz;
  std::vector<float>       *FourL_pt, *FourL_eta, *FourL_phi;
  std::vector<float>       *FourL_VtxProb;
  std::vector<float>       *FourL_PVx,  *FourL_PVy,  *FourL_PVz;
  std::vector<float>       *FourL_PVxError,  *FourL_PVyError,  *FourL_PVzError;

  std::vector<float>        *B_Z_dca;
  std::vector<bool>        *B_Z_TriggerPath;
  std::vector<float>       *B_Z_TriggerPt1, *B_Z_TriggerEta1, *B_Z_TriggerPhi1;
  std::vector<float>       *B_Z_TriggerPt2, *B_Z_TriggerEta2, *B_Z_TriggerPhi2;
  std::vector<float>       *B_Z_TriggerPt3, *B_Z_TriggerEta3, *B_Z_TriggerPhi3;
  std::vector<float>       *B_Z_TriggerPt4, *B_Z_TriggerEta4, *B_Z_TriggerPhi4;
  std::vector<float>       *B_Z_TriggerPt5, *B_Z_TriggerEta5, *B_Z_TriggerPhi5;
  std::vector<float>       *B_Z_mass, *B_Z_VtxProb, *B_Z_px, *B_Z_py, *B_Z_pz;
  std::vector<float>       *B_Z_pt,*B_Z_eta,*B_Z_phi, *B_Z_rapidity;
  std::vector<float>       *B_Z_VtxPx, *B_Z_VtxPy, *B_Z_VtxPz;
  std::vector<float>       *B_Z_VtxPt, *B_Z_VtxEta, *B_Z_VtxPhi,  *B_Z_VtxRapidity, *B_Z_VtxMass;
  std::vector<float>       *B_Z_PVx,  *B_Z_PVy,  *B_Z_PVz;
  std::vector<float>       *B_Z_PVxError,  *B_Z_PVyError,  *B_Z_PVzError;

  std::vector<float>       *B_Z_px1, *B_Z_py1, *B_Z_pz1;
  std::vector<float>       *B_Z_pt1,*B_Z_eta1,*B_Z_phi1;
  
  std::vector<float>       *B_Z_ecalIso1, *B_Z_hcalIso1, *B_Z_trackIso1;
  std::vector<bool>        *B_Z_looseMva1; 
  std::vector<float>       *B_Z_ID_BDT1;
  std::vector<bool>        *B_Z_loose1;
  std::vector<float>       *B_Z_px2, *B_Z_py2, *B_Z_pz2;
  std::vector<float>       *B_Z_pt2,*B_Z_eta2,*B_Z_phi2;
  std::vector<float>       *B_Z_ecalIso2, *B_Z_hcalIso2, *B_Z_trackIso2;
  std::vector<bool>        *B_Z_looseMva2;
  std::vector<float>       *B_Z_ID_BDT2;
  std::vector<bool>        *B_Z_loose2;
  std::vector<int>         *B_Z_charge1, *B_Z_charge2;
  std::vector<float>       *B_Z_dxy1,  *B_Z_dxy2,  *B_Z_dz1,  *B_Z_dz2;

  std::vector<float>        *B_J_dca;

  std::vector<float>       *B_J_mass, *B_J_px, *B_J_py, *B_J_pz;
  std::vector<float>       *B_J_pt, *B_J_eta, *B_J_phi, *B_J_rapidity;
  std::vector<float>       *B_J_VtxPx, *B_J_VtxPy, *B_J_VtxPz;
  std::vector<float>       *B_J_VtxPt, *B_J_VtxEta, *B_J_VtxPhi,  *B_J_VtxRapidity, *B_J_VtxMass;
  std::vector<float>       *B_J_PVx,  *B_J_PVy,  *B_J_PVz;
  std::vector<float>       *B_J_PVxError,  *B_J_PVyError,  *B_J_PVzError;

  std::vector<float>       *B_J_px1, *B_J_py1, *B_J_pz1;
  std::vector<float>       *B_J_pt1, *B_J_eta1, *B_J_phi1;
  std::vector<bool>       *B_J_soft1,  *B_J_tight1,  *B_J_loose1;


  std::vector<float>       *B_J_px2, *B_J_py2, *B_J_pz2;
  std::vector<float>       *B_J_pt2, *B_J_eta2, *B_J_phi2;
  std::vector<int>         *B_J_charge1, *B_J_charge2;
  std::vector<bool>       *B_J_soft2,  *B_J_tight2,  *B_J_loose2;
  std::vector<float>       *B_J_VtxProb;
  std::vector<float>       *B_J_xyP,  *B_J_xyM,  *B_J_zP,  *B_J_zM;

  std::vector<float>       *mumC2;
  std::vector<int>         *mumNHits, *mumNPHits; 
  std::vector<float>       *mupC2;
  std::vector<int>         *mupNHits, *mupNPHits;
 
  unsigned int             nB;
    
};
#endif
