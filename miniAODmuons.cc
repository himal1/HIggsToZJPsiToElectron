// -*- C++ -*-
//
// Package:    miniAODmuons
// Class:      miniAODmuons
// 

//=================================================
// original author:  Jhovanny Andres Mejia        |
//         created:  Monday Aug 28 (2017)         |
//         <jhovanny.andres.mejia.guisao@cern.ch> | 
//=================================================

// system include files
#include <memory>

#include "myAnalyzers/JPsiKsPAT/src/miniAODmuons.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
//kalman vertexing
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
//GenInfo
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//trigger                                                                                                                                                                                                   
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TROOT.h"
#include "TH2F.h"

//
// constants, enums and typedefs
//

typedef math::Error<3>::type CovarianceMatrix;

//
// static data member definitions
//

//
// constructors and destructor
//

miniAODmuons::miniAODmuons(const edm::ParameterSet& iConfig)
  :
  dimuon_Label(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("dimuons"))),
  dielectron_Label(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("dielectron"))),
  trakCollection_label(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("Trak"))),
  primaryVertices_Label(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  //trigger
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  //GenLevel Info
  prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
  
  isMC_(iConfig.getParameter<bool>("isMC")),


  tree_(0), 
  Run(0),LumiBlock(0),Event(0),
  FourL_mass(0),
  FourL_px(0),FourL_py(0),FourL_pz(0),
  FourL_pt(0),FourL_eta(0),FourL_phi(0),
  FourL_VtxProb(0),
  FourL_PVx(0), FourL_PVy(0), FourL_PVz(0),
  FourL_PVxError(0), FourL_PVyError(0), FourL_PVzError(0),

  //Z 
  B_Z_TriggerPath(0),
  B_Z_TriggerPt1(0),B_Z_TriggerEta1(0),B_Z_TriggerPhi1(0),
  B_Z_TriggerPt2(0),B_Z_TriggerEta2(0),B_Z_TriggerPhi2(0),
  B_Z_TriggerPt3(0),B_Z_TriggerEta3(0),B_Z_TriggerPhi3(0),
  B_Z_TriggerPt4(0),B_Z_TriggerEta4(0),B_Z_TriggerPhi4(0),
  B_Z_TriggerPt5(0),B_Z_TriggerEta5(0),B_Z_TriggerPhi5(0),
  
  B_Z_mass(0), B_Z_VtxProb(0), 
  B_Z_px(0), B_Z_py(0), B_Z_pz(0),
  B_Z_pt(0),B_Z_eta(0),B_Z_phi(0),
  B_Z_rapidity(0), B_Z_VtxPx(0), B_Z_VtxPy(0), B_Z_VtxPz(0),
  B_Z_VtxPt(0),B_Z_VtxEta(0),B_Z_VtxPhi(0),B_Z_VtxRapidity(0),B_Z_VtxMass(0),
  
  B_Z_PVx(0), B_Z_PVy(0), B_Z_PVz(0),
  B_Z_PVxError(0), B_Z_PVyError(0), B_Z_PVzError(0),
  B_Z_px1(0), B_Z_py1(0), B_Z_pz1(0),
  B_Z_pt1(0), B_Z_eta1(0), B_Z_phi1(0),
  B_Z_ecalIso1(0),B_Z_hcalIso1(0),B_Z_trackIso1(0), B_Z_looseMva1(0), 
  B_Z_ID_BDT1(0), B_Z_loose1(0),
  B_Z_px2(0), B_Z_py2(0), B_Z_pz2(0),
  B_Z_pt2(0), B_Z_eta2(0), B_Z_phi2(0),
  B_Z_ecalIso2(0),B_Z_hcalIso2(0),B_Z_trackIso2(0), B_Z_looseMva2(0), 
  B_Z_ID_BDT2(0), B_Z_loose2(0),
  B_Z_charge1(0), B_Z_charge2(0), B_Z_dxy1(0), B_Z_dxy2(0),
  B_Z_dz1(0),B_Z_dz2(0),
  B_J_dca(0),
  B_J_mass(0), B_J_px(0), B_J_py(0), B_J_pz(0),
  B_J_pt(0),B_J_eta(0), B_J_phi(0),
  B_J_rapidity(0), B_J_VtxPx(0), B_J_VtxPy(0), B_J_VtxPz(0),
  B_J_VtxPt(0),B_J_VtxEta(0),B_J_VtxPhi(0),B_J_VtxRapidity(0),B_J_VtxMass(0),
  B_J_PVx(0), B_J_PVy(0), B_J_PVz(0),
  B_J_PVxError(0), B_J_PVyError(0), B_J_PVzError(0),
  B_J_px1(0), B_J_py1(0), B_J_pz1(0),
  B_J_pt1(0),B_J_eta1(0), B_J_phi1(0),
  B_J_soft1(0),B_J_tight1(0),B_J_loose1(0),
  B_J_px2(0), B_J_py2(0), B_J_pz2(0),
  B_J_pt2(0),B_J_eta2(0), B_J_phi2(0),
  B_J_charge1(0), B_J_charge2(0),
  B_J_soft2(0),B_J_tight2(0),B_J_loose2(0),
  B_J_VtxProb(0),
  B_J_xyP(0),  B_J_xyM(0),  B_J_zP(0),
  B_J_zM(0),
  
  mumC2(0), mumNHits(0), mumNPHits(0),
  mupC2(0), mupNHits(0), mupNPHits(0),
   
  nB(0)
  
{
   //now do what ever initialization is needed
}


miniAODmuons::~miniAODmuons()
{

}


//
// member functions
//

// ------------ method called to for each event  ------------
void miniAODmuons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;
  
  //*********************************
  // Get event content information
  //*********************************  

  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

  edm::Handle< View<pat::PackedCandidate> > thePATTrackHandle;
  iEvent.getByToken(trakCollection_label,thePATTrackHandle);


  edm::Handle< View<pat::Muon> > thePATMuonHandle;
  iEvent.getByToken(dimuon_Label,thePATMuonHandle);
  
  edm::Handle< View<pat::Electron> > thePATElectronHandle;
  iEvent.getByToken(dielectron_Label,thePATElectronHandle);
  
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);
  //gen Particles
  // Pruned particles are the one containing "important" stuff
  edm::Handle<edm::View<reco::GenParticle> > pruned;
  iEvent.getByToken(prunedGenToken_,pruned);
  //some  cross checks
  if (!theB.isValid()) {
    edm::LogWarning("miniAODmuons") << "no Transient Track in event";
    return;
  }

  if (!thePATElectronHandle.isValid()) {
    edm::LogWarning("miniAODmuons") << "no pat::Electrons in event";
    return;
  }
  if (!thePATMuonHandle.isValid()) {
    edm::LogWarning("miniAODmuons") << "no pat::Muons in event";
    return;
  }
  if (!triggerBits.isValid()) {
    edm::LogWarning("miniAODmuons") << "no Trigger path in event";
    return;
  }
  if (!triggerObjects.isValid()) {
    edm::LogWarning("miniAODmuons") << "no Trigger Object in event";
    return;
  }
  
  //Before Begining lets get Gen level Information
  int h=0; int pass = 0; int passTrig = 0;
  bool EET=false;
  float EET_pt[5] = {-999,-999,-999,-999,-999};
  float EET_eta[5] = {-999,-999,-999,-999,-999};
  float EET_phi[5] = {-999,-999,-999,-999,-999};
  int nG1=0;
  int nG2=0;
  bool GenInfo=false;
  if (GenInfo==true){
    //Gen Level Info
    /*
    float B_J_GenMuon_pt=-999;
    float B_J_GenMuon_eta=-999;
    float B_J_GenMuon_phi=-999;
    float B_Z_GenMuon_pt=-999;
    float B_Z_GenMuon_eta=-999;
    float B_Z_GenMuon_phi=-999;
    */
    //cout<<"Start Gen Work"<<endl;
    for(size_t i=0; i<pruned->size();i++){
      //if( (*pruned)[i].isPromptFinalState()  && abs((*pruned)[i].pdgId() ==13) ){
      if( abs((*pruned)[i].pdgId()) ==13 ){
	//cout<<"Found gen level muon"<<endl;
	if ( (*pruned)[i].mother()->pdgId()==443 ) {
	  //B_J_GenMuon_pt = (*pruned)[i].pt();
	  //B_J_GenMuon_eta = (*pruned)[i].eta();
	  //B_J_GenMuon_phi = (*pruned)[i].phi();
	  //B_J_GenMuonPt->push_back(B_J_GenMuon_pt);
	  //B_J_GenMuonEta->push_back(B_J_GenMuon_eta);
	  //B_J_GenMuonPhi->push_back(B_J_GenMuon_phi);
	  nG1++;
	}
      }
      if( abs((*pruned)[i].pdgId()) ==11 ){
	//cout<<"Found gen level muon"<<endl;
	if ( (*pruned)[i].mother()->pdgId()==23 ) {
	  //B_Z_GenMuon_pt = (*pruned)[i].pt();
	  //B_Z_GenMuon_eta = (*pruned)[i].eta();
	  //B_Z_GenMuon_phi = (*pruned)[i].phi();
	  //B_Z_GenMuonPt->push_back(B_Z_GenMuon_pt);
	  //B_Z_GenMuonEta->push_back(B_Z_GenMuon_eta);
	  //B_Z_GenMuonPhi->push_back(B_Z_GenMuon_phi);
	  nG2++;
	}
      }
    }
    cout<<"End Gen Work"<<endl;
  }
  
  //First Work on Trigger Information//
  //**********************************
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  //std::cout << "\n === TRIGGER PATHS === " << std::endl;
  
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    bool acceptE = triggerBits->accept(i);
    if (names.triggerName(i).find("HLT_Ele35_WPTight_Gsf_v")!=string::npos)  {
      if (triggerBits->accept(i)) {
        //std::cout << "Trigger " << names.triggerName(i) <<
        // ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
        //        << std::endl;
        pass=acceptE;
        if (pass==1){
          passTrig++;
        }
      }
    }
    
  }
  //****************************************************
  //***************Trigger Object***********************
  //****************************************************      
  if (passTrig>0){
    EET=true;
    //float ElectronTriggerPt;
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
      obj.unpackPathNames(names);
      //for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
      //bool isBoth = obj.hasPathName("HLT_Ele35_WPTight_Gsf_v*", true, true );
      //bool isL3 = obj.hasFilterLabel("HLTEle32WPTightGsfSequence");
      //cout<<"The L3 Filter is : "<<isL3<<endl;
      //bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
      //bool isL3   = obj.hasPathName( "HLT_Ele35_WPTight_Gsf_v*", false, true );
      //bool isLF   = obj.hasPathName( "HLT_Ele35_WPTight_Gsf_v*", true, false );
      //bool isNone = obj.hasPathName( "HLT_Ele35_WPTight_Gsf_v*", false, false );
      //**************************************************************************************************//
      //****************Definition of hasPathName*********************************************************//
      //**************************************************************************************************//
      //bool hasPathName(const std::string &pathName, 
      //               bool pathLastFilterAccepted = false,
      //               bool pathL3FilterAccepted = true) const {
      //return hasPathOrAlgorithm(pathName, pathLastFilterAccepted, pathL3FilterAccepted);
      //};
      //**************************************************************************************************//
      //cout<< "obj.hasPathL3FilterAccepted() "<<obj.hasPathL3FilterAccepted()<<endl; 
      //std::cout << "   " << pathNamesAll[h]; 
      if (obj.hasPathName("HLT_Ele35_WPTight_Gsf_v*", true, true )>0) {
        //std::cout << "\tTrigger objectisBoth:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
        //open from here
	EET_pt[h]=obj.pt();
        EET_eta[h]=obj.eta();
        EET_phi[h]=obj.phi();
	
      }
      //if (isL3>0) { 
      //std::cout << "\tTrigger objectisL3:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      //}
      //if (isLF>0) {  
      //      std::cout << "\tTrigger objectisLF:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
      //}
      //if (isNone>0) {
      //std::cout << "\tTrigger objectNone :  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;                                                                            
      //}
      h++;
    }  
  }                           
  
  //*********************************
  //Now we get the primary vertex 
  //*********************************
  
  reco::Vertex bestVtx;
  edm::Handle<reco::VertexCollection> primaryVertices_handle;
  iEvent.getByToken(primaryVertices_Label, primaryVertices_handle);

  bestVtx = *(primaryVertices_handle->begin());

  //nVtx = primaryVertices_handle->size(); 
  
  //*****************************************
  //Let's begin by looking for J/psi

  //unsigned int nMu_tmp = thePATMuonHandle->size();

  for(View<pat::Muon>::const_iterator iMuon1 = thePATMuonHandle->begin(); iMuon1 != thePATMuonHandle->end(); ++iMuon1) 
    {
      
      for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != thePATMuonHandle->end(); ++iMuon2) 
	{
	  if(iMuon1==iMuon2) continue;
	  
	  //opposite charge 
	  if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;
	  
	  for(View<pat::Electron>::const_iterator iEle1 = thePATElectronHandle->begin(); iEle1 != thePATElectronHandle->end(); ++iEle1) {
	    //cout<<"Begining of cut on Ele 1"<<endl;
	    for(View<pat::Electron>::const_iterator iEle2 = iEle1+1; iEle2 != thePATElectronHandle->end(); ++iEle2){
	      
	      if(iEle1==iEle2) continue;
	      //cout<<"Begining of cut on Ele 2"<<endl;
	      if( (iEle1->charge())*(iEle2->charge()) == 1) continue;
	      TrackRef glbTrackP;	  
	      TrackRef glbTrackM;	  
	      //cout<<"Electron Track Pt"<<iEle1->gsfTrack()->pt()<<endl;
	      
	      if(iMuon1->charge() == 1){ glbTrackP = iMuon1->track();}
	      if(iMuon1->charge() == -1){ glbTrackM = iMuon1->track();}
	      
	      if(iMuon2->charge() == 1) { glbTrackP = iMuon2->track();}
	      if(iMuon2->charge() == -1){ glbTrackM = iMuon2->track();}
	      
	      if( glbTrackP.isNull() || glbTrackM.isNull() ) 
		{
		  //std::cout << "continue due to no track ref" << endl;
		  continue;
		}
	      
	      
	      reco::TrackRef kfTrackRefP;
	      reco::TrackRef kfTrackRefM;
	      //cout<<"Start Electron Analysis"<<endl;
	      if (! (abs(iEle1->charge())==1)) continue;
	      if (! (abs(iEle2->charge())==1)) continue;
	      //if (iEle->charge()==0) continue;
	      /*
	      if(iEle1->charge() == 1){ kfTrackRefP = iEle1->closestCtfTrackRef();}
	      if(iEle1->charge() == -1){ kfTrackRefM = iEle1->closestCtfTrackRef();}
	     
	      if(iEle2->charge() == 1){ kfTrackRefP = iEle2->closestCtfTrackRef();}
	      if(iEle2->charge() == -1){ kfTrackRefM = iEle2->closestCtfTrackRef();}
	      */
	      
	      TLorentzVector M1,M2,E1,E2,MM,EE,EEMM;
	      float mu_mass = 0.1056583745;//[PDG mass]
	      float ele_mass =  0.000510998928;//PDG mass
	      M1.SetXYZM(iMuon1->track()->px(),iMuon1->track()->py(),iMuon1->track()->pz(),mu_mass);
	      M2.SetXYZM(iMuon2->track()->px(),iMuon2->track()->py(),iMuon2->track()->pz(),mu_mass);
	      //E1.SetXYZM(iEle1->gsfTrack()->px(),iEle1->gsfTrack()->py(),iEle1->gsfTrack()->pz(),ele_mass);
	      //E2.SetXYZM(iEle2->gsfTrack()->px(),iEle2->gsfTrack()->py(),iEle2->gsfTrack()->pz(),ele_mass);
	      E1.SetPtEtaPhiE(iEle1->pt(),iEle1->eta(),iEle1->phi(),iEle1->energy());
	      E2.SetPtEtaPhiE(iEle2->pt(),iEle2->eta(),iEle2->phi(),iEle2->energy());
	      MM=M1+M2;
	      EE=E1+E2;
	      EEMM=MM+EE;
	      //if (MM.M()<2.8) continue;
	      //if (MM.M()>3.4) continue;
	      //if (EE.M()<70) continue;
	      //if (EE.M()>110) continue;
	      
	      int tkquality=0;
	      /*	      
	      if (kfTrackRefP.isAvailable() && kfTrackRefP.isNonnull()) {
		//cout<<"TRACK REFERENCE IS AVILABLE FOR FIRST ELECTRON"<<endl;
		if (kfTrackRefM.isAvailable() && kfTrackRefM.isNonnull()) {
		  tkquality++;
		}
	      }
	      */
	      if (iEle1->gsfTrack().isAvailable()&&iEle1->gsfTrack().isNonnull()) {
		//cout<<"TRACK REFERENCE IS AVILABLE FOR FIRST ELECTRON"<<endl;
		if (iEle2->gsfTrack().isAvailable() && iEle2->gsfTrack().isNonnull()) {
		  tkquality++;
		}
	      }
	      
	      if (tkquality==0) continue;
	      //cout<<"CTF track quality "<<tkquality<<endl;
	      //if (tkquality<1)continue;
	      // cout<<"CTF track quality "<<tkquality<<endl;
	      if(iMuon1->track()->pt()<3.0) continue;
	      if(iMuon2->track()->pt()<3.0) continue;
	      	      
	      //cout<<"Start looking muon track quality"<<endl;
	      if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
	      if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;	 
	      //cout<<"Start Building Track"<<endl;
	      
	      reco::TransientTrack muon1TT((*theB).build(glbTrackP));
	      reco::TransientTrack muon2TT((*theB).build(glbTrackM));
	      //reco::TransientTrack electron1TT((*theB).build(kfTrackRefP));//working
	      //reco::TransientTrack electron2TT((*theB).build(kfTrackRefM));//working
	      reco::TransientTrack electron1TT((*theB).build(iEle1->gsfTrack()));
	      reco::TransientTrack electron2TT((*theB).build(iEle2->gsfTrack()));
	      
	      // *****  Trajectory states to calculate DCA for the 2 muons *********************
	      FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
	      FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();
	      //FreeTrajectoryState electron1State = muon1TT.impactPointTSCP().theState();
	      //FreeTrajectoryState electron2State = muon2TT.impactPointTSCP().theState();
	      //cout<<"Start validating impact point"<<endl;
	      if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;
	      //cout<<"Validated  impact point for muon"<<endl;
	      if( !electron1TT.impactPointTSCP().isValid() || !electron2TT.impactPointTSCP().isValid() ) continue;
	      //cout<<"Start validating impact point for electron"<<endl;
	      // Measure distance between tracks at their closest approach
	      ClosestApproachInRPhi cApp;
	      cApp.calculate(mu1State, mu2State);
	      if( !cApp.status() ) continue;
	      float dca = fabs( cApp.distance() );	  
	      //if (dca < 0. || dca > 0.5) continue;
	      //cout<<"dca"<<dca<<endl;
	      //cout<<" closest approach  "<<dca<<endl;
	      //ClosestApproachInRPhi cAppE;
	      //cApp.calculate(electron1State, electron2State);
	      //if( !cAppE.status() ) continue;
	      //float dcaE = fabs( cAppE.distance() );	  
	      //if (dcaE < 0. || dcaE > 0.5) continue;
	      //cout<<"dca E "<<dcaE<<endl;
	      // ******  Methods to check to which category of muon candidates a given pat::Muon object belongs ****
	      //Kalman Vtx----------------------//
	      
	      vector <TransientTrack> ele_tks; 
	      vector <TransientTrack> mu_tks;
	      vector <TransientTrack> mmee_tks;
	      KalmanVertexFitter kvf(true);
	      ele_tks.clear();
	      ele_tks.push_back(electron1TT); 
	      ele_tks.push_back(electron2TT);
	  
	      TransientVertex Z_candi = kvf.vertex(ele_tks);
	      reco::Vertex Z_Vtx = Z_candi;
	      const math::XYZTLorentzVectorD Z_mom = Z_Vtx.p4(ele_mass,0.0);
	      if (!Z_candi.isValid()){
		//cout<<"Z candidate non validated by kalman fitter"<<endl; 
		continue; 
	      }
	      float B_Prob_tmp1 = TMath::Prob(Z_candi.totalChiSquared(),Z_candi.degreesOfFreedom());
	      if (B_Prob_tmp1 < 0.001) continue;                                                                            
	      
	      KalmanVertexFitter kvfM(true);
	      mu_tks.clear();
	      mu_tks.push_back(muon1TT); 
	      mu_tks.push_back(muon2TT);
	      TransientVertex J_candi = kvfM.vertex(mu_tks);
	      if (!J_candi.isValid()){
		//cout<<"J candidate non validated by kalman fitter"<<endl;
		continue;
	      }
	      reco::Vertex JPsi_Vtx = J_candi;
	      float B_Prob_tmp = TMath::Prob(J_candi.totalChiSquared(),J_candi.degreesOfFreedom());
	      const math::XYZTLorentzVectorD JPsi_mom = JPsi_Vtx.p4(mu_mass,0.0);
	      
	      // const math::XYZTLorentzVectorD jpsi_mom = jpsi_vtx.p4(0.1056583,0.0);

	      if (B_Prob_tmp < 0.001) continue;
	      KalmanVertexFitter kvfEM(true);
	      mmee_tks.clear();
	      mmee_tks.push_back(electron1TT);
	      mmee_tks.push_back(electron2TT);
	      mmee_tks.push_back(muon1TT);
	      mmee_tks.push_back(muon2TT);
	      TransientVertex FourL_candi = kvfEM.vertex(mmee_tks);
	      if (!FourL_candi.isValid()){
		//cout<<"J candidate non validated by kalman fitter"<<endl;
		continue;
	      }         
	      float B_Prob_tmp4L = TMath::Prob(FourL_candi.totalChiSquared(),FourL_candi.degreesOfFreedom());
	      reco::Vertex FourL_Vtx = FourL_candi;
     
	      if (B_Prob_tmp4L < 0.001) continue;
	      	  
	      /*
	      //if (iMuon1->isTrackerMuon() || iMuon2->isTrackerMuon())
	      //if (muon::isHighPtMuon(*iMuon1,bestVtx) || muon::isHighPtMuon(*iMuon2,bestVtx))
	      if (muon::isGoodMuon(*iMuon1,muon::TMLastStationAngTight) || muon::isGoodMuon(*iMuon2,muon::TMLastStationAngTight))
	      {
	      cout<<" is category muon  "<<endl;
	      }
	      else
	      {
	      cout<<" it is not category muon  "<<endl;
	      }
	      */
	      
	      // ******   Let's check the vertex and mass ****
	      
	      
	      //The mass of a muon and the insignificant mass sigma 
	      //to avoid singularities in the covariance matrix.
	      
	      bool KinFit=false;
	      
	      //cout<<"Start Kin Loop"<<endl;
	      
	      if (KinFit==true){
		ParticleMass muon_mass = 0.10565837; //pdg mass
		//ParticleMass psi_mass = 3.096916;
		float muon_sigma = muon_mass*1.e-6;
		//float psi_sigma = psi_mass*1.e-6;
		vector<RefCountedKinematicParticle> muonParticles;
		//Creating a KinematicParticleFactory
		KinematicParticleFactoryFromTransientTrack pFactory;
		
		//initial chi2 and ndf before kinematic fits.
		float chi = 0.;
		float ndf = 0.;
		//vector<RefCountedKinematicParticle> muonParticles;
		try {
		  muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
		  muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
		}
		catch(...) { 
		  std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
		  continue;
		}
		
		KinematicParticleVertexFitter fitter;   
		
		RefCountedKinematicTree psiVertexFitTree;
		try {
		  psiVertexFitTree = fitter.fit(muonParticles); 
		}
		catch (...) { 
		  std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
		  continue;
		}
		
		if (!psiVertexFitTree->isValid()) 
		  {
		//std::cout << "caught an exception in the psi vertex fit" << std::endl;
		    continue; 
		  }
		
		psiVertexFitTree->movePointerToTheTop();
		
		RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
		RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();
		
		if( psi_vFit_vertex_noMC->chiSquared() < 0 )
		  {
		    //std::cout << "negative chisq from psi fit" << endl;
		    continue;
		  }
		
		//some loose cuts go here
		
		if(psi_vFit_vertex_noMC->chiSquared()>50.) continue;
		if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.3) continue;
		
		//fill variables?iMuon1->track()->pt()
		
		//B_J_mass->push_back( psi_vFit_noMC->currentState().mass() );
		//B_J_px->push_back( psi_vFit_noMC->currentState().globalMomentum().x() );
		//B_J_py->push_back( psi_vFit_noMC->currentState().globalMomentum().y() );
		//B_J_pz->push_back( psi_vFit_noMC->currentState().globalMomentum().z() );
	      }
	      //cout<<JPsi_mom.mass()<<" : Is JPsi Vtx Mass "<<MM.M()<<"Is JPsi mass"<<endl;
	      //cout<<Z_mom.mass()<<" : Is Z Vtx Mass "<<EE.M()<<"Is Z mass"<<endl;
	      //Event Information
	      Run->push_back(iEvent.id().run());
	      LumiBlock->push_back(iEvent.luminosityBlock());
	      Event->push_back(iEvent.id().event());
	      
	      //Four Muon Information
	      
	      FourL_mass->push_back(EEMM.M());
	      FourL_px->push_back(EEMM.Px());
	      FourL_py->push_back(EEMM.Py());
	      FourL_pz->push_back(EEMM.Pz());
	      FourL_pt->push_back(EEMM.Pt());
	      FourL_eta->push_back(EEMM.Eta());
	      FourL_phi->push_back(EEMM.Phi());
	      FourL_VtxProb->push_back(B_Prob_tmp4L);
	      FourL_PVx->push_back(FourL_Vtx.x());
	      FourL_PVy->push_back(FourL_Vtx.y());
	      FourL_PVz->push_back(FourL_Vtx.z());
	      FourL_PVxError->push_back(FourL_Vtx.xError());
	      FourL_PVyError->push_back(FourL_Vtx.yError());
	      FourL_PVzError->push_back(FourL_Vtx.zError());
	      //muonParticles.clear();//open for kinematicFit                         
	      
	      //B_Z_dca->push_back(dcaM);
	      B_Z_TriggerPath->push_back(EET);
	      
	      B_Z_TriggerPt1->push_back(EET_pt[0]);
	      B_Z_TriggerEta1->push_back(EET_eta[0]);
	      B_Z_TriggerPhi1->push_back(EET_phi[0]);
	      B_Z_TriggerPt2->push_back(EET_pt[1]);
	      B_Z_TriggerEta2->push_back(EET_eta[1]);
	      B_Z_TriggerPhi2->push_back(EET_phi[1]);
	      B_Z_TriggerPt3->push_back(EET_pt[2]);
	      B_Z_TriggerEta3->push_back(EET_eta[2]);
	      B_Z_TriggerPhi3->push_back(EET_phi[2]);
	      B_Z_TriggerPt4->push_back(EET_pt[3]);
	      B_Z_TriggerEta4->push_back(EET_eta[3]);
	      B_Z_TriggerPhi4->push_back(EET_phi[3]);
	      B_Z_TriggerPt5->push_back(EET_pt[4]);
	      B_Z_TriggerEta5->push_back(EET_eta[4]);
	      B_Z_TriggerPhi5->push_back(EET_phi[4]);
	      
	      
	      B_Z_mass->push_back(EE.M() );
	      B_Z_VtxProb->push_back(B_Prob_tmp1);
	      
	      B_Z_px->push_back(EE.Px() );
	      B_Z_py->push_back(EE.Py() );
	      B_Z_pz->push_back(EE.Pz() );
	      B_Z_pt->push_back(EE.Pt() );
	      B_Z_eta->push_back(EE.Eta() );
	      B_Z_phi->push_back(EE.Phi() );
	      B_Z_VtxPx->push_back(Z_mom.Px() );
	      B_Z_VtxPy->push_back(Z_mom.Py() );
	      B_Z_VtxPz->push_back(Z_mom.Pz() );
	      B_Z_VtxPt->push_back(Z_mom.Pt() );
	      B_Z_VtxEta->push_back(Z_mom.Eta() );
	      B_Z_VtxPhi->push_back(Z_mom.Phi() );
	      B_Z_VtxRapidity->push_back(Z_mom.Rapidity() );
	      B_Z_VtxMass->push_back(Z_mom.mass() );
	      
	      B_Z_PVx->push_back(Z_Vtx.x());
	      B_Z_PVy->push_back(Z_Vtx.y());
	      B_Z_PVz->push_back(Z_Vtx.z());
	      B_Z_PVxError->push_back(Z_Vtx.xError());
	      B_Z_PVyError->push_back(Z_Vtx.yError());
	      B_Z_PVzError->push_back(Z_Vtx.zError());
	      
	      B_Z_px1->push_back(iEle1->px());
	      B_Z_py1->push_back(iEle1->py());
	      B_Z_pz1->push_back(iEle1->pz());
	      B_Z_pt1->push_back(iEle1->pt() );
	      B_Z_eta1->push_back(iEle1->eta() );
	      B_Z_phi1->push_back(iEle1->phi() );
	      B_Z_ecalIso1->push_back(iEle1->ecalIso());
	      B_Z_hcalIso1->push_back(iEle1->hcalIso());
	      B_Z_trackIso1->push_back(iEle1->trackIso());
	      B_Z_looseMva1->push_back(iEle1->electronID("mvaEleID-Spring15-25ns-nonTrig-V1-wpLoose"));
	      B_Z_ID_BDT1->push_back(iEle1->userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"));
	      B_Z_loose1->push_back(iEle1->electronID("eidLoose"));
	      B_Z_charge1->push_back(iEle1->charge());
	      
	      B_Z_px2->push_back(iEle2->px());
	      B_Z_py2->push_back(iEle2->py());
	      B_Z_pz2->push_back(iEle2->pz());

	      B_Z_pt2->push_back(iEle2->pt() );
	      B_Z_eta2->push_back(iEle2->eta() );
	      B_Z_phi2->push_back(iEle2->phi() );
	      B_Z_ecalIso2->push_back(iEle2->ecalIso());
	      B_Z_hcalIso2->push_back(iEle2->hcalIso());
	      B_Z_trackIso2->push_back(iEle2->trackIso());
	      B_Z_looseMva2->push_back(iEle2->electronID("mvaEleID-Spring15-25ns-nonTrig-V1-wpLoose"));
	      B_Z_ID_BDT2->push_back(iEle2->userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"));
	      
	      //cout<<"BDT value of the Ele"<<iEle2->userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values")<<endl;
	      B_Z_loose2->push_back(iEle2->electronID("eidLoose"));
	      B_Z_charge2->push_back(iEle2->charge());
	      B_Z_dxy1->push_back(iEle1->gsfTrack()->dxy(bestVtx.position()));
	      B_Z_dz1->push_back(iEle1->gsfTrack()->dz(bestVtx.position()));
	      B_Z_dxy2->push_back(iEle2->gsfTrack()->dxy(bestVtx.position()));
	      B_Z_dz2->push_back(iEle2->gsfTrack()->dz(bestVtx.position()));
	      
	      B_J_dca->push_back(dca);
	      B_J_mass->push_back(MM.M() );
	      B_J_px->push_back(MM.Px() );
	      B_J_py->push_back(MM.Py() );
	      B_J_pz->push_back(MM.Pz() );

	      B_J_pt->push_back(MM.Pt() );
	      B_J_eta->push_back(MM.Eta() );
	      B_J_phi->push_back(MM.Phi() );
	      B_J_rapidity->push_back(MM.Rapidity() );
	      
	      B_J_VtxPx->push_back(JPsi_mom.Px() );
	      B_J_VtxPy->push_back(JPsi_mom.Py() );
	      B_J_VtxPz->push_back(JPsi_mom.Pz() );
	      B_J_VtxPt->push_back(JPsi_mom.Pt() );
	      B_J_VtxEta->push_back(JPsi_mom.Eta() );
	      B_J_VtxPhi->push_back(JPsi_mom.Phi() );
	      B_J_VtxRapidity->push_back(JPsi_mom.Rapidity() );
	      B_J_VtxMass->push_back(JPsi_mom.mass() );
	      B_J_PVx->push_back(JPsi_Vtx.x());
	      B_J_PVy->push_back(JPsi_Vtx.y());
	      B_J_PVz->push_back(JPsi_Vtx.z());
	      B_J_PVxError->push_back(JPsi_Vtx.xError());
	      B_J_PVyError->push_back(JPsi_Vtx.yError());
	      B_J_PVzError->push_back(JPsi_Vtx.zError());


	      B_J_px1->push_back(iMuon1->track()->px());
	      B_J_py1->push_back(iMuon1->track()->py());
	      B_J_pz1->push_back(iMuon1->track()->pz());
	      B_J_pt1->push_back(iMuon1->track()->pt() );
	      B_J_eta1->push_back(iMuon1->track()->eta() );
	      B_J_phi1->push_back(iMuon1->track()->phi() );
	      B_J_charge1->push_back(iMuon1->charge());
	      B_J_soft1->push_back(iMuon1->isSoftMuon(bestVtx) );
	      B_J_tight1->push_back(iMuon1->isTightMuon(bestVtx) );
	      B_J_loose1->push_back(muon::isLooseMuon(*iMuon1) );

	      B_J_px2->push_back(iMuon2->track()->px());
	      B_J_py2->push_back(iMuon2->track()->py());
	      B_J_pz2->push_back(iMuon2->track()->pz());
	      B_J_pt2->push_back(iMuon2->track()->pt() );
	      B_J_eta2->push_back(iMuon2->track()->eta() );
	      B_J_phi2->push_back(iMuon2->track()->phi() );
	      B_J_charge2->push_back(iMuon2->charge());
	      B_J_VtxProb->push_back(B_Prob_tmp);
	      B_J_soft2->push_back(iMuon2->isSoftMuon(bestVtx) );
	      B_J_tight2->push_back(iMuon2->isTightMuon(bestVtx) );
	      B_J_loose2->push_back(muon::isLooseMuon(*iMuon2) );
	      B_J_xyP->push_back(glbTrackP->dxy(bestVtx.position()) );
	      B_J_xyM->push_back(glbTrackM->dxy(bestVtx.position()) );
	      B_J_zP->push_back(glbTrackM->dz(bestVtx.position()) );
	      B_J_zM->push_back(glbTrackP->dz(bestVtx.position()) );
	      

	      
	      //cout<<"End of all loop"<<endl;
	      	      
	      // ************
	  
	      mumC2->push_back( glbTrackP->normalizedChi2() );
	      //mumAngT->push_back( muon::isGoodMuon(*iMuon1,muon::TMLastStationAngTight) ); // 
	      mumNHits->push_back( glbTrackP->numberOfValidHits() );
	      mumNPHits->push_back( glbTrackP->hitPattern().numberOfValidPixelHits() );	       
	      mupC2->push_back( glbTrackM->normalizedChi2() );
	      //mupAngT->push_back( muon::isGoodMuon(*iMuon2,muon::TMLastStationAngTight) );  // 
	      mupNHits->push_back( glbTrackM->numberOfValidHits() );
	      mupNPHits->push_back( glbTrackM->hitPattern().numberOfValidPixelHits() );
	      
	      nB++;	       
	      if (KinFit==true){
		//muonParticles.clear();
	      }
	    }
	  }
	}
    }
  
  
  if (nB > 0 ) 
    {
      
      //std::cout << "filling tree" << endl;
      tree_->Fill();
    }
  
  nB = 0; 
  Run->clear();LumiBlock->clear();Event->clear();
  FourL_mass->clear(); FourL_px->clear();FourL_py->clear();FourL_pz->clear();
  FourL_pt->clear();FourL_eta->clear(); FourL_phi->clear(); FourL_VtxProb->clear();
  FourL_PVx->clear(); FourL_PVy->clear(); FourL_PVz->clear();
  FourL_PVxError->clear(); FourL_PVyError->clear(); FourL_PVzError->clear();
  
  //B_Z_dca->clear();
  B_Z_TriggerPath->clear();
  B_Z_TriggerPt1->clear();B_Z_TriggerEta1->clear();B_Z_TriggerPhi1->clear();
  B_Z_TriggerPt2->clear();B_Z_TriggerEta2->clear();B_Z_TriggerPhi2->clear();
  B_Z_TriggerPt3->clear();B_Z_TriggerEta3->clear();B_Z_TriggerPhi3->clear();
  B_Z_TriggerPt4->clear();B_Z_TriggerEta4->clear();B_Z_TriggerPhi4->clear();
  B_Z_TriggerPt5->clear();B_Z_TriggerEta5->clear();B_Z_TriggerPhi5->clear();
  
  B_Z_mass->clear();  B_Z_VtxProb->clear(); 
  B_Z_px->clear();
  B_Z_py->clear();  B_Z_pz->clear();
  B_Z_pt->clear();   B_Z_eta->clear();  B_Z_phi->clear(); 
  B_Z_rapidity->clear(); B_Z_VtxPx->clear();
  B_Z_VtxPy->clear();  B_Z_VtxPz->clear();
  B_Z_VtxPt->clear();   B_Z_VtxEta->clear();  B_Z_VtxPhi->clear(); B_Z_VtxMass->clear();
  B_Z_PVx->clear(); B_Z_PVy->clear(); B_Z_PVz->clear();
  B_Z_PVxError->clear(); B_Z_PVyError->clear(); B_Z_PVzError->clear();
  B_Z_VtxRapidity->clear();
  B_Z_px1->clear();  B_Z_py1->clear();  B_Z_pz1->clear(); B_Z_charge1->clear();
  B_Z_pt1->clear();  B_Z_eta1->clear(); B_Z_phi1->clear();
  B_Z_ecalIso1->clear(); B_Z_hcalIso1->clear(); B_Z_trackIso1->clear(); B_Z_looseMva1->clear(); 
  B_Z_ID_BDT1->clear(); B_Z_loose1->clear();
  B_Z_px2->clear();  B_Z_py2->clear();  B_Z_pz2->clear(); B_Z_charge2->clear();
  B_Z_pt2->clear();  B_Z_eta2->clear(); B_Z_phi2->clear();
  B_Z_ecalIso2->clear(); B_Z_hcalIso2->clear(); B_Z_trackIso2->clear(); B_Z_looseMva2->clear(); 
  B_Z_ID_BDT2->clear(); B_Z_loose2->clear();
  B_Z_dxy1->clear(); B_Z_dxy2->clear();B_Z_dz1->clear(); B_Z_dz2->clear();

  B_J_dca->clear();
  B_J_mass->clear(); B_J_px->clear();   B_J_py->clear();  B_J_pz->clear();
  B_J_pt->clear();   B_J_eta->clear();  B_J_phi->clear();
  B_J_rapidity->clear(); B_J_VtxPx->clear();
  B_J_VtxPy->clear();  B_J_VtxPz->clear();
  B_J_VtxPt->clear();   B_J_VtxEta->clear();  B_J_VtxPhi->clear(); B_J_VtxMass->clear();
  B_J_PVx->clear(); B_J_PVy->clear(); B_J_PVz->clear();
  B_J_PVxError->clear(); B_J_PVyError->clear(); B_J_PVzError->clear();
  B_J_px1->clear();  B_J_py1->clear();  B_J_pz1->clear(); B_J_charge1->clear();
  B_J_pt1->clear();  B_J_eta1->clear(); B_J_phi1->clear();
  B_J_soft1->clear(); B_J_tight1->clear(); B_J_loose1->clear();
  B_J_px2->clear();  B_J_py2->clear();  B_J_pz2->clear(); B_J_charge2->clear();
  B_J_pt2->clear();  B_J_eta2->clear(); B_J_phi2->clear(); B_J_VtxProb->clear();
  B_J_soft2->clear(); B_J_tight2->clear(); B_J_loose2->clear();
  B_J_xyP->clear(); B_J_xyM->clear(); B_J_zP->clear(); B_J_zM->clear();
  
  mumC2->clear();
  mumNHits->clear(); mumNPHits->clear();
  mupC2->clear();
  mupNHits->clear(); mupNPHits->clear();
  
}


// ------------ method called once each job just before starting event loop  ------------

void 
miniAODmuons::beginJob()
{

  std::cout << "Beginning analyzer job with value of isMC= " << isMC_ << std::endl;
  
  //edm::Service<TFileService> fs;
  //tree_ = fs->make<TTree>("ntuple"," J/psi ntuple");
  
  //tree_->Branch("nB",&nB,"nB/i");
  tree_ = new TTree("ntuple", "ntuple");

  tree_->Branch("nB",&nB,"nB/i");
  //electron channels 
  tree_->Branch("Run", &Run);
  tree_->Branch("LumiBlock", &LumiBlock);
  tree_->Branch("Event", &Event);
  tree_->Branch("FourL_mass", &FourL_mass);
  tree_->Branch("FourL_px", &FourL_px);
  tree_->Branch("FourL_py", &FourL_py);
  tree_->Branch("FourL_pz", &FourL_pz);
  tree_->Branch("FourL_pt", &FourL_pt);
  tree_->Branch("FourL_eta", &FourL_eta);
  tree_->Branch("FourL_phi", &FourL_phi);
  tree_->Branch("FourL_VtxProb", &FourL_VtxProb);
  tree_->Branch("FourL_PVx", &FourL_PVx);
  tree_->Branch("FourL_PVy", &FourL_PVy);
  tree_->Branch("FourL_PVz", &FourL_PVz);
  tree_->Branch("FourL_PVxError", &FourL_PVxError);
  tree_->Branch("FourL_PVyError", &FourL_PVyError);
  tree_->Branch("FourL_PVzError", &FourL_PVzError);

  //tree_->Branch("B_Z_dca", &B_Z_dca);
  tree_->Branch("B_Z_TriggerPath", &B_Z_TriggerPath);
  tree_->Branch("B_Z_TriggerPt1", &B_Z_TriggerPt1);
  tree_->Branch("B_Z_TriggerEta1", &B_Z_TriggerEta1);
  tree_->Branch("B_Z_TriggerPhi1", &B_Z_TriggerPhi1);
  tree_->Branch("B_Z_TriggerPt2", &B_Z_TriggerPt2);
  tree_->Branch("B_Z_TriggerEta2", &B_Z_TriggerEta2);
  tree_->Branch("B_Z_TriggerPhi2", &B_Z_TriggerPhi2);
  tree_->Branch("B_Z_TriggerPt3", &B_Z_TriggerPt3);
  tree_->Branch("B_Z_TriggerEta3", &B_Z_TriggerEta3);
  tree_->Branch("B_Z_TriggerPhi3", &B_Z_TriggerPhi3);
  tree_->Branch("B_Z_TriggerPt4", &B_Z_TriggerPt4);
  tree_->Branch("B_Z_TriggerEta4", &B_Z_TriggerEta4);
  tree_->Branch("B_Z_TriggerPhi4", &B_Z_TriggerPhi4);
  tree_->Branch("B_Z_TriggerPt5", &B_Z_TriggerPt5);
  tree_->Branch("B_Z_TriggerEta5", &B_Z_TriggerEta5);
  tree_->Branch("B_Z_TriggerPhi5", &B_Z_TriggerPhi5);
  
  tree_->Branch("B_Z_mass", &B_Z_mass);
  tree_->Branch("B_Z_VtxProb", &B_Z_VtxProb);
  
  tree_->Branch("B_Z_px", &B_Z_px);
  tree_->Branch("B_Z_py", &B_Z_py);
  tree_->Branch("B_Z_pz", &B_Z_pz);
  tree_->Branch("B_Z_pt", &B_Z_pt);
  tree_->Branch("B_Z_eta", &B_Z_eta);
  tree_->Branch("B_Z_phi", &B_Z_phi);
  tree_->Branch("B_Z_rapidity", &B_Z_rapidity);
  tree_->Branch("B_Z_VtxPx", &B_Z_VtxPx);
  tree_->Branch("B_Z_VtxPy", &B_Z_VtxPy);
  tree_->Branch("B_Z_VtxPz", &B_Z_VtxPz);
  tree_->Branch("B_Z_VtxPt", &B_Z_VtxPt);
  tree_->Branch("B_Z_VtxEta", &B_Z_VtxEta);
  tree_->Branch("B_Z_VtxPhi", &B_Z_VtxPhi);
  tree_->Branch("B_Z_VtxRapidity", &B_Z_VtxRapidity);
  tree_->Branch("B_Z_VtxMass", &B_Z_VtxMass);

  
  tree_->Branch("B_Z_PVx", &B_Z_PVx);
  tree_->Branch("B_Z_PVy", &B_Z_PVy);
  tree_->Branch("B_Z_PVz", &B_Z_PVz);  
  tree_->Branch("B_Z_PVxError", &B_Z_PVxError);
  tree_->Branch("B_Z_PVyError", &B_Z_PVyError);
  tree_->Branch("B_Z_PVzError", &B_Z_PVzError);

  tree_->Branch("B_Z_px1", &B_Z_px1);
  tree_->Branch("B_Z_py1", &B_Z_py1);
  tree_->Branch("B_Z_pz1", &B_Z_pz1);
  tree_->Branch("B_Z_pt1", &B_Z_pt1);
  tree_->Branch("B_Z_eta1", &B_Z_eta1);
  tree_->Branch("B_Z_phi1", &B_Z_phi1);
  tree_->Branch("B_Z_ecalIso1", &B_Z_ecalIso1);
  tree_->Branch("B_Z_hcalIso1", &B_Z_hcalIso1);
  tree_->Branch("B_Z_trackIso1", &B_Z_trackIso1);
  tree_->Branch("B_Z_looseMva1", &B_Z_looseMva1);
  tree_->Branch("B_Z_ID_BDT1", &B_Z_ID_BDT1);
  tree_->Branch("B_Z_loose1", &B_Z_loose1);
  tree_->Branch("B_Z_charge1", &B_Z_charge1);
  tree_->Branch("B_Z_px2", &B_Z_px2);
  tree_->Branch("B_Z_py2", &B_Z_py2);
  tree_->Branch("B_Z_pz2", &B_Z_pz2);
  tree_->Branch("B_Z_pt2", &B_Z_pt2);
  tree_->Branch("B_Z_eta2", &B_Z_eta2);
  tree_->Branch("B_Z_phi2", &B_Z_phi2);
  tree_->Branch("B_Z_ecalIso2", &B_Z_ecalIso2);
  tree_->Branch("B_Z_hcalIso2", &B_Z_hcalIso2);
  tree_->Branch("B_Z_trackIso2", &B_Z_trackIso2);
  tree_->Branch("B_Z_looseMva2", &B_Z_looseMva2);
  tree_->Branch("B_Z_ID_BDT2", &B_Z_ID_BDT2);
  tree_->Branch("B_Z_loose2", &B_Z_loose2);
  tree_->Branch("B_Z_charge2", &B_Z_charge2);
  tree_->Branch("B_Z_dxy1", &B_Z_dxy1);
  tree_->Branch("B_Z_dxy2", &B_Z_dxy2);
  tree_->Branch("B_Z_dz1", &B_Z_dz1);
  tree_->Branch("B_Z_dz2", &B_Z_dz2);
  
  tree_->Branch("B_J_dca", &B_J_dca);
  tree_->Branch("B_J_mass", &B_J_mass);
  tree_->Branch("B_J_px", &B_J_px);
  tree_->Branch("B_J_py", &B_J_py);
  tree_->Branch("B_J_pz", &B_J_pz);
  tree_->Branch("B_J_pt", &B_J_pt);
  tree_->Branch("B_J_eta", &B_J_eta);
  tree_->Branch("B_J_phi", &B_J_phi);
  tree_->Branch("B_J_rapidity", &B_J_rapidity);
  tree_->Branch("B_J_VtxPx", &B_J_VtxPx);
  tree_->Branch("B_J_VtxPy", &B_J_VtxPy);
  tree_->Branch("B_J_VtxPz", &B_J_VtxPz);
  tree_->Branch("B_J_VtxPt", &B_J_VtxPt);
  tree_->Branch("B_J_VtxEta", &B_J_VtxEta);
  tree_->Branch("B_J_VtxPhi", &B_J_VtxPhi);
  tree_->Branch("B_J_VtxRapidity", &B_J_VtxRapidity);
  tree_->Branch("B_J_VtxMass", &B_J_VtxMass);
  tree_->Branch("B_J_PVx", &B_J_PVx);
  tree_->Branch("B_J_PVy", &B_J_PVy);
  tree_->Branch("B_J_PVz", &B_J_PVz);
  tree_->Branch("B_J_PVxError", &B_J_PVxError);
  tree_->Branch("B_J_PVyError", &B_J_PVyError);
  tree_->Branch("B_J_PVzError", &B_J_PVzError);

  tree_->Branch("B_J_px1", &B_J_px1);
  tree_->Branch("B_J_py1", &B_J_py1);
  tree_->Branch("B_J_pz1", &B_J_pz1);
  tree_->Branch("B_J_pt1", &B_J_pt1);
  tree_->Branch("B_J_eta1", &B_J_eta1);
  tree_->Branch("B_J_phi1", &B_J_phi1);
  tree_->Branch("B_J_charge1", &B_J_charge1);
  tree_->Branch("B_J_soft1", &B_J_soft1);
  tree_->Branch("B_J_tight1", &B_J_tight1);
  tree_->Branch("B_J_loose1", &B_J_loose1);

  tree_->Branch("B_J_px2", &B_J_px2);
  tree_->Branch("B_J_py2", &B_J_py2);
  tree_->Branch("B_J_pz2", &B_J_pz2);
  tree_->Branch("B_J_pt2", &B_J_pt2);
  tree_->Branch("B_J_eta2", &B_J_eta2);
  tree_->Branch("B_J_phi2", &B_J_phi2);
  tree_->Branch("B_J_charge2", &B_J_charge2);
  tree_->Branch("B_J_soft2", &B_J_soft2);
  tree_->Branch("B_J_tight2", &B_J_tight2);
  tree_->Branch("B_J_loose2", &B_J_loose2);
  tree_->Branch("B_J_VtxProb", &B_J_VtxProb);
  tree_->Branch("B_J_xyP", &B_J_xyP);
  tree_->Branch("B_J_xyM", &B_J_xyM);
  tree_->Branch("B_J_zP", &B_J_zP);
  tree_->Branch("B_J_zM", &B_J_zM);
  

 
  tree_->Branch("mumC2",&mumC2);  
  tree_->Branch("mumNHits",&mumNHits);
  tree_->Branch("mumNPHits",&mumNPHits);
  tree_->Branch("mupC2",&mupC2);  
  tree_->Branch("mupNHits",&mupNHits);
  tree_->Branch("mupNPHits",&mupNPHits);
}




// ------------ method called once each job just after ending the event loop  ------------
void miniAODmuons::endJob() {
  tree_->GetDirectory()->cd();
  tree_->Write();
}

//define this as a plug-in
DEFINE_FWK_MODULE(miniAODmuons);
