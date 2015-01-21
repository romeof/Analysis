// -*- C++ -*-
//
// Package:    Analysis/FlatTreer
// Class:      FlatTreer
// 
/**\class FlatTreer FlatTreer.cc Analysis/FlatTreer/plugins/FlatTreer.cc

 Description:
 This class aims to map miniAOD infos into flat tree for stand-alone analyses (using only root)

 Implementation:
 The implementation follows some rules. 
 If you modify part of the code, you are kindly invited to be coherent with the style it is written.
 For instance, pay attention to  
  the way you write comments
  the indentation
  the use of methods and functions in the main code 

 Usefull reference:
 https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
 Variables are stored in a tree described by CTree, implemented in TreeVariables.h
*/
//
// Original Author:  Francesco Romeo
//         Created:  Mon, 19 Jan 2015 08:54:10 GMT
//
//
/////
//   Headers
/////
//System and event
#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TStopwatch.h"
//Gen info
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//Trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
//Vertex
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//Muon
#include "DataFormats/PatCandidates/interface/Muon.h"
//Electron
#include "DataFormats/PatCandidates/interface/Electron.h"
//Tau
#include "DataFormats/PatCandidates/interface/Tau.h"
//Photon
#include "DataFormats/PatCandidates/interface/Photon.h"
//Jet and Met
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
//Packed
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
//Track builder infos
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
//Math
#include "DataFormats/Math/interface/deltaR.h"
//Store info
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
//User defined classes
#include "Analysis/FlatTreer/interface/TreeVariables.h"
#include "Analysis/FlatTreer/interface/ObjEvtFunctions.h"
/////
//   Namespace
/////
using namespace reco;
using namespace edm;
using namespace std;
/////
//   Class declaration
/////
class FlatTreer : public edm::EDAnalyzer {
 public:
 explicit FlatTreer(const edm::ParameterSet&);
 ~FlatTreer();
 static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
 private:
 //Default methods
 virtual void beginJob() override;
 virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
 virtual void endJob() override;
 //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
 //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
 //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
 //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
 //Collections of objects
 edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
 edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;
 edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
 edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
 edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
 edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
 edm::EDGetTokenT<pat::MuonCollection> muonToken_;
 edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
 edm::EDGetTokenT<pat::TauCollection> tauToken_;
 edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
 edm::EDGetTokenT<pat::JetCollection> jetToken_;
 edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
 edm::EDGetTokenT<pat::METCollection> metToken_;
 edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
 edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken_;
 //Values for begin job
 int evt_totnum;
 //Values for the whole analysis
 const double mindr_genobjmatch;
 //Watch time and cpu for the analysis
 TStopwatch* stopwatch;
 //Tree
 CTree *tree;
 const edm::Service<TFileService> fs;
};
/////
//   Constructors and destructor
/////
FlatTreer::FlatTreer(const edm::ParameterSet& iConfig):
 //Collections of objects
 prunedGenToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("pruned"))),
 packedGenToken_(consumes<edm::View<pat::PackedGenParticle> >(iConfig.getParameter<edm::InputTag>("packed"))),
 triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
 triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
 triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
 vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
 muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
 electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
 tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
 photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
 jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
 fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
 metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets"))),
 pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
 lostTracksToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("lostTracks"))),
 //Values for the whole analysis
 mindr_genobjmatch(iConfig.getUntrackedParameter<double>("mindr_genobjmatch")),
 //TTree 
 tree(new CTree(fs->make<TTree>("tree", "tree")))   
{
 //Now do what ever initialization is needed
 tree->make_branches();
 stopwatch = new TStopwatch();
}

FlatTreer::~FlatTreer()
{
 // do anything here that needs to be done at desctruction time
 // (e.g. close files, deallocate resources etc.)
 delete stopwatch;
}
/////
//   Member functions
/////
// ------------ method called for each event  ------------
void FlatTreer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
 using namespace edm;
 evt_totnum++; 
 /////
 //   Handle iEvent.getByToken
 /////
 Handle<edm::View<reco::GenParticle> > pruned;
 iEvent.getByToken(prunedGenToken_,pruned);
 Handle<edm::View<pat::PackedGenParticle> > packed;
 iEvent.getByToken(packedGenToken_,packed);
 edm::Handle<edm::TriggerResults> triggerBits;
 iEvent.getByToken(triggerBits_, triggerBits);
 edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
 iEvent.getByToken(triggerObjects_, triggerObjects);
 edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
 iEvent.getByToken(triggerPrescales_, triggerPrescales);
 edm::Handle<reco::VertexCollection> vertices;
 iEvent.getByToken(vtxToken_, vertices);
 edm::Handle<pat::MuonCollection> muons;
 iEvent.getByToken(muonToken_, muons);
 edm::Handle<pat::ElectronCollection> electrons;
 iEvent.getByToken(electronToken_, electrons);
 edm::Handle<pat::TauCollection> taus;
 iEvent.getByToken(tauToken_, taus);
 edm::Handle<pat::PhotonCollection> photons;
 iEvent.getByToken(photonToken_, photons);
 edm::Handle<pat::JetCollection> jets;
 iEvent.getByToken(jetToken_, jets);
 PFJetIDSelectionFunctor pfLooseJetID(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE),
                         pfTightJetID(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT);
 pat::strbitset passLooseCuts(pfLooseJetID.getBitTemplate()),
                passTightCuts(pfTightJetID.getBitTemplate());
 edm::Handle<pat::JetCollection> fatjets;
 iEvent.getByToken(fatjetToken_, fatjets);
 edm::Handle<pat::METCollection> mets;
 iEvent.getByToken(metToken_, mets);
 edm::Handle<pat::PackedCandidateCollection> pfs;
 iEvent.getByToken(pfToken_, pfs);
 edm::Handle<pat::PackedCandidateCollection> lostrks;
 iEvent.getByToken(lostTracksToken_, lostrks);
 /////
 //   Trigger 
 /////
 //bool istriggeredevent = false;
 //const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
 //for(unsigned int i = 0; i<triggerBits->size(); ++i){
 // if(names.triggerName(i)=="mytriggername" && triggerBits->accept(i)) istriggeredevent = true;
 //} 
 //if(!istriggeredevent) return;
 /////
 //   Primary vertex
 /////
 if(vertices->empty()) return; // skip the event if no PV found
 const reco::Vertex &PV = vertices->front();
 /////
 //   Initial skim
 /////
 int nmuevt = 0;
 for(const pat::Muon &mu : *muons) if(is_tth_mu(mu,PV)) nmuevt++; 
 int neleevt = 0;
 for(const pat::Electron &ele : *electrons) if(is_tth_ele(ele,PV)) neleevt++;
 int njetevt   = 0;
 int nlbjetevt = 0;
 int nmbjetevt = 0;
 for(const pat::Jet &j : *jets){
  if(is_tth_jet(j)){
   njetevt++;
   if(j.bDiscriminator("combinedSecondaryVertexBJetTags")>0.244) nlbjetevt++;
   if(j.bDiscriminator("combinedSecondaryVertexBJetTags")>0.679) nmbjetevt++;
  }
 }
 if((nmuevt+neleevt)<2 || !(njetevt>=4 && (nlbjetevt>=2 || nmbjetevt>=1))) return;
 /////
 //   Take relevant info 
 /////
 tree->loop_initialize();
 tree->evt_id   = (unsigned int)iEvent.id().event(); 
 //tree->evt_json = //Fill when json will be used
 tree->evt_lumi = (unsigned int)iEvent.id().luminosityBlock();
 tree->evt_run  = (unsigned int)iEvent.id().run();
 double prova = 3;
 tree->prova = prova;
 /////
 //   Fill tree
 /////
 tree->tree->Fill();
 /////
 //   End of event
 /////
#ifdef THIS_IS_AN_EVENT_EXAMPLE
 Handle<ExampleData> pIn;
 iEvent.getByLabel("example",pIn);
#endif
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
 ESHandle<SetupData> pSetup;
 iSetup.get<SetupRecord>().get(pSetup);
#endif
}
// ------------ method called once each job just before starting event loop  ------------
void 
FlatTreer::beginJob()
{
 stopwatch->Start();
 evt_totnum = 0;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FlatTreer::endJob() 
{
 stopwatch->Stop(); 
 cout<<endl;
 cout<<"Rapid job summary "<<endl;
 cout<<evt_totnum<<" events analysed in "<<stopwatch->RealTime()<<" seconds"<<endl;
 cout<<endl;
}
// ------------ method called when starting to processes a run  ------------
/*
void 
FlatTreer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
FlatTreer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
FlatTreer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
FlatTreer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FlatTreer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(FlatTreer);
