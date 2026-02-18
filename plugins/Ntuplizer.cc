// -*- C++ -*-
//
// Package:    HLTScoutingAnalyzer/HLTScoutingAnalyzer
// Class:      HLTScoutingAnalyzer
//
/**\class HLTScoutingAnalyzer Ntuplizer.cc HLTScoutingAnalyzer/HLTScoutingAnalyzer/plugins/Ntuplizer.cc

 Description: [one line class summary]
 Designed for reading ntuples with rerun HLT 

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Prijith Pradeep
//         Created:  Wed, 05 Jun 2024 21:53:24 GMT
//
//

// ROOT includes
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TLorentzVector.h"
// Uses regexp for trigger names?
#include "TPRegexp.h"

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"

// dataformats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/fillCovariance.h"

// For vertex distance
#include "DataFormats/Math/interface/Error.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

// L1 trigger dataformats
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/Common/interface/AssociationMap.h"

// PAT dataformats
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// Scouting dataformats
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingHitPatternPOD.h"

// Trigger results
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// Gen info
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

// Track builder
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"

// CMSSW utils
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"

typedef math::XYZPoint Point;
typedef math::Error<3>::type Error3;
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class Ntuplizer : public
edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit Ntuplizer(const edm::ParameterSet&);
  ~Ntuplizer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;
  virtual void endRun(const edm::Run&, const edm::EventSetup&) override;
  virtual void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;
  virtual void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;
  virtual void clearVars();

  // ----------member data ---------------------------
  // Tokens
  const edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> muonsNoVtxToken_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingMuon>> muonsVtxToken_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> PVToken_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> SVNoVtxToken_;
  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex>> SVVtxToken_;
  const edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticleToken_;
  // Extra vertexing 
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbESToken_;

  // HLT 
  std::vector<std::string> hltPaths_;
  std::map<std::string, int> hltPathIndex_;
  std::vector<UChar_t> hltDecisions_;

  // L1
  bool doL1_;
  edm::InputTag algInputTag_;
  edm::InputTag extInputTag_;
  edm::EDGetToken algToken_;
  std::shared_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
  std::vector<std::string> l1Seeds_;
  std::vector<UChar_t> l1Decisions_;

  // Ntuple branches
  UInt_t nScoutingMuonNoVtx;
  std::vector<Float16_t> ScoutingMuonNoVtx_pt;
  std::vector<Float16_t> ScoutingMuonNoVtx_eta;
  std::vector<Float16_t> ScoutingMuonNoVtx_phi;
  std::vector<Float16_t> ScoutingMuonNoVtx_phiCorr;
  std::vector<Float16_t> ScoutingMuonNoVtx_m;
  std::vector<Int_t> ScoutingMuonNoVtx_charge;
  std::vector<Float16_t> ScoutingMuonNoVtx_normalizedChi2;
  std::vector<Float16_t> ScoutingMuonNoVtx_trkchi2;
  std::vector<Float16_t> ScoutingMuonNoVtx_trkndof;
  std::vector<Float16_t> ScoutingMuonNoVtx_trkdxy;
  std::vector<Float16_t> ScoutingMuonNoVtx_trkdz;
  std::vector<Float16_t> ScoutingMuonNoVtx_trkqoverp;
  std::vector<Float16_t> ScoutingMuonNoVtx_trklambda;
  std::vector<Float16_t> ScoutingMuonNoVtx_trkpt;
  std::vector<Float16_t> ScoutingMuonNoVtx_trkphi;
  std::vector<Float16_t> ScoutingMuonNoVtx_trketa;
  std::vector<Float16_t> ScoutingMuonNoVtx_trkqoverpError;
  std::vector<Float16_t> ScoutingMuonNoVtx_trklambdaError;
  std::vector<Float16_t> ScoutingMuonNoVtx_trkdxyError;
  std::vector<Float16_t> ScoutingMuonNoVtx_trkdzError;
  std::vector<Float16_t> ScoutingMuonNoVtx_trkphiError;
  std::vector<Float16_t> ScoutingMuonNoVtx_trkdsz;
  std::vector<Float16_t> ScoutingMuonNoVtx_trkdszError;
  std::vector<Float16_t> ScoutingMuonNoVtx_trkvx;
  std::vector<Float16_t> ScoutingMuonNoVtx_trkvy;
  std::vector<Float16_t> ScoutingMuonNoVtx_trkvz;
  std::vector<std::vector<Int_t>> ScoutingMuonNoVtx_vtxIndx;
  std::vector<Bool_t> ScoutingMuonNoVtx_isGlobal;
  std::vector<Bool_t> ScoutingMuonNoVtx_isTracker;
  std::vector<Bool_t> ScoutingMuonNoVtx_isStandalone;


  UInt_t nScoutingMuonVtx;
  std::vector<Float16_t> ScoutingMuonVtx_pt;
  std::vector<Float16_t> ScoutingMuonVtx_eta;
  std::vector<Float16_t> ScoutingMuonVtx_phi;
  std::vector<Float16_t> ScoutingMuonVtx_phiCorr;
  std::vector<Float16_t> ScoutingMuonVtx_m;
  std::vector<Int_t> ScoutingMuonVtx_charge;
  std::vector<Float16_t> ScoutingMuonVtx_normalizedChi2;
  std::vector<Float16_t> ScoutingMuonVtx_trkchi2;
  std::vector<Float16_t> ScoutingMuonVtx_trkndof;
  std::vector<Float16_t> ScoutingMuonVtx_trkdxy;
  std::vector<Float16_t> ScoutingMuonVtx_trkdz;
  std::vector<Float16_t> ScoutingMuonVtx_trkqoverp;
  std::vector<Float16_t> ScoutingMuonVtx_trklambda;
  std::vector<Float16_t> ScoutingMuonVtx_trkpt;
  std::vector<Float16_t> ScoutingMuonVtx_trkphi;
  std::vector<Float16_t> ScoutingMuonVtx_trketa;
  std::vector<Float16_t> ScoutingMuonVtx_trkqoverpError;
  std::vector<Float16_t> ScoutingMuonVtx_trklambdaError;
  std::vector<Float16_t> ScoutingMuonVtx_trkdxyError;
  std::vector<Float16_t> ScoutingMuonVtx_trkdzError;
  std::vector<Float16_t> ScoutingMuonVtx_trkphiError;
  std::vector<Float16_t> ScoutingMuonVtx_trkdsz;
  std::vector<Float16_t> ScoutingMuonVtx_trkdszError;
  std::vector<Float16_t> ScoutingMuonVtx_trkvx;
  std::vector<Float16_t> ScoutingMuonVtx_trkvy;
  std::vector<Float16_t> ScoutingMuonVtx_trkvz;
  std::vector<std::vector<Int_t>> ScoutingMuonVtx_vtxIndx;
  std::vector<Bool_t> ScoutingMuonVtx_isGlobal;
  std::vector<Bool_t> ScoutingMuonVtx_isTracker;
  std::vector<Bool_t> ScoutingMuonVtx_isStandalone;

  // Vertices
  UInt_t nPV;
  Float16_t PV_x;
  Float16_t PV_y;
  Float16_t PV_z;
  Float16_t PV_xError;
  Float16_t PV_yError;
  Float16_t PV_zError;
  Int_t PV_trksize;
  Float16_t PV_chi2;
  Float16_t PV_ndof;
  Bool_t PV_isvalidvtx;

  UInt_t nSVNoVtx;
  std::vector<Float16_t> SVNoVtx_x;
  std::vector<Float16_t> SVNoVtx_y;
  std::vector<Float16_t> SVNoVtx_z;
  std::vector<Float16_t> SVNoVtx_xError;
  std::vector<Float16_t> SVNoVtx_yError;
  std::vector<Float16_t> SVNoVtx_zError;
  std::vector<Int_t> SVNoVtx_trksize;
  std::vector<Float16_t> SVNoVtx_chi2;
  std::vector<Float16_t> SVNoVtx_ndof;
  std::vector<Bool_t> SVNoVtx_isvalidvtx;
  // Calculated
  std::vector<Float16_t> SVNoVtx_dxy;
  std::vector<Float16_t> SVNoVtx_dxySig;
  std::vector<Float16_t> SVNoVtx_dlen;
  std::vector<Float16_t> SVNoVtx_dlenSig;
  std::vector<Float16_t> SVNoVtx_mass;
  std::vector<Int_t> SVNoVtx_nMuon;

  UInt_t nSVVtx;
  std::vector<Float16_t> SVVtx_x;
  std::vector<Float16_t> SVVtx_y;
  std::vector<Float16_t> SVVtx_z;
  std::vector<Float16_t> SVVtx_xError;
  std::vector<Float16_t> SVVtx_yError;
  std::vector<Float16_t> SVVtx_zError;
  std::vector<Int_t> SVVtx_trksize;
  std::vector<Float16_t> SVVtx_chi2;
  std::vector<Float16_t> SVVtx_ndof;
  std::vector<Bool_t> SVVtx_isvalidvtx;
  // Calculated
  std::vector<Float16_t> SVVtx_dxy;
  std::vector<Float16_t> SVVtx_dxySig;
  std::vector<Float16_t> SVVtx_dlen;
  std::vector<Float16_t> SVVtx_dlenSig;
  std::vector<Float16_t> SVVtx_mass;
  std::vector<Int_t> SVVtx_nMuon;

  // GenPart
  UInt_t nGenPart;
  std::vector<Float16_t> GenPart_pt;
  std::vector<Float16_t> GenPart_eta;
  std::vector<Float16_t> GenPart_phi;
  std::vector<Float16_t> GenPart_m;
  std::vector<Int_t> GenPart_pdgId;
  std::vector<Int_t> GenPart_status;
  std::vector<Int_t> GenPart_genPartIdxMother;

  TTree* Events;

  Int_t run;
  Int_t event;
  Int_t lumi;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig):
  triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
  muonsNoVtxToken_(consumes<std::vector<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("muonsNoVtx"))),
  muonsVtxToken_(consumes<std::vector<Run3ScoutingMuon>>(iConfig.getParameter<edm::InputTag>("muonsVtx"))),
  PVToken_(consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("PV"))),
  SVNoVtxToken_(consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("SVNoVtx"))),
  SVVtxToken_(consumes<std::vector<Run3ScoutingVertex>>(iConfig.getParameter<edm::InputTag>("SVVtx"))),
  genParticleToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"))),
  ttbESToken_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
  hltPaths_(iConfig.getParameter<std::vector<std::string>>("hltPaths")),
  doL1_(iConfig.getParameter<bool>("doL1")),
  l1Seeds_(iConfig.getParameter<std::vector<std::string>>("l1Seeds")){
  
  // Necessary for L1 seeds
  algToken_ = consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("AlgInputTag"));
  l1GtUtils_ = std::make_shared<l1t::L1TGlobalUtil>(iConfig, consumesCollector(), l1t::UseEventSetupIn::RunAndEvent);

  usesResource("TFileService");
  edm::Service<TFileService> fs;
  Events = fs->make<TTree>("Events", "Events");

  Events->Branch("run", &run, "run/I");
  Events->Branch("event", &event, "event/I");
  Events->Branch("lumi", &lumi, "lumi/I");

  // Trigger paths
  // HLT
  hltDecisions_.assign(hltPaths_.size(), 0);
  for (unsigned i = 0; i < hltPaths_.size(); ++i) {
    std::string br = hltPaths_[i];
    Events->Branch(br.c_str(), &hltDecisions_[i], (br + "/b").c_str());
  }

  if (doL1_) {
    l1Decisions_.assign(l1Seeds_.size(), 0);
    for (unsigned i = 0; i < l1Seeds_.size(); ++i) {
      std::string br = l1Seeds_[i];
      Events->Branch(br.c_str(), &l1Decisions_[i], (br + "/b").c_str());
    }
  }

  // Collections
  Events->Branch("nScoutingMuonNoVtx", &nScoutingMuonNoVtx, "nScoutingMuonNoVtx/i");
  Events->Branch("ScoutingMuonNoVtx_pt", &ScoutingMuonNoVtx_pt);
  Events->Branch("ScoutingMuonNoVtx_eta", &ScoutingMuonNoVtx_eta);
  Events->Branch("ScoutingMuonNoVtx_phi", &ScoutingMuonNoVtx_phi);
  Events->Branch("ScoutingMuonNoVtx_phiCorr", &ScoutingMuonNoVtx_phiCorr);
  Events->Branch("ScoutingMuonNoVtx_m", &ScoutingMuonNoVtx_m);
  Events->Branch("ScoutingMuonNoVtx_charge", &ScoutingMuonNoVtx_charge);
  Events->Branch("ScoutingMuonNoVtx_normalizedChi2", &ScoutingMuonNoVtx_normalizedChi2);
  Events->Branch("ScoutingMuonNoVtx_trkchi2", &ScoutingMuonNoVtx_trkchi2);
  Events->Branch("ScoutingMuonNoVtx_trkndof", &ScoutingMuonNoVtx_trkndof);
  Events->Branch("ScoutingMuonNoVtx_trkdxy", &ScoutingMuonNoVtx_trkdxy);
  Events->Branch("ScoutingMuonNoVtx_trkdz", &ScoutingMuonNoVtx_trkdz);
  Events->Branch("ScoutingMuonNoVtx_trkqoverp", &ScoutingMuonNoVtx_trkqoverp);
  Events->Branch("ScoutingMuonNoVtx_trklambda", &ScoutingMuonNoVtx_trklambda);
  Events->Branch("ScoutingMuonNoVtx_trkpt", &ScoutingMuonNoVtx_trkpt);
  Events->Branch("ScoutingMuonNoVtx_trkphi", &ScoutingMuonNoVtx_trkphi);
  Events->Branch("ScoutingMuonNoVtx_trketa", &ScoutingMuonNoVtx_trketa);
  Events->Branch("ScoutingMuonNoVtx_trkqoverpError", &ScoutingMuonNoVtx_trkqoverpError);
  Events->Branch("ScoutingMuonNoVtx_trklambdaError", &ScoutingMuonNoVtx_trklambdaError);
  Events->Branch("ScoutingMuonNoVtx_trkdxyError", &ScoutingMuonNoVtx_trkdxyError);
  Events->Branch("ScoutingMuonNoVtx_trkdzError", &ScoutingMuonNoVtx_trkdzError);
  Events->Branch("ScoutingMuonNoVtx_trkphiError", &ScoutingMuonNoVtx_trkphiError);
  Events->Branch("ScoutingMuonNoVtx_trkdsz", &ScoutingMuonNoVtx_trkdsz);
  Events->Branch("ScoutingMuonNoVtx_trkdszError", &ScoutingMuonNoVtx_trkdszError);
  Events->Branch("ScoutingMuonNoVtx_trkvx", &ScoutingMuonNoVtx_trkvx);
  Events->Branch("ScoutingMuonNoVtx_trkvy", &ScoutingMuonNoVtx_trkvy);
  Events->Branch("ScoutingMuonNoVtx_trkvz", &ScoutingMuonNoVtx_trkvz);
  Events->Branch("ScoutingMuonNoVtx_vtxIndx", &ScoutingMuonNoVtx_vtxIndx);
  Events->Branch("ScoutingMuonNoVtx_isGlobal", &ScoutingMuonNoVtx_isGlobal);
  Events->Branch("ScoutingMuonNoVtx_isTracker", &ScoutingMuonNoVtx_isTracker);
  Events->Branch("ScoutingMuonNoVtx_isStandalone", &ScoutingMuonNoVtx_isStandalone);

  Events->Branch("nScoutingMuonVtx", &nScoutingMuonVtx, "nScoutingMuonVtx/i");
  Events->Branch("ScoutingMuonVtx_pt", &ScoutingMuonVtx_pt);
  Events->Branch("ScoutingMuonVtx_eta", &ScoutingMuonVtx_eta);
  Events->Branch("ScoutingMuonVtx_phi", &ScoutingMuonVtx_phi);
  Events->Branch("ScoutingMuonVtx_phiCorr", &ScoutingMuonVtx_phiCorr);
  Events->Branch("ScoutingMuonVtx_m", &ScoutingMuonVtx_m);
  Events->Branch("ScoutingMuonVtx_charge", &ScoutingMuonVtx_charge);
  Events->Branch("ScoutingMuonVtx_normalizedChi2", &ScoutingMuonVtx_normalizedChi2);
  Events->Branch("ScoutingMuonVtx_trkchi2", &ScoutingMuonVtx_trkchi2);
  Events->Branch("ScoutingMuonVtx_trkndof", &ScoutingMuonVtx_trkndof);
  Events->Branch("ScoutingMuonVtx_trkdxy", &ScoutingMuonVtx_trkdxy);
  Events->Branch("ScoutingMuonVtx_trkdz", &ScoutingMuonVtx_trkdz);
  Events->Branch("ScoutingMuonVtx_trkqoverp", &ScoutingMuonVtx_trkqoverp);
  Events->Branch("ScoutingMuonVtx_trklambda", &ScoutingMuonVtx_trklambda);
  Events->Branch("ScoutingMuonVtx_trkpt", &ScoutingMuonVtx_trkpt);
  Events->Branch("ScoutingMuonVtx_trkphi", &ScoutingMuonVtx_trkphi);
  Events->Branch("ScoutingMuonVtx_trketa", &ScoutingMuonVtx_trketa);
  Events->Branch("ScoutingMuonVtx_trkqoverpError", &ScoutingMuonVtx_trkqoverpError);
  Events->Branch("ScoutingMuonVtx_trklambdaError", &ScoutingMuonVtx_trklambdaError);
  Events->Branch("ScoutingMuonVtx_trkdxyError", &ScoutingMuonVtx_trkdxyError);
  Events->Branch("ScoutingMuonVtx_trkdzError", &ScoutingMuonVtx_trkdzError);
  Events->Branch("ScoutingMuonVtx_trkphiError", &ScoutingMuonVtx_trkphiError);
  Events->Branch("ScoutingMuonVtx_trkdsz", &ScoutingMuonVtx_trkdsz);
  Events->Branch("ScoutingMuonVtx_trkdszError", &ScoutingMuonVtx_trkdszError);
  Events->Branch("ScoutingMuonVtx_trkvx", &ScoutingMuonVtx_trkvx);
  Events->Branch("ScoutingMuonVtx_trkvy", &ScoutingMuonVtx_trkvy);
  Events->Branch("ScoutingMuonVtx_trkvz", &ScoutingMuonVtx_trkvz);
  Events->Branch("ScoutingMuonVtx_vtxIndx", &ScoutingMuonVtx_vtxIndx);
  Events->Branch("ScoutingMuonVtx_isGlobal", &ScoutingMuonVtx_isGlobal);
  Events->Branch("ScoutingMuonVtx_isTracker", &ScoutingMuonVtx_isTracker);
  Events->Branch("ScoutingMuonVtx_isStandalone", &ScoutingMuonVtx_isStandalone);

  // Vertices
  Events->Branch("nPV", &nPV, "nPV/i");
  Events->Branch("PV_x", &PV_x, "PV_x/F");
  Events->Branch("PV_y", &PV_y, "PV_y/F");
  Events->Branch("PV_z", &PV_z, "PV_z/F");
  Events->Branch("PV_xError", &PV_xError, "PV_xError/F");
  Events->Branch("PV_yError", &PV_yError, "PV_yError/F");
  Events->Branch("PV_zError", &PV_zError, "PV_zError/F");
  Events->Branch("PV_trksize", &PV_trksize, "PV_trksize/I");
  Events->Branch("PV_chi2", &PV_chi2, "PV_chi2/F");
  Events->Branch("PV_ndof", &PV_ndof, "PV_ndof/F");
  Events->Branch("PV_isvalidvtx", &PV_isvalidvtx);

  Events->Branch("nSVNoVtx", &nSVNoVtx, "nSVNoVtx/i");
  Events->Branch("SVNoVtx_x", &SVNoVtx_x);
  Events->Branch("SVNoVtx_y", &SVNoVtx_y);
  Events->Branch("SVNoVtx_z", &SVNoVtx_z);
  Events->Branch("SVNoVtx_xError", &SVNoVtx_xError);
  Events->Branch("SVNoVtx_yError", &SVNoVtx_yError);
  Events->Branch("SVNoVtx_zError", &SVNoVtx_zError);
  Events->Branch("SVNoVtx_trksize", &SVNoVtx_trksize);
  Events->Branch("SVNoVtx_chi2", &SVNoVtx_chi2);
  Events->Branch("SVNoVtx_ndof", &SVNoVtx_ndof);
  Events->Branch("SVNoVtx_isvalidvtx", &SVNoVtx_isvalidvtx);
  Events->Branch("SVNoVtx_dxy", &SVNoVtx_dxy);
  Events->Branch("SVNoVtx_dxySig", &SVNoVtx_dxySig);
  Events->Branch("SVNoVtx_dlen", &SVNoVtx_dlen);
  Events->Branch("SVNoVtx_dlenSig", &SVNoVtx_dlenSig);
  Events->Branch("SVNoVtx_mass", &SVNoVtx_mass);
  Events->Branch("SVNoVtx_nMuon", &SVNoVtx_nMuon);

  Events->Branch("nSVVtx", &nSVVtx, "nSVVtx/i");
  Events->Branch("SVVtx_x", &SVVtx_x);
  Events->Branch("SVVtx_y", &SVVtx_y);
  Events->Branch("SVVtx_z", &SVVtx_z);
  Events->Branch("SVVtx_xError", &SVVtx_xError);
  Events->Branch("SVVtx_yError", &SVVtx_yError);
  Events->Branch("SVVtx_zError", &SVVtx_zError);
  Events->Branch("SVVtx_trksize", &SVVtx_trksize);
  Events->Branch("SVVtx_chi2", &SVVtx_chi2);
  Events->Branch("SVVtx_ndof", &SVVtx_ndof);
  Events->Branch("SVVtx_isvalidvtx", &SVVtx_isvalidvtx);
  Events->Branch("SVVtx_dxy", &SVVtx_dxy);
  Events->Branch("SVVtx_dxySig", &SVVtx_dxySig);
  Events->Branch("SVVtx_dlen", &SVVtx_dlen);
  Events->Branch("SVVtx_dlenSig", &SVVtx_dlenSig);
  Events->Branch("SVVtx_mass", &SVVtx_mass);
  Events->Branch("SVVtx_nMuon", &SVVtx_nMuon);

  Events->Branch("nGenPart", &nGenPart);
  Events->Branch("GenPart_pt", &GenPart_pt);
  Events->Branch("GenPart_eta", &GenPart_eta);
  Events->Branch("GenPart_phi", &GenPart_phi);
  Events->Branch("GenPart_m", &GenPart_m);
  Events->Branch("GenPart_pdgId", &GenPart_pdgId);
  Events->Branch("GenPart_status", &GenPart_status);
  Events->Branch("GenPart_genPartIdxMother", &GenPart_genPartIdxMother);
};

Ntuplizer::~Ntuplizer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByToken(triggerResultsToken_, triggerResults);
  
  edm::Handle<std::vector<Run3ScoutingMuon>> muonsNoVtx;
  iEvent.getByToken(muonsNoVtxToken_, muonsNoVtx);

  edm::Handle<std::vector<Run3ScoutingMuon>> muonsVtx;
  iEvent.getByToken(muonsVtxToken_, muonsVtx);

  edm::Handle<std::vector<Run3ScoutingVertex>> PV;
  iEvent.getByToken(PVToken_, PV);

  edm::Handle<std::vector<Run3ScoutingVertex>> SVNoVtx;
  iEvent.getByToken(SVNoVtxToken_, SVNoVtx);

  edm::Handle<std::vector<Run3ScoutingVertex>> SVVtx;
  iEvent.getByToken(SVVtxToken_, SVVtx);

  edm::Handle<std::vector<reco::GenParticle>> genParticles;
  iEvent.getByToken(genParticleToken_, genParticles);

  const TransientTrackBuilder* theB = &iSetup.getData(ttbESToken_);

  run = iEvent.eventAuxiliary().run();
  event = iEvent.eventAuxiliary().event();
  lumi = iEvent.eventAuxiliary().luminosityBlock();

  std::fill(hltDecisions_.begin(), hltDecisions_.end(), 0);
  if (doL1_) std::fill(l1Decisions_.begin(), l1Decisions_.end(), 0);

  // Read out trigger decisions
  // Fill with L1 decision
  if (doL1_){
    l1GtUtils_->retrieveL1(iEvent,iSetup,algToken_);
    for (unsigned i = 0; i < l1Seeds_.size(); ++i) {
      bool pass = false;
      // If seed name not found, pass stays false
      l1GtUtils_->getFinalDecisionByName(l1Seeds_[i], pass);
      l1Decisions_[i] = pass ? 1 : 0;
    }
  }

  if (triggerResults.isValid()) {
    for (unsigned i = 0; i < hltPaths_.size(); ++i) {
      auto it = hltPathIndex_.find(hltPaths_[i]);
      if (it == hltPathIndex_.end()) continue;
      int idx = it->second;
      if (idx < 0) continue;
      hltDecisions_[i] = triggerResults->accept(idx) ? 1 : 0;
    }
  }

  // Fill collections
  nScoutingMuonNoVtx = 0;
  if(muonsNoVtx.isValid()){
    for(auto &muonNoVtx: *muonsNoVtx){
      ScoutingMuonNoVtx_pt.push_back(muonNoVtx.pt());
      ScoutingMuonNoVtx_eta.push_back(muonNoVtx.eta());
      ScoutingMuonNoVtx_phi.push_back(muonNoVtx.phi());
      ScoutingMuonNoVtx_m.push_back(muonNoVtx.m());
      ScoutingMuonNoVtx_charge.push_back(muonNoVtx.charge());
      ScoutingMuonNoVtx_normalizedChi2.push_back(muonNoVtx.normalizedChi2());
      ScoutingMuonNoVtx_trkchi2.push_back(muonNoVtx.trk_chi2());
      ScoutingMuonNoVtx_trkndof.push_back(muonNoVtx.trk_ndof());
      ScoutingMuonNoVtx_trkdxy.push_back(muonNoVtx.trk_dxy());
      ScoutingMuonNoVtx_trkdz.push_back(muonNoVtx.trk_dz());
      ScoutingMuonNoVtx_trkqoverp.push_back(muonNoVtx.trk_qoverp());
      ScoutingMuonNoVtx_trklambda.push_back(muonNoVtx.trk_lambda());
      ScoutingMuonNoVtx_trkpt.push_back(muonNoVtx.trk_pt());
      ScoutingMuonNoVtx_trkphi.push_back(muonNoVtx.trk_phi());
      ScoutingMuonNoVtx_trketa.push_back(muonNoVtx.trk_eta());
      ScoutingMuonNoVtx_trkqoverpError.push_back(muonNoVtx.trk_qoverpError());
      ScoutingMuonNoVtx_trklambdaError.push_back(muonNoVtx.trk_lambdaError());
      ScoutingMuonNoVtx_trkdxyError.push_back(muonNoVtx.trk_dxyError());
      ScoutingMuonNoVtx_trkdzError.push_back(muonNoVtx.trk_dzError());
      ScoutingMuonNoVtx_trkphiError.push_back(muonNoVtx.trk_phiError());
      ScoutingMuonNoVtx_trkdsz.push_back(muonNoVtx.trk_dsz());
      ScoutingMuonNoVtx_trkdszError.push_back(muonNoVtx.trk_dszError());
      ScoutingMuonNoVtx_trkvx.push_back(muonNoVtx.trk_vx());
      ScoutingMuonNoVtx_trkvy.push_back(muonNoVtx.trk_vy());
      ScoutingMuonNoVtx_trkvz.push_back(muonNoVtx.trk_vz());
      ScoutingMuonNoVtx_vtxIndx.push_back(muonNoVtx.vtxIndx());
      ScoutingMuonNoVtx_isGlobal.push_back(muonNoVtx.isGlobalMuon());
      ScoutingMuonNoVtx_isTracker.push_back(muonNoVtx.isTrackerMuon());
      ScoutingMuonNoVtx_isStandalone.push_back(muonNoVtx.type() & 1 << 3);

      // Perform track extrapolation to the first SV in vtxIndx has nonzero size
      float muonNoVtx_phiCorr = muonNoVtx.phi(); 
      //Build the muon track and then make it transient
      reco::Track::Point v(muonNoVtx.trk_vx(), muonNoVtx.trk_vy(), muonNoVtx.trk_vz());
      reco::Track::Vector p(
        muonNoVtx.trk_pt()*std::cos(muonNoVtx.trk_phi()),
        muonNoVtx.trk_pt()*std::sin(muonNoVtx.trk_phi()), 
        muonNoVtx.trk_pt()*std::sinh(muonNoVtx.trk_eta())
      );

      double vec[15];
      for (auto i = 0; i < 15; i++) vec[i] = 1.;
      reco::TrackBase::CovarianceMatrix cov(vec, vec + 15);
      cov(0, 0) = std::pow(muonNoVtx.trk_qoverpError(),2);
      cov(0, 1) = muonNoVtx.trk_qoverp_lambda_cov();
      cov(0, 2) = muonNoVtx.trk_qoverp_phi_cov();
      cov(0, 3) = muonNoVtx.trk_qoverp_dxy_cov();
      cov(0, 4) = muonNoVtx.trk_qoverp_dsz_cov();
      cov(1, 1) = std::pow(muonNoVtx.trk_lambdaError(),2);
      cov(1, 2) = muonNoVtx.trk_lambda_phi_cov();
      cov(1, 3) = muonNoVtx.trk_lambda_dxy_cov();
      cov(1, 4) = muonNoVtx.trk_lambda_dsz_cov();
      cov(2, 2) = std::pow(muonNoVtx.trk_phiError(),2);
      cov(2, 3) = muonNoVtx.trk_phi_dxy_cov();
      cov(2, 4) = muonNoVtx.trk_phi_dsz_cov();
      cov(3, 3) = std::pow(muonNoVtx.trk_dxyError(),2);
      cov(3, 4) = muonNoVtx.trk_dxy_dsz_cov();
      cov(4, 4) = std::pow(muonNoVtx.trk_dszError(),2);
      reco::Track trk(muonNoVtx.trk_chi2(), muonNoVtx.trk_ndof(), v, p, muonNoVtx.charge(), cov);

      reco::TransientTrack trans = theB->build(trk);

      if(muonNoVtx.vtxIndx().size()>0){
        auto vtxIndx = (muonNoVtx.vtxIndx())[0];
        auto svNoVtx = (*SVNoVtx)[vtxIndx];
        ROOT::Math::PtEtaPhiMVector muonNoVtx_p4(0, 0, 0, 0);
        GlobalPoint svNoVtx_pos(svNoVtx.x(), svNoVtx.y(), svNoVtx.z());
        TrajectoryStateClosestToPoint traj = trans.trajectoryStateClosestToPoint(svNoVtx_pos);
        GlobalVector muonNoVtx_p3prop = traj.momentum();
        muonNoVtx_phiCorr = muonNoVtx_p3prop.phi();
      }

      ScoutingMuonNoVtx_phiCorr.push_back(muonNoVtx_phiCorr);
      nScoutingMuonNoVtx++;
    }  
  }

  nScoutingMuonVtx = 0;
  if(muonsVtx.isValid()){
    for(auto &muonVtx: *muonsVtx){
      ScoutingMuonVtx_pt.push_back(muonVtx.pt());
      ScoutingMuonVtx_eta.push_back(muonVtx.eta());
      ScoutingMuonVtx_phi.push_back(muonVtx.phi());
      ScoutingMuonVtx_m.push_back(muonVtx.m());
      ScoutingMuonVtx_charge.push_back(muonVtx.charge());
      ScoutingMuonVtx_normalizedChi2.push_back(muonVtx.normalizedChi2());
      ScoutingMuonVtx_trkchi2.push_back(muonVtx.trk_chi2());
      ScoutingMuonVtx_trkndof.push_back(muonVtx.trk_ndof());
      ScoutingMuonVtx_trkdxy.push_back(muonVtx.trk_dxy());
      ScoutingMuonVtx_trkdz.push_back(muonVtx.trk_dz());
      ScoutingMuonVtx_trkqoverp.push_back(muonVtx.trk_qoverp());
      ScoutingMuonVtx_trklambda.push_back(muonVtx.trk_lambda());
      ScoutingMuonVtx_trkpt.push_back(muonVtx.trk_pt());
      ScoutingMuonVtx_trkphi.push_back(muonVtx.trk_phi());
      ScoutingMuonVtx_trketa.push_back(muonVtx.trk_eta());
      ScoutingMuonVtx_trkqoverpError.push_back(muonVtx.trk_qoverpError());
      ScoutingMuonVtx_trklambdaError.push_back(muonVtx.trk_lambdaError());
      ScoutingMuonVtx_trkdxyError.push_back(muonVtx.trk_dxyError());
      ScoutingMuonVtx_trkdzError.push_back(muonVtx.trk_dzError());
      ScoutingMuonVtx_trkphiError.push_back(muonVtx.trk_phiError());
      ScoutingMuonVtx_trkdsz.push_back(muonVtx.trk_dsz());
      ScoutingMuonVtx_trkdszError.push_back(muonVtx.trk_dszError());
      ScoutingMuonVtx_trkvx.push_back(muonVtx.trk_vx());
      ScoutingMuonVtx_trkvy.push_back(muonVtx.trk_vy());
      ScoutingMuonVtx_trkvz.push_back(muonVtx.trk_vz());
      ScoutingMuonVtx_vtxIndx.push_back(muonVtx.vtxIndx());
      ScoutingMuonVtx_isGlobal.push_back(muonVtx.isGlobalMuon());
      ScoutingMuonVtx_isTracker.push_back(muonVtx.isTrackerMuon());
      ScoutingMuonVtx_isStandalone.push_back(muonVtx.type() & 1 << 3);

      // Perform track extrapolation to the first SV in vtxIndx has nonzero size
      float muonVtx_phiCorr = muonVtx.phi(); 
      //Build the muon track and then make it transient
      reco::Track::Point v(muonVtx.trk_vx(), muonVtx.trk_vy(), muonVtx.trk_vz());
      reco::Track::Vector p(
        muonVtx.trk_pt()*std::cos(muonVtx.trk_phi()),
        muonVtx.trk_pt()*std::sin(muonVtx.trk_phi()), 
        muonVtx.trk_pt()*std::sinh(muonVtx.trk_eta())
      );

      double vec[15];
      for (auto i = 0; i < 15; i++) vec[i] = 1.;
      reco::TrackBase::CovarianceMatrix cov(vec, vec + 15);
      cov(0, 0) = std::pow(muonVtx.trk_qoverpError(),2);
      cov(0, 1) = muonVtx.trk_qoverp_lambda_cov();
      cov(0, 2) = muonVtx.trk_qoverp_phi_cov();
      cov(0, 3) = muonVtx.trk_qoverp_dxy_cov();
      cov(0, 4) = muonVtx.trk_qoverp_dsz_cov();
      cov(1, 1) = std::pow(muonVtx.trk_lambdaError(),2);
      cov(1, 2) = muonVtx.trk_lambda_phi_cov();
      cov(1, 3) = muonVtx.trk_lambda_dxy_cov();
      cov(1, 4) = muonVtx.trk_lambda_dsz_cov();
      cov(2, 2) = std::pow(muonVtx.trk_phiError(),2);
      cov(2, 3) = muonVtx.trk_phi_dxy_cov();
      cov(2, 4) = muonVtx.trk_phi_dsz_cov();
      cov(3, 3) = std::pow(muonVtx.trk_dxyError(),2);
      cov(3, 4) = muonVtx.trk_dxy_dsz_cov();
      cov(4, 4) = std::pow(muonVtx.trk_dszError(),2);
      reco::Track trk(muonVtx.trk_chi2(), muonVtx.trk_ndof(), v, p, muonVtx.charge(), cov);

      reco::TransientTrack trans = theB->build(trk);

      if(muonVtx.vtxIndx().size()>0){
        auto vtxIndx = (muonVtx.vtxIndx())[0];
        auto svVtx = (*SVVtx)[vtxIndx];
        ROOT::Math::PtEtaPhiMVector muonVtx_p4(0, 0, 0, 0);
        GlobalPoint svVtx_pos(svVtx.x(), svVtx.y(), svVtx.z());
        TrajectoryStateClosestToPoint traj = trans.trajectoryStateClosestToPoint(svVtx_pos);
        GlobalVector muonVtx_p3prop = traj.momentum();
        muonVtx_phiCorr = muonVtx_p3prop.phi();
      }

      ScoutingMuonVtx_phiCorr.push_back(muonVtx_phiCorr);
      nScoutingMuonVtx++;
    }  
  }

  // Keep a reco vertex object of the first PV
  nPV = 0;
  auto PV0Ptr = PV->begin();
  Point PV0_pos(PV0Ptr->x(), PV0Ptr->y(), PV0Ptr->z());
  Error3 PV0_err;
  PV0_err(0, 0) = std::pow(PV0Ptr->xError(), 2);
  PV0_err(1, 1) = std::pow(PV0Ptr->yError(), 2);
  PV0_err(2, 2) = std::pow(PV0Ptr->zError(), 2);
  reco::Vertex PV0(PV0_pos, PV0_err, PV0Ptr->chi2(), PV0Ptr->ndof(), PV0Ptr->tracksSize());

  //Fill branches
  PV_x = PV0Ptr->x();
  PV_y = PV0Ptr->y();
  PV_z = PV0Ptr->z();
  PV_xError = PV0Ptr->xError();
  PV_yError = PV0Ptr->yError();
  PV_zError = PV0Ptr->zError();
  PV_trksize = PV0Ptr->tracksSize();
  PV_chi2 = PV0Ptr->chi2();
  PV_ndof = PV0Ptr->ndof();
  PV_isvalidvtx = PV0Ptr->isValidVtx();

  for(auto &PViter: *PV){
    nPV++;
  }

  nSVNoVtx = 0;
  VertexDistance3D vdist;
  VertexDistanceXY vdistXY;

  if(SVNoVtx.isValid()){
    for(auto &svNoVtx: *SVNoVtx){
      SVNoVtx_x.push_back(svNoVtx.x());
      SVNoVtx_y.push_back(svNoVtx.y());
      SVNoVtx_z.push_back(svNoVtx.z());
      SVNoVtx_xError.push_back(svNoVtx.xError());
      SVNoVtx_yError.push_back(svNoVtx.yError());
      SVNoVtx_zError.push_back(svNoVtx.zError());
      SVNoVtx_trksize.push_back(svNoVtx.tracksSize());
      SVNoVtx_chi2.push_back(svNoVtx.chi2());
      SVNoVtx_ndof.push_back(svNoVtx.ndof());
      SVNoVtx_isvalidvtx.push_back(svNoVtx.isValidVtx());
      
      //Calculated for PV0
      Point SVNoVtx_pos(svNoVtx.x(), svNoVtx.y(), svNoVtx.z());
      Error3 SVNoVtx_err;
      SVNoVtx_err(0, 0) = std::pow(svNoVtx.xError(), 2);
      SVNoVtx_err(1, 1) = std::pow(svNoVtx.yError(), 2);
      SVNoVtx_err(2, 2) = std::pow(svNoVtx.zError(), 2);
      reco::Vertex SVNoVtx_cand(SVNoVtx_pos, SVNoVtx_err, svNoVtx.chi2(), svNoVtx.ndof(), svNoVtx.tracksSize());
      Measurement1D dxy = vdistXY.distance(PV0, SVNoVtx_cand);
      Measurement1D dlen = vdist.distance(PV0, SVNoVtx_cand);
      SVNoVtx_dxy.push_back(dxy.value());
      SVNoVtx_dxySig.push_back(dxy.significance());
      SVNoVtx_dlen.push_back(dlen.value());
      SVNoVtx_dlenSig.push_back(dlen.significance());

      //Iterate over muons to find SV match
      int nMuonMatch = 0;
      float svNoVtx_mass = -1.;
      TLorentzVector svNoVtx_p4;
      for(UInt_t i=0; i<nScoutingMuonNoVtx; i++){
        std::vector<int> scoutingmuonNoVtx_vtxIndx = ScoutingMuonNoVtx_vtxIndx[i];
        if(std::find(scoutingmuonNoVtx_vtxIndx.begin(), scoutingmuonNoVtx_vtxIndx.end(), nSVNoVtx) != scoutingmuonNoVtx_vtxIndx.end()){
          //Muon matched
          nMuonMatch++;
          TLorentzVector scoutingmuonNoVtx_p4;
          scoutingmuonNoVtx_p4.SetPtEtaPhiM(ScoutingMuonNoVtx_pt[i], ScoutingMuonNoVtx_eta[i], ScoutingMuonNoVtx_phiCorr[i], 0.10566);
          svNoVtx_p4 += scoutingmuonNoVtx_p4;
        }
      }
      if(nMuonMatch>0) svNoVtx_mass = svNoVtx_p4.M();

      SVNoVtx_mass.push_back(svNoVtx_mass);
      SVNoVtx_nMuon.push_back(nMuonMatch);

      nSVNoVtx++;
    }
  }

  nSVVtx = 0;
  if(SVVtx.isValid()){
    for(auto &svVtx: *SVVtx){
      SVVtx_x.push_back(svVtx.x());
      SVVtx_y.push_back(svVtx.y());
      SVVtx_z.push_back(svVtx.z());
      SVVtx_xError.push_back(svVtx.xError());
      SVVtx_yError.push_back(svVtx.yError());
      SVVtx_zError.push_back(svVtx.zError());
      SVVtx_trksize.push_back(svVtx.tracksSize());
      SVVtx_chi2.push_back(svVtx.chi2());
      SVVtx_ndof.push_back(svVtx.ndof());
      SVVtx_isvalidvtx.push_back(svVtx.isValidVtx());
      
      //Calculated for PV0
      Point SVVtx_pos(svVtx.x(), svVtx.y(), svVtx.z());
      Error3 SVVtx_err;
      SVVtx_err(0, 0) = std::pow(svVtx.xError(), 2);
      SVVtx_err(1, 1) = std::pow(svVtx.yError(), 2);
      SVVtx_err(2, 2) = std::pow(svVtx.zError(), 2);
      reco::Vertex SVVtx_cand(SVVtx_pos, SVVtx_err, svVtx.chi2(), svVtx.ndof(), svVtx.tracksSize());
      Measurement1D dxy = vdistXY.distance(PV0, SVVtx_cand);
      Measurement1D dlen = vdist.distance(PV0, SVVtx_cand);
      SVVtx_dxy.push_back(dxy.value());
      SVVtx_dxySig.push_back(dxy.significance());
      SVVtx_dlen.push_back(dlen.value());
      SVVtx_dlenSig.push_back(dlen.significance());

      //Iterate over muons to find SV match
      int nMuonMatch = 0;
      float svVtx_mass = -1.;
      TLorentzVector svVtx_p4;
      for(UInt_t i=0; i<nScoutingMuonVtx; i++){
        std::vector<int> scoutingmuonVtx_vtxIndx = ScoutingMuonVtx_vtxIndx[i];
        if(std::find(scoutingmuonVtx_vtxIndx.begin(), scoutingmuonVtx_vtxIndx.end(), nSVVtx) != scoutingmuonVtx_vtxIndx.end()){
          //Muon matched
          nMuonMatch++;
          TLorentzVector scoutingmuonVtx_p4;
          scoutingmuonVtx_p4.SetPtEtaPhiM(ScoutingMuonVtx_pt[i], ScoutingMuonVtx_eta[i], ScoutingMuonVtx_phiCorr[i], 0.10566);
          svVtx_p4 += scoutingmuonVtx_p4;
        }
      }
      if(nMuonMatch>0) svVtx_mass = svVtx_p4.M();

      SVVtx_mass.push_back(svVtx_mass);
      SVVtx_nMuon.push_back(nMuonMatch);

      nSVVtx++;
    }
  }

  // Warning: Large gen multiplicity since no pruning done
  nGenPart = 0;
  if(genParticles.isValid()){
    for(auto &genParticle: *genParticles){
      GenPart_pt.push_back(genParticle.pt());
      GenPart_eta.push_back(genParticle.eta());
      GenPart_phi.push_back(genParticle.phi());
      TLorentzVector genParticle_p4;
      genParticle_p4.SetPtEtaPhiE(genParticle.pt(), genParticle.eta(), genParticle.phi(), genParticle.energy());
      GenPart_m.push_back(genParticle_p4.M());
      GenPart_pdgId.push_back(genParticle.pdgId());
      GenPart_status.push_back(genParticle.status());
      int GenPart_genPartIdxMother_temp = -1;
      if(genParticle.numberOfMothers() > 0) GenPart_genPartIdxMother_temp = static_cast<int>(genParticle.mother(0)->pdgId());
      GenPart_genPartIdxMother.push_back(GenPart_genPartIdxMother_temp);

      nGenPart++;
    }
  }

  Events->Fill();
  clearVars();
}

// ------------ Clear vectors  ------------
void Ntuplizer::clearVars(){
  ScoutingMuonNoVtx_pt.clear();
  ScoutingMuonNoVtx_eta.clear();
  ScoutingMuonNoVtx_phi.clear();
  ScoutingMuonNoVtx_phiCorr.clear();
  ScoutingMuonNoVtx_m.clear();
  ScoutingMuonNoVtx_charge.clear();
  ScoutingMuonNoVtx_normalizedChi2.clear();
  ScoutingMuonNoVtx_trkchi2.clear();
  ScoutingMuonNoVtx_trkndof.clear();
  ScoutingMuonNoVtx_trkdxy.clear();
  ScoutingMuonNoVtx_trkdz.clear();
  ScoutingMuonNoVtx_trkqoverp.clear();
  ScoutingMuonNoVtx_trklambda.clear();
  ScoutingMuonNoVtx_trkpt.clear();
  ScoutingMuonNoVtx_trkphi.clear();
  ScoutingMuonNoVtx_trketa.clear();
  ScoutingMuonNoVtx_trkqoverpError.clear();
  ScoutingMuonNoVtx_trklambdaError.clear();
  ScoutingMuonNoVtx_trkdxyError.clear();
  ScoutingMuonNoVtx_trkdzError.clear();
  ScoutingMuonNoVtx_trkphiError.clear();
  ScoutingMuonNoVtx_trkdsz.clear();
  ScoutingMuonNoVtx_trkdszError.clear();
  ScoutingMuonNoVtx_trkvx.clear();
  ScoutingMuonNoVtx_trkvy.clear();
  ScoutingMuonNoVtx_trkvz.clear();
  ScoutingMuonNoVtx_vtxIndx.clear();
  ScoutingMuonNoVtx_isGlobal.clear();
  ScoutingMuonNoVtx_isTracker.clear();
  ScoutingMuonNoVtx_isStandalone.clear();


  ScoutingMuonVtx_pt.clear();
  ScoutingMuonVtx_eta.clear();
  ScoutingMuonVtx_phi.clear();
  ScoutingMuonVtx_phiCorr.clear();
  ScoutingMuonVtx_m.clear();
  ScoutingMuonVtx_charge.clear();
  ScoutingMuonVtx_normalizedChi2.clear();
  ScoutingMuonVtx_trkchi2.clear();
  ScoutingMuonVtx_trkndof.clear();
  ScoutingMuonVtx_trkdxy.clear();
  ScoutingMuonVtx_trkdz.clear();
  ScoutingMuonVtx_trkqoverp.clear();
  ScoutingMuonVtx_trklambda.clear();
  ScoutingMuonVtx_trkpt.clear();
  ScoutingMuonVtx_trkphi.clear();
  ScoutingMuonVtx_trketa.clear();
  ScoutingMuonVtx_trkqoverpError.clear();
  ScoutingMuonVtx_trklambdaError.clear();
  ScoutingMuonVtx_trkdxyError.clear();
  ScoutingMuonVtx_trkdzError.clear();
  ScoutingMuonVtx_trkphiError.clear();
  ScoutingMuonVtx_trkdsz.clear();
  ScoutingMuonVtx_trkdszError.clear();  
  ScoutingMuonVtx_trkvx.clear();
  ScoutingMuonVtx_trkvy.clear();
  ScoutingMuonVtx_trkvz.clear();
  ScoutingMuonVtx_vtxIndx.clear();  
  ScoutingMuonVtx_isGlobal.clear();
  ScoutingMuonVtx_isTracker.clear();
  ScoutingMuonVtx_isStandalone.clear();

  SVNoVtx_x.clear();
  SVNoVtx_y.clear(); 
  SVNoVtx_z.clear();
  SVNoVtx_xError.clear();
  SVNoVtx_yError.clear();
  SVNoVtx_zError.clear();
  SVNoVtx_trksize.clear();
  SVNoVtx_chi2.clear();
  SVNoVtx_ndof.clear();
  SVNoVtx_isvalidvtx.clear();
  SVNoVtx_dxy.clear();
  SVNoVtx_dxySig.clear();
  SVNoVtx_dlen.clear();
  SVNoVtx_dlenSig.clear();
  SVNoVtx_mass.clear();
  SVNoVtx_nMuon.clear();

  SVVtx_x.clear();
  SVVtx_y.clear();
  SVVtx_z.clear();
  SVVtx_xError.clear();
  SVVtx_yError.clear(); 
  SVVtx_zError.clear();
  SVVtx_trksize.clear();
  SVVtx_chi2.clear();
  SVVtx_ndof.clear();
  SVVtx_isvalidvtx.clear();
  SVVtx_dxy.clear();
  SVVtx_dxySig.clear();
  SVVtx_dlen.clear();
  SVVtx_dlenSig.clear();
  SVVtx_mass.clear();
  SVVtx_nMuon.clear();

  GenPart_pt.clear();
  GenPart_eta.clear();
  GenPart_phi.clear();
  GenPart_m.clear();
  GenPart_pdgId.clear();
  GenPart_status.clear();
  GenPart_genPartIdxMother.clear();
}

// ------------ method called once each job just before starting event loop  ------------
void Ntuplizer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void Ntuplizer::endJob() {
  // please remove this method if not needed
}

void Ntuplizer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {

  HLTConfigProvider hltConfig;
  bool changedConfig = false;
  hltConfig.init(iRun, iSetup, "HLT", changedConfig);

  hltPathIndex_.clear();
  for (const auto &p : hltPaths_) hltPathIndex_[p] = -1;

  for (const auto &patternStr : hltPaths_) {
    TString pattern(patternStr);
    for (size_t j = 0; j < hltConfig.triggerNames().size(); ++j) {
      const std::string &pathName = hltConfig.triggerNames()[j];
      if (TString(pathName).Contains(pattern)) {
        edm::LogPrint("Ntuplizer") << "Found HLT path " << pathName
                                  << " matching " << patternStr
                                  << " with index " << j;
        hltPathIndex_[patternStr] = static_cast<int>(j);
        break; // stop at first match
      }
    }
  }
}

void Ntuplizer::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {
  // please remove this method if not needed
}

void Ntuplizer::beginLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup) {
  // please remove this method if not needed
}

void Ntuplizer::endLuminosityBlock(const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup) {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("ctfWithMaterialTracks"));
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntuplizer);