// -*- C++ -*-
//
// Package:    DiHadronCorrelator/DiHadronCorrAnalyzer
// Class:      DiHadronCorrAnalyzer
//
/**\class DiHadronCorrAnalyzer DiHadronCorrAnalyzer.cc DiHadronCorrelator/DiHadronCorrAnalyzer/plugins/DiHadronCorrAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yousen Zhang
//         Created:  Mon, 27 Jan 2020 20:48:24 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <utility>      // std::move

// user include files

#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "FWCore/Common/interface/TriggerNames.h"

#include "TH1F.h"
#include "TH2F.h"

//
// class declaration
//
#include "myVtxObj.h"

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;
using reco::VertexCollection;
using std::vector;

using myVtx::vtxAtt;
using myVtx::hadronAtt;

class DiHadronCorrAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit DiHadronCorrAnalyzer(const edm::ParameterSet&);
      ~DiHadronCorrAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<reco::VertexCollection> verticesToken_; // used to select what primary vertices to read from config

      edm::Service<TFileService> fs;

      vector<double> ptBinEdges_;
      vector<double> zVtxBinEdges_;
      vector<unsigned> nTrkOfflineBinEdges_;

      unsigned int nMixEvts_;

      vector<TH1*> vec_hPt_;
      vector<TH1*> vec_hEta_;
      vector<TH1*> vec_hPhi_;
      vector<TH1*> vec_nTrkOffline_;
      vector<TH2*> vec_nTrkOfflineVsVtxZ_;
      vector<vector<TH2*>> vec_vec_h2DSignal;
      vector<vector<TH2*>> vec_vec_h2DBackground;

      vector<vector<vtxAtt>> vec_vec_vtxAtt; // event<vtx<hadron>>
};

//
// constants, enums and typedefs
//
const double PI = TMath::Pi();
const unsigned int nEtaBins_ = 33;
const unsigned int nPhiBins_ = 32;

const double etaBegin_ = -4.95;
const double etaEnd_ = 4.95;

const double phiBegin_ = -(0.5+1.0/32)*PI;
const double phiEnd_ = (1.5-1.0/32)*PI;

//
// static data member definitions
//

//
// constructors and destructor
//
DiHadronCorrAnalyzer::DiHadronCorrAnalyzer(const edm::ParameterSet& iConfig)
 :
  tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
  //tracksToken_(consumes<TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  verticesToken_(consumes<VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("primaryVertices"))),
  //verticesToken_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertices"))),
  ptBinEdges_(iConfig.getParameter< vector<double> >("ptBinEdges")),
  zVtxBinEdges_(iConfig.getParameter< vector<double> >("zVtxBinEdges")),
  nTrkOfflineBinEdges_(iConfig.getParameter< vector<unsigned> >("nTrkOfflineBinEdges")),
  nMixEvts_(iConfig.getParameter<unsigned int>("nMixEvts"))
{
   //now do what ever initialization is needed
  const unsigned int nPt = ptBinEdges_.size() - 1;
  const unsigned int nzVtx = zVtxBinEdges_.size() - 1;
  const unsigned int nu_nTrkOffline = nTrkOfflineBinEdges_.size() - 1;

  vec_vec_h2DSignal.resize(nu_nTrkOffline);
  vec_vec_h2DBackground.resize(nu_nTrkOffline);

  vec_hPt_.resize(nu_nTrkOffline);
  vec_hEta_.resize(nu_nTrkOffline);
  vec_hPhi_.resize(nu_nTrkOffline);

  if(nPt<1) throw cms::Exception("DiHadronCorrAnalyzer") << "length of ptBinEdges is less than 2 !!!" << std::endl;
  if(nzVtx<1) throw cms::Exception("DiHadronCorrAnalyzer") << "length of zVtxBinEdges is less than 2 !!!" << std::endl;
  if(nu_nTrkOffline<1) throw cms::Exception("DiHadronCorrAnalyzer") << "length of nTrkOfflineBinEdges is less than 2 !!!" << std::endl;

  for(unsigned int itrk=0; itrk<nu_nTrkOffline; ++itrk){

    vec_vec_h2DSignal.at(itrk).resize(nPt);
    vec_vec_h2DBackground.at(itrk).resize(nPt);

    vec_hPt_.at(itrk) = fs->make<TH1F>(Form("hPt_trk%d", itrk), Form("Ntrk_%u_%u;p_{T} (GeV);Entries", 
          nTrkOfflineBinEdges_[itrk], nTrkOfflineBinEdges_[itrk+1]), 100, 0, 10);
    vec_hEta_.at(itrk) = fs->make<TH1F>(Form("hEta_trk%d", itrk), Form("Ntrk_%u_%u;#eta (GeV);Entries",
          nTrkOfflineBinEdges_[itrk], nTrkOfflineBinEdges_[itrk+1]), 100, -5, 5);
    vec_hPhi_.at(itrk) = fs->make<TH1F>(Form("hPhi_trk%d", itrk), Form("Ntrk_%u_%u;#phi (GeV);Entries",
          nTrkOfflineBinEdges_[itrk], nTrkOfflineBinEdges_[itrk+1]), 100, -PI, PI);

    for(unsigned int ipt=0; ipt<nPt; ++ipt){

      vec_vec_h2DSignal.at(itrk).at(ipt) =
          fs->make<TH2F>(Form("hSignalNtrk%dPt%d", itrk, ipt),
            Form("Ntrk_%u_%u_pT_%.1f_%.1fGeV;#Delta#eta;#Delta#phi", nTrkOfflineBinEdges_[itrk], 
              nTrkOfflineBinEdges_[itrk+1], ptBinEdges_[ipt], ptBinEdges_[ipt+1]),
              nEtaBins_, etaBegin_, etaEnd_, nPhiBins_, phiBegin_, phiEnd_);

      vec_vec_h2DBackground.at(itrk).at(ipt) =
          fs->make<TH2F>(Form("hBackgroundNtrk%dPt%d", itrk, ipt),
            Form("Ntrk_%u_%u_pT_%.1f_%.1fGeV;#Delta#eta;#Delta#phi", nTrkOfflineBinEdges_[itrk], 
              nTrkOfflineBinEdges_[itrk+1], ptBinEdges_[ipt], ptBinEdges_[ipt+1]),
              nEtaBins_, etaBegin_, etaEnd_, nPhiBins_, phiBegin_, phiEnd_);

    }
  }
}


DiHadronCorrAnalyzer::~DiHadronCorrAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DiHadronCorrAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  edm::Handle< reco::VertexCollection > vertices;
  iEvent.getByToken(verticesToken_, vertices);
  if(!vertices->size()) { 
    //std::cout << "Invalid or empty vertex collection!" << std::endl;
    return;
  }

  edm::Handle< reco::TrackCollection > tracks;
  iEvent.getByToken(tracksToken_, tracks);
  if(!tracks->size()) {
  //  std::cout << "Invalid or empty track collection!" << std::endl;
    return;
  }

  unsigned int nVtx = 0;
  for(unsigned int iv=0; iv<vertices->size(); iv++) {
    const reco::Vertex & vtx = (*vertices)[iv];
    if(!vtx.isFake() && vtx.tracksSize()>=2 && fabs(vtx.z())<15) 
      nVtx++;
  }

  if(!nVtx) return;

  vector<vtxAtt> vec_vtx;
  vec_vtx.reserve(nVtx);

  vector<hadronAtt> vec_trig_tmp;
  vector<hadronAtt> vec_assc_tmp;
  vec_trig_tmp.reserve(100);
  vec_assc_tmp.reserve(100);
  for(unsigned int iv=0; iv<vertices->size(); ++iv){
    const reco::Vertex & vtx = (*vertices)[iv];
    // must contain more than 2 tracks and z within |15cm|
    bool passVtx = !vtx.isFake() && vtx.tracksSize()>=2 && fabs(vtx.z())<15;
    if(!passVtx) continue;
    int iz = 0;
    int itrk = 0;

    iz = myVtx::getIdx(zVtxBinEdges_, vtx.z());

    // loop tracks, to be done
    //
    //
    unsigned int Ntrkoffline = 0;

    unsigned int trk_id = 0;
    for (reco::Vertex::trackRef_iterator iTrack = vtx.tracks_begin(); 
         iTrack != vtx.tracks_end(); iTrack++) {
     
       reco::TrackRef track = iTrack->castTo<reco::TrackRef>();

       trk_id++;

       const float dxy = track->dxy(vtx.position());
       const float dxyoerr = dxy/sqrt(track->d0Error()*track->d0Error()+ vtx.xError()*vtx.yError());
       const float dz = track->dz(vtx.position());
       const float dzoerr = dz/sqrt(track->dzError()*track->dzError()+ vtx.zError()*vtx.zError());

       // require high purity tracks
       // require pTerr/pT<0.1
       // require zDCA significance w.r.t. the vertex < 3
       // require xyDCA significance w.r.t. the vertex < 3
       // pt >= 0.3, |eta| < 2.4
       bool isPrimaryTrack = track->quality( reco::TrackBase::highPurity ) &&
           track->ptError()/track->pt() < 0.1 && dzoerr < 3 && dxyoerr < 3 &&
           track->pt() >= 0.3 && std::fabs(track->eta()) < 2.4;

       if(!isPrimaryTrack) continue;

       // primary track and pt > 0.4
       if(isPrimaryTrack && track->pt() >= 0.4) 
         Ntrkoffline++;

       if(isPrimaryTrack && track->pt() < 3.) 
         vec_assc_tmp.push_back(hadronAtt(trk_id, -99,
               track->eta(), track->phi()) );

       auto ipt = myVtx::getIdx(ptBinEdges_, track->pt());
       if(ipt < 0) continue;
       if(isPrimaryTrack)
         vec_trig_tmp.push_back(hadronAtt(trk_id, ipt,
               track->eta(), track->phi()) );

    }

    itrk = myVtx::getIdx(nTrkOfflineBinEdges_, Ntrkoffline);


    if(iz<0) continue;
    if(itrk<0) continue;

    for(auto& itrack : vec_trig_tmp){
      vec_hEta_.at(itrk)->Fill(itrack.eta());
      vec_hPhi_.at(itrk)->Fill(itrack.phi());
    }

    vtxAtt tmp(iv, iz, itrk);
    tmp.addVecHadron( vec_assc_tmp, "assc");
    tmp.addVecHadron( vec_trig_tmp, "trig");

    vec_trig_tmp.clear();
    vec_assc_tmp.clear();

    vec_vtx.push_back( std::move(tmp) );
  }

  // start correlating trig and assc
  for(auto& ivtx : vec_vtx){
    for(auto& itrig : ivtx.vec_trig()){
      for(auto& iassc : ivtx.vec_assc()){
        // do not use the same track
        if(itrig.id() == iassc.id()) continue;

        float deltaPhi = reco::deltaPhi(iassc.phi(), itrig.phi());
        float deltaEta = iassc.eta() - itrig.eta();

        // reject close tracks
        //if(std::fabs(deltaEta) < 0.03 && std::fabs(deltaPhi) < 0.03) continue;

        // rerange deltaPhi to -0.5 PI to 1.5 PI
        if(deltaPhi < phiBegin_) deltaPhi += 2*PI;

        vec_vec_h2DSignal.at( ivtx.itrkBin() ).at( itrig.ipt() )->Fill(
            deltaEta, deltaPhi);
      }
    }
  }

  // push to the event vector for background study
  vec_vec_vtxAtt.push_back( std::move(vec_vtx) );
}


// ------------ method called once each job just before starting event loop  ------------
void
DiHadronCorrAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
DiHadronCorrAnalyzer::endJob()
{
  // need more time to do event mixing
  // for each event, let us find 10 events (vertices) to evalute mixed background
  for(unsigned long int ievt=0; ievt< vec_vec_vtxAtt.size(); ievt++){
    // read vertex collection;
    auto& ivec_vtx = vec_vec_vtxAtt[ievt];
    // read each vertex
    for(auto& ivtx : ivec_vtx){
      // read the collection of trigger particles in the current event
      for(auto& itrig : ivtx.vec_trig()){
        // find the associated particles in mixed events
        unsigned int imix = 0;
        for(unsigned int ievt_mix = ievt+1; 
            ievt_mix < vec_vec_vtxAtt.size(); ievt_mix++){
          for(auto& ivtx_mix : vec_vec_vtxAtt[ievt_mix]){

            // match the same z range and ntrk range
            bool isFound= ivtx.iz() == ivtx_mix.iz() && ivtx.itrkBin() == ivtx_mix.itrkBin();
            if(!isFound) continue;

            // loop asscoiated particles in mixed events
            for(auto& iassc : ivtx_mix.vec_assc()){

              float deltaPhi = reco::deltaPhi(iassc.phi(), itrig.phi());
              float deltaEta = iassc.eta() - itrig.eta();

              // reject close tracks
              //if(std::fabs(deltaEta) < 0.03 && std::fabs(deltaPhi) < 0.03) continue;

              // rerange deltaPhi to -0.5 PI to 1.5 PI
              if(deltaPhi < phiBegin_) deltaPhi += 2*PI;

              vec_vec_h2DBackground.at( ivtx.itrkBin() ).at( itrig.ipt() )->Fill(
                  deltaEta, deltaPhi);
            }

            // mix more than ten times, break
            if(++imix == nMixEvts_) break;
          }
          if(imix == nMixEvts_) break;
        }
      }
    }
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DiHadronCorrAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  //edm::ParameterSetDescription desc;
  //desc.setUnknown();
  //descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  edm::ParameterSetDescription desc;
  desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.addUntracked<edm::InputTag>("primaryVertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<vector<unsigned>>("nTrkOfflineBinEdges", {0, 20, 40, 80, 100, 10000000});
  desc.add<vector<double>>("ptBinEdges", {0.3, 0.6, 0.9, 1.2,
      1.5, 1.8, 2.1, 2.4, 2.7, 3.0});
  desc.add<vector<double>>("zVtxBinEdges", { -15, -13, -11, -9, -7, -5, -3, -1,
                   1, 3, 5, 7, 9, 11, 13, 15 });
  desc.add<unsigned int>("nMixEvts", 10);
  //descriptions.addDefault(desc);
  descriptions.add("dihadroncorrelator", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DiHadronCorrAnalyzer);
