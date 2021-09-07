// -*- C++ -*-
//
// Package:    PhysicsTools/NanoAOD
// Class:      TauTrackTableProducer
//
/**\class TauTrackTableProducer TauTrackTableProducer.cc PhysicsTools/TauTrackTableProducer/plugins/TauTrackTableProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michal Bluj
//         Created:  Mon, 28 Aug 2017 09:26:39 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

//
// helper data structure
//

namespace {
  class TrackData {
  public:
    std::vector<float> pt_, eta_, phi_, dxy_, dz_, normChi2_;
    std::vector<unsigned int> nPxlHits_, nHits_, pvAssocQ_;
    std::vector<int> genTauIdx_, tauIdx_, vtxIdx_, charge_;
    std::vector<bool> isLost_;

  public:
    size_t size() { return pt_.size(); };
    void fill(const reco::Track& t,
              int genTauIdx,
              int tauIdx = -1,
              math::XYZPoint pvPos = math::XYZPoint(),
              unsigned int pvAssocQ = 0,
              int vtxIdx = -1,
              bool isLost = false) {
      pt_.push_back(t.pt());
      eta_.push_back(t.eta());
      phi_.push_back(t.phi());
      dxy_.push_back(t.dxy(pvPos));
      dz_.push_back(t.dz(pvPos));
      normChi2_.push_back(t.normalizedChi2());
      charge_.push_back(t.charge());
      nPxlHits_.push_back(t.hitPattern().numberOfValidPixelHits());
      nHits_.push_back(t.hitPattern().numberOfValidHits());
      genTauIdx_.push_back(genTauIdx);
      tauIdx_.push_back(tauIdx);
      vtxIdx_.push_back(vtxIdx);
      pvAssocQ_.push_back(pvAssocQ);
      isLost_.push_back(isLost);
    };
    void fill(const pat::PackedCandidate& p,
              int genTauIdx,
              int tauIdx = -1,
              math::XYZPoint pvPos = math::XYZPoint(),
              bool isLost = false) {
      if (!p.hasTrackDetails())
        return;
      fill(*p.bestTrack(), genTauIdx, tauIdx, pvPos, p.pvAssociationQuality(), p.vertexRef().key(), isLost);
    };
  };
}  // namespace

//
// class declaration
//

class TauTrackTableProducer : public edm::stream::EDProducer<> {
public:
  explicit TauTrackTableProducer(const edm::ParameterSet&);
  ~TauTrackTableProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------

  const edm::EDGetTokenT<std::vector<reco::Vertex>> pvs_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> pfCands_, lostTracks_;
  const edm::EDGetTokenT<edm::RefVector<std::vector<pat::Tau>>> taus_;
  const edm::EDGetTokenT<reco::GenParticleCollection> genTaus_;
  const std::string trkName_;
  const double dR2max_;
};

//
// constructors and destructor
//
TauTrackTableProducer::TauTrackTableProducer(const edm::ParameterSet& params)
    : pvs_(consumes<std::vector<reco::Vertex>>(params.getParameter<edm::InputTag>("pvSrc"))),
      pfCands_(consumes<pat::PackedCandidateCollection>(params.getParameter<edm::InputTag>("pfCandidatesSrc"))),
      lostTracks_(consumes<pat::PackedCandidateCollection>(params.getParameter<edm::InputTag>("lostTracksSrc"))),
      taus_(consumes<edm::RefVector<std::vector<pat::Tau>>>(params.getParameter<edm::InputTag>("tauSrc"))),
      genTaus_(consumes<reco::GenParticleCollection>(params.getParameter<edm::InputTag>("genTauSrc"))),
      trkName_(params.getParameter<std::string>("trkName")),
      dR2max_(std::pow(params.getParameter<double>("dRMatch"), 2))

{
  produces<nanoaod::FlatTable>("");
}

TauTrackTableProducer::~TauTrackTableProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to produce the data  ------------

void TauTrackTableProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //using namespace edm;
  edm::Handle<std::vector<reco::Vertex>> pvsIn;
  iEvent.getByToken(pvs_, pvsIn);
  edm::Handle<pat::PackedCandidateCollection> pfCandsIn, lostTracksIn;
  iEvent.getByToken(pfCands_, pfCandsIn);
  iEvent.getByToken(lostTracks_, lostTracksIn);
  edm::Handle<edm::RefVector<std::vector<pat::Tau>>> tausIn;
  iEvent.getByToken(taus_, tausIn);
  edm::Handle<reco::GenParticleCollection> genTausIn;
  iEvent.getByToken(genTaus_, genTausIn);

  TrackData trackData;
  //pfCands' tracks
  for (size_t iCand = 0; iCand < pfCandsIn->size(); ++iCand) {
    const pat::PackedCandidate& cand = (*pfCandsIn)[iCand];
    if (!(cand.pt() > 0.5 && cand.hasTrackDetails()))
      continue;
    const reco::Track& track = *cand.bestTrack();
    if (!(track.pt() > 0.5))
      continue;
    float dR2match = dR2max_;
    int genTauIdx = -1;
    for (size_t iGenTau = 0; iGenTau < genTausIn->size(); ++iGenTau) {
      const reco::GenParticle& genTau = (*genTausIn)[iGenTau];
      float dR2 = deltaR2(track, genTau);
      if (dR2 > dR2match)
        continue;
      dR2match = dR2;
      genTauIdx = iGenTau;
    }
    if (genTauIdx >= 0) {
      int tauIdx = -1;
      //matching to recoTaus ref
      for (size_t iTau = 0; iTau < tausIn->size(); ++iTau) {
        const pat::Tau& tau = *(*tausIn)[iTau];
        for (const reco::CandidatePtr& tauSignalCand : tau.signalChargedHadrCands()) {
          if (tauSignalCand.key() == iCand) {  //same collection assumed
            tauIdx = iTau;
            break;
          }
        }
        if (tauIdx >= 0)
          break;
      }

      trackData.fill(cand, genTauIdx, tauIdx, (*pvsIn)[0].position(), false);
    }
  }
  //lost-tracks
  for (size_t iCand = 0; iCand < lostTracksIn->size(); ++iCand) {
    const pat::PackedCandidate& cand = (*lostTracksIn)[iCand];
    if (!(cand.pt() > 0.5 && cand.hasTrackDetails()))
      continue;
    const reco::Track& track = *cand.bestTrack();
    if (!(track.pt() > 0.5))
      continue;
    float dR2match = dR2max_;
    int genTauIdx = -1;
    for (size_t iGenTau = 0; iGenTau < genTausIn->size(); ++iGenTau) {
      const reco::GenParticle& genTau = (*genTausIn)[iGenTau];
      float dR2 = deltaR2(track, genTau);
      if (dR2 > dR2match)
        continue;
      dR2match = dR2;
      genTauIdx = iGenTau;
    }
    if (genTauIdx >= 0) {
      int tauIdx = -1;
      //matching to recoTaus by dR
      for (size_t iTau = 0; iTau < tausIn->size(); ++iTau) {
        const pat::Tau& tau = *(*tausIn)[iTau];
        float signalConeR2 = std::pow(std::clamp(3. / std::max(1., tau.pt()), 0.05, 0.1), 2);
        if (deltaR2(cand, tau) < signalConeR2) {
          tauIdx = iTau;
          break;
        }
      }

      trackData.fill(cand, genTauIdx, tauIdx, (*pvsIn)[0].position(), true);
    }
  }

  auto trkTable = std::make_unique<nanoaod::FlatTable>(trackData.size(), trkName_, false);
  trkTable->addColumn<float>("pt", trackData.pt_, "track pt", -1);  //12
  trkTable->addColumn<float>("eta", trackData.eta_, "track eta", 12);
  trkTable->addColumn<float>("phi", trackData.phi_, "track phi", 12);
  trkTable->addColumn<int>("charge", trackData.charge_, "track charge", 10);
  trkTable->addColumn<float>("dz", trackData.dz_, "track dz wrt PV", 14);
  trkTable->addColumn<float>("dxy", trackData.dxy_, "track dxy wrt PV", 10);
  trkTable->addColumn<float>("normChi2", trackData.normChi2_, "track chi2/ndof", 8);
  trkTable->addColumn<int>("nPxlHits", trackData.nPxlHits_, "no. of pixel hits", -1);
  trkTable->addColumn<int>("nHits", trackData.nHits_, "no. of hits", -1);
  trkTable->addColumn<int>("pvAssocQ", trackData.pvAssocQ_, "association quality to vertex", -1);
  trkTable->addColumn<int>("genTauIdx", trackData.genTauIdx_, "index of matched genTau", -1);
  trkTable->addColumn<int>("tauIdx", trackData.tauIdx_, "index of matched tau", -1);
  trkTable->addColumn<int>("vtxIdx", trackData.vtxIdx_, "index of associated vertex", -1);
  trkTable->addColumn<bool>("isLost", trackData.isLost_, "is lost-track", -1);

  iEvent.put(std::move(trkTable));
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void TauTrackTableProducer::beginStream(edm::StreamID) {}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void TauTrackTableProducer::endStream() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TauTrackTableProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TauTrackTableProducer);
