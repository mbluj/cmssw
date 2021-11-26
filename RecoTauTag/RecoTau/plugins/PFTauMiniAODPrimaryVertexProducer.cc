#include "RecoTauTag/RecoTau/interface/PFTauPrimaryVertexProducerBase.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

/// MiniAOD implementation of the PFTauPrimaryVertexProducer plugin
class PFTauMiniAODPrimaryVertexProducer final : public PFTauPrimaryVertexProducerBase {
public:
  explicit PFTauMiniAODPrimaryVertexProducer(const edm::ParameterSet &iConfig);
  ~PFTauMiniAODPrimaryVertexProducer() override;

  void beginEvent(const edm::Event &, const edm::EventSetup &) override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

protected:
  void nonTauTracksInPV(const reco::VertexRef &,
                        const std::vector<edm::Ptr<reco::TrackBase> > &,
                        std::vector<const reco::Track *> &) override;

private:
  void nonTauTracksInPVFromPackedCands(const size_t &,
                                       const pat::PackedCandidateCollection &,
                                       const std::vector<edm::Ptr<reco::TrackBase> > &,
                                       std::vector<const reco::Track *> &,
                                       bool checkPt = true);

  edm::EDGetTokenT<pat::PackedCandidateCollection> packedCandsToken_, lostCandsToken_, eleKFTracksToken_;
  edm::Handle<pat::PackedCandidateCollection> packedCands_, lostCands_, eleKFTracks_;
  bool useEleKFTracks_;
};

PFTauMiniAODPrimaryVertexProducer::PFTauMiniAODPrimaryVertexProducer(const edm::ParameterSet &iConfig)
    : PFTauPrimaryVertexProducerBase::PFTauPrimaryVertexProducerBase(iConfig),
      packedCandsToken_(
          consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("packedCandidatesTag"))),
      lostCandsToken_(
          consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("lostCandidatesTag"))),
      eleKFTracksToken_(
          consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("eleKFTracksTag"))),
      useEleKFTracks_(iConfig.getParameter<bool>("useEleKFTracks")) {}

PFTauMiniAODPrimaryVertexProducer::~PFTauMiniAODPrimaryVertexProducer() {}

void PFTauMiniAODPrimaryVertexProducer::beginEvent(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  //Get candidate collections
  iEvent.getByToken(packedCandsToken_, packedCands_);
  iEvent.getByToken(lostCandsToken_, lostCands_);
  if (useEleKFTracks_)
    iEvent.getByToken(eleKFTracksToken_, eleKFTracks_);
}

void PFTauMiniAODPrimaryVertexProducer::nonTauTracksInPV(const reco::VertexRef &thePVRef,
                                                         const std::vector<edm::Ptr<reco::TrackBase> > &tauTracks,
                                                         std::vector<const reco::Track *> &nonTauTracks) {
  //Find non-tau tracks associated to thePV
  //PackedCandidates first...
  if (packedCands_.isValid()) {
    nonTauTracksInPVFromPackedCands(thePVRef.key(), *packedCands_, tauTracks, nonTauTracks);
  }
  //then lostCandidates
  if (lostCands_.isValid()) {
    nonTauTracksInPVFromPackedCands(thePVRef.key(), *lostCands_, tauTracks, nonTauTracks);
  }
  //finally electron KF-track lostCandidates
  if (useEleKFTracks_ && eleKFTracks_.isValid()) {
    nonTauTracksInPVFromPackedCands(thePVRef.key(), *eleKFTracks_, tauTracks, nonTauTracks, false);
  }
}

void PFTauMiniAODPrimaryVertexProducer::nonTauTracksInPVFromPackedCands(
    const size_t &thePVkey,
    const pat::PackedCandidateCollection &cands,
    const std::vector<edm::Ptr<reco::TrackBase> > &tauTracks,
    std::vector<const reco::Track *> &nonTauTracks,
    bool checkPt) {
  //Find candidates/tracks associated to thePV
  for (const auto &cand : cands) {
    //MB: Skip candidates with ill-defined momentum as they return ill-defined tracks (why it happens?)
    if (!std::isfinite(cand.pt()))  //MB: it is enough to check just pt (?)
      continue;
    if (cand.vertexRef().isNull())
      continue;
    int quality = cand.pvAssociationQuality();
    if (cand.vertexRef().key() != thePVkey ||
        (quality != pat::PackedCandidate::UsedInFitTight && quality != pat::PackedCandidate::UsedInFitLoose))
      continue;
    const reco::Track *track = cand.bestTrack();
    if (track == nullptr)
      continue;
    //Remove signal (tau) tracks
    //MB: Only deltaR deltaPt overlap removal possible (?)
    //MB: It should be fine as pat objects stores same track info with same presision
    bool matched = false;
    for (const auto &tauTrack : tauTracks) {
      if (std::abs(tauTrack->eta() - track->eta()) < 0.005 &&
          std::abs(deltaPhi(tauTrack->phi(), track->phi())) < 0.005 &&
          (!checkPt || std::abs(tauTrack->pt() / track->pt() - 1.) < 0.005)) {
        matched = true;
        break;
      }
    }
    if (useEleKFTracks_ && std::abs(cand.pdgId()) == 11) {
      matched = true;
    }
    if (!matched)
      nonTauTracks.push_back(track);
  }
}

void PFTauMiniAODPrimaryVertexProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  auto desc = PFTauPrimaryVertexProducerBase::getDescriptionsBase();
  desc.add<edm::InputTag>("lostCandidatesTag", edm::InputTag("lostTracks"));
  desc.add<edm::InputTag>("packedCandidatesTag", edm::InputTag("packedPFCandidates"));
  desc.add<bool>("useEleKFTracks", false);
  desc.add<edm::InputTag>("eleKFTracksTag", edm::InputTag("lostTracks:eleTracks"));

  descriptions.add("pfTauMiniAODPrimaryVertexProducer", desc);
}

DEFINE_FWK_MODULE(PFTauMiniAODPrimaryVertexProducer);
