/**
  \class    PATLeptonTimeLifeInfoUpdater
  \brief    Produces updated collection of pat::Electron's or pat::Muon's 

   The PATLeptonInfoUpdater produces analysis-level pat::Electron's or 
   pat::Muon's with updated time-life information.

  \author   Michal Bluj, NCBJ, Warsaw
*/

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalTrajectoryExtrapolatorToLine.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalImpactPointExtrapolator.h"

#include <cstring>

template <typename T>
class PATLeptonTimeLifeInfoUpdater : public edm::stream::EDProducer<> {
public:
  explicit PATLeptonTimeLifeInfoUpdater(const edm::ParameterSet&);
  ~PATLeptonTimeLifeInfoUpdater() override{};

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  //--- private utility methods
  const reco::Track* getTrack(const T&);
  void produceAndFillIPInfo(T&, const TransientTrackBuilder&, const reco::Vertex&);
  void produceAndFillSVInfo(T&, const TransientTrackBuilder&, const reco::Vertex&);

  //--- configuration parameters
  edm::EDGetTokenT<std::vector<T>> leptonsToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transTrackBuilderToken_;
  int pvChoice_;

  enum PVChoice { useFront = 0, useClosestInDz };
};

template <typename T>
PATLeptonTimeLifeInfoUpdater<T>::PATLeptonTimeLifeInfoUpdater(const edm::ParameterSet& cfg)
    : leptonsToken_(consumes<std::vector<T>>(cfg.getParameter<edm::InputTag>("src"))),
      pvToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("pvSource"))),
      transTrackBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      pvChoice_(cfg.getParameter<int>("pvChoice")) {
  produces<std::vector<T>>();
}

template <typename T>
void PATLeptonTimeLifeInfoUpdater<T>::produce(edm::Event& evt, const edm::EventSetup& es) {
  // Get leptons
  edm::Handle<std::vector<T>> inputLeptons;
  evt.getByToken(leptonsToken_, inputLeptons);

  auto outputLeptons = std::make_unique<std::vector<T>>();
  outputLeptons->reserve(inputLeptons->size());

  // Get the vertices
  edm::Handle<reco::VertexCollection> vertices;
  evt.getByToken(pvToken_, vertices);

  // Get transient track builder
  const TransientTrackBuilder& transTrackBuilder = es.getData(transTrackBuilderToken_);

  for (const auto& inputLepton : *inputLeptons) {
    T outputLepton(inputLepton);
    size_t pv_idx = 0;
    if (pvChoice_ == useClosestInDz && getTrack(inputLepton) != nullptr) {
      float dz_min = 999;
      size_t vtx_idx = 0;
      for (const auto& vtx : *vertices) {
        float dz_tmp = std::abs(getTrack(inputLepton)->dz(vtx.position()));
        if (dz_tmp < dz_min) {
          dz_min = dz_tmp;
          pv_idx = vtx_idx;
        }
        vtx_idx++;
      }
    }
    const reco::Vertex& pv = !vertices->empty() ? (*vertices)[pv_idx] : reco::Vertex();
    // Obtain IP vector and set related info into lepton
    produceAndFillIPInfo(outputLepton, transTrackBuilder, pv);

    // Fit SV and set related info in taus or do nothing for other lepton types
    produceAndFillSVInfo(outputLepton, transTrackBuilder, pv);

    // Set updated lepton to output collection
    outputLeptons->push_back(outputLepton);
  }  // end of lepton loop

  evt.put(std::move(outputLeptons));
}

template <>
const reco::Track* PATLeptonTimeLifeInfoUpdater<pat::Electron>::getTrack(const pat::Electron& electron) {
  return electron.gsfTrack().isNonnull() ? electron.gsfTrack().get() : nullptr;
}

template <>
const reco::Track* PATLeptonTimeLifeInfoUpdater<pat::Muon>::getTrack(const pat::Muon& muon) {
  return muon.innerTrack().isNonnull() ? muon.innerTrack().get() : nullptr;
}

template <>
const reco::Track* PATLeptonTimeLifeInfoUpdater<pat::Tau>::getTrack(const pat::Tau& tau) {
  const reco::Track* track = nullptr;
  if (tau.leadChargedHadrCand().isNonnull())
    track = tau.leadChargedHadrCand()->bestTrack();
  return track;
}

template <typename T>
void PATLeptonTimeLifeInfoUpdater<T>::produceAndFillIPInfo(T& lepton,
                                                           const TransientTrackBuilder& transTrackBuilder,
                                                           const reco::Vertex& pv) {
  GlobalPoint pca;
  GlobalError pca_cov;
  GlobalVector ip_vec;
  GlobalError ip_cov;
  float ipLength = -1;
  float ipLength_sig = -1;
  float bField_z = 0;
  const reco::Track* track = getTrack(lepton);
  if (track != nullptr) {
    bField_z = transTrackBuilder.field()->inInverseGeV(GlobalPoint(track->vx(), track->vy(), track->vz())).z();
    // Extrapolate track to the point closest to PV
    reco::TransientTrack transTrack = transTrackBuilder.build(track);
    AnalyticalImpactPointExtrapolator extrapolator(transTrack.field());
    TrajectoryStateOnSurface closestState =
        extrapolator.extrapolate(transTrack.impactPointState(), RecoVertex::convertPos(pv.position()));
    pca = closestState.globalPosition();
    pca_cov = closestState.cartesianError().position();
    ip_vec = GlobalVector(pca.x() - pv.x(), pca.y() - pv.y(), pca.z() - pv.z());
    ip_cov = pca_cov + GlobalError(pv.covariance());
    VertexDistance3D pca_dist;
    Measurement1D ip_mes = pca_dist.distance(pv, VertexState(pca, pca_cov));
    if (ip_vec.dot(GlobalVector(lepton.px(), lepton.py(), lepton.pz())) < 0)
      ip_mes = Measurement1D(-1. * ip_mes.value(), ip_mes.error());

    ipLength = ip_mes.value();
    ipLength_sig = ip_mes.significance();
  }
  // Store PCA info
  lepton.addUserFloat("refitPCA_x", pca.x());
  lepton.addUserFloat("refitPCA_y", pca.y());
  lepton.addUserFloat("refitPCA_z", pca.z());
  lepton.addUserFloat("refitPCACov_xx", pca_cov.cxx());
  lepton.addUserFloat("refitPCACov_yx", pca_cov.cyx());
  lepton.addUserFloat("refitPCACov_zx", pca_cov.czx());
  lepton.addUserFloat("refitPCACov_yy", pca_cov.cyy());
  lepton.addUserFloat("refitPCACov_zy", pca_cov.czy());
  lepton.addUserFloat("refitPCACov_zz", pca_cov.czz());
  // Store IP info
  lepton.addUserFloat("IP", ipLength);
  lepton.addUserFloat("IP_sig", ipLength_sig);
  lepton.addUserFloat("IP_x", ip_vec.x());
  lepton.addUserFloat("IP_y", ip_vec.y());
  lepton.addUserFloat("IP_z", ip_vec.z());
  lepton.addUserFloat("IPCov_xx", ip_cov.cxx());
  lepton.addUserFloat("IPCov_yx", ip_cov.cyx());
  lepton.addUserFloat("IPCov_zx", ip_cov.czx());
  lepton.addUserFloat("IPCov_yy", ip_cov.cyy());
  lepton.addUserFloat("IPCov_zy", ip_cov.czy());
  lepton.addUserFloat("IPCov_zz", ip_cov.czz());
  // Store leading track parameters and covariance at reference point
  if (track != nullptr) {
    lepton.addUserFloat("track_qoverp", track->parameter(reco::TrackBase::i_qoverp));
    lepton.addUserFloat("track_lambda", track->parameter(reco::TrackBase::i_lambda));
    lepton.addUserFloat("track_phi", track->parameter(reco::TrackBase::i_phi));
    lepton.addUserFloat("track_dxy", track->parameter(reco::TrackBase::i_dxy));
    lepton.addUserFloat("track_dsz", track->parameter(reco::TrackBase::i_dsz));
    for (size_t i = 0; i < reco::TrackBase::dimension; ++i)
      for (size_t j = i; j < reco::TrackBase::dimension; ++j)
        lepton.addUserFloat("trackCov_" + std::to_string(j + i), track->covariance(j, i));
  } else {
    lepton.addUserFloat("track_qoverp", 0);
    lepton.addUserFloat("track_lambda", 0);
    lepton.addUserFloat("track_phi", 0);
    lepton.addUserFloat("track_dxy", 0);
    lepton.addUserFloat("track_dsz", 0);
    for (size_t i = 0; i < reco::TrackBase::dimension; ++i)
      for (size_t j = i; j < reco::TrackBase::dimension; ++j)
        lepton.addUserFloat("trackCov_" + std::to_string(j + i), 0);
  }
  // ... and z coordinate of magnetic field at track reference point
  lepton.addUserFloat("bField_z", bField_z);
}

template <typename T>
void PATLeptonTimeLifeInfoUpdater<T>::produceAndFillSVInfo(T& lepton,
                                                           const TransientTrackBuilder& transTrackBuilder,
                                                           const reco::Vertex& pv) {}

template <>
void PATLeptonTimeLifeInfoUpdater<pat::Tau>::produceAndFillSVInfo(pat::Tau& tau,
                                                                  const TransientTrackBuilder& transTrackBuilder,
                                                                  const reco::Vertex& pv) {
  //MB: Set initially SV to the tau vertex (tau.vertex()) or the PV choosen?
  reco::Vertex sv;
  GlobalVector flight_vec;
  GlobalError flight_cov;
  float flightLength = -1;
  float flightLength_sig = -1;

  // Fit SV with tracks of charged tau decay products
  int fitOK = 0;
  if (tau.signalChargedHadrCands().size() + tau.signalLostTracks().size() > 1) {
    // Get tracks form tau signal charged candidates
    std::vector<reco::TransientTrack> transTrks;
    TransientVertex transVtx;
    for (const auto& cand : tau.signalChargedHadrCands()) {
      if (cand.isNull())
        continue;
      const reco::Track* track = cand->bestTrack();
      if (track != nullptr)
        transTrks.push_back(transTrackBuilder.build(track));
    }
    for (const auto& cand : tau.signalLostTracks()) {
      if (cand.isNull())
        continue;
      const reco::Track* track = cand->bestTrack();
      if (track != nullptr)
        transTrks.push_back(transTrackBuilder.build(track));
    }
    // Fit SV with KalmanVertexFitter
    fitOK = 1;
    KalmanVertexFitter kvf(true);
    if (transTrks.size() > 1) {
      try {
        transVtx = kvf.vertex(transTrks);
      } catch (...) {
        fitOK = -1;
      }
    } else {
      fitOK = -1;
    }
    if (!transVtx.hasRefittedTracks() || transVtx.refittedTracks().size() != transTrks.size())
      fitOK = -1;
    if (fitOK > 0) {
      sv = transVtx;
      // Get flight-length
      // Full PV->SV flight vector with its covariance
      flight_vec = GlobalVector(sv.x() - pv.x(), sv.y() - pv.y(), sv.z() - pv.z());
      flight_cov = transVtx.positionError() + GlobalError(pv.covariance());
      //MB: can be taken from tau itself (but with different fit of PV and SV) as follows:
      //flightLength = sqrt(tau.flightLength().mag2());
      //flightLength_sig = tau.flightLengthSig();
      VertexDistance3D sv_dist;
      Measurement1D flightLength_mes = sv_dist.signedDistance(pv, sv, GlobalVector(tau.px(), tau.py(), tau.pz()));
      flightLength = flightLength_mes.value();
      flightLength_sig = flightLength_mes.significance();
    }
  }
  // Embed SV info in tau
  // Store SV info
  tau.addUserInt("refitSV_status", fitOK);
  tau.addUserFloat("refitSV_x", sv.x());
  tau.addUserFloat("refitSV_y", sv.y());
  tau.addUserFloat("refitSV_z", sv.z());
  tau.addUserFloat("refitSV_chi2", sv.chi2());
  tau.addUserFloat("refitSV_ndof", sv.ndof());
  tau.addUserFloat("refitSVCov_xx", sv.covariance(0, 0));
  tau.addUserFloat("refitSVCov_yx", sv.covariance(1, 0));
  tau.addUserFloat("refitSVCov_zx", sv.covariance(2, 0));
  tau.addUserFloat("refitSVCov_yy", sv.covariance(1, 1));
  tau.addUserFloat("refitSVCov_zy", sv.covariance(2, 1));
  tau.addUserFloat("refitSVCov_zz", sv.covariance(2, 2));
  // Store flight-length info
  tau.addUserFloat("flightLength", flightLength);
  tau.addUserFloat("flightLength_sig", flightLength_sig);
  tau.addUserFloat("flight_x", flight_vec.x());
  tau.addUserFloat("flight_y", flight_vec.y());
  tau.addUserFloat("flight_z", flight_vec.z());
  tau.addUserFloat("flightCov_xx", flight_cov.cxx());
  tau.addUserFloat("flightCov_yx", flight_cov.cyx());
  tau.addUserFloat("flightCov_zx", flight_cov.czx());
  tau.addUserFloat("flightCov_yy", flight_cov.cyy());
  tau.addUserFloat("flightCov_zy", flight_cov.czy());
  tau.addUserFloat("flightCov_zz", flight_cov.czz());
}

template <typename T>
void PATLeptonTimeLifeInfoUpdater<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // pat{Electron,Muon,Tau}TimeLifeInfoUpdater
  edm::ParameterSetDescription desc;

  std::string lepCollName;
  if (typeid(T) == typeid(pat::Electron))
    lepCollName = "slimmedElectrons";
  else if (typeid(T) == typeid(pat::Muon))
    lepCollName = "slimmedMuons";
  else if (typeid(T) == typeid(pat::Tau))
    lepCollName = "slimmedTaus";
  desc.add<edm::InputTag>("src", edm::InputTag(lepCollName));
  desc.add<edm::InputTag>("pvSource", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<int>("pvChoice", useFront)
      ->setComment(
          "Define PV to compute IP: 0: first PV, 1: PV with the smallest dz of the tau leading track (default: 0)");

  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
typedef PATLeptonTimeLifeInfoUpdater<pat::Electron> PATElectronTimeLifeInfoUpdater;
DEFINE_FWK_MODULE(PATElectronTimeLifeInfoUpdater);
typedef PATLeptonTimeLifeInfoUpdater<pat::Muon> PATMuonTimeLifeInfoUpdater;
DEFINE_FWK_MODULE(PATMuonTimeLifeInfoUpdater);
typedef PATLeptonTimeLifeInfoUpdater<pat::Tau> PATTauTimeLifeInfoUpdater;
DEFINE_FWK_MODULE(PATTauTimeLifeInfoUpdater);
