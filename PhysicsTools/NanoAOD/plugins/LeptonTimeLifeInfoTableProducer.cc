/**
  \class    LeptonTimeLifeInfoTableProducer
  \brief    Produces FlatTable with lepton life-time information

  \author   Michal Bluj, NCBJ, Warsaw
*/

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/allowedValues.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
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
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"

#include "Utilities/General/interface/ClassName.h"
#include "PhysicsTools/NanoAOD/interface/DumpedVariable.h"
#include "DataFormats/VertexReco/interface/TrackTimeLifeInfo.h"

#include <cstring>

template <typename T>
class LeptonTimeLifeInfoTableProducer : public edm::stream::EDProducer<> {
public:
  explicit LeptonTimeLifeInfoTableProducer(const edm::ParameterSet&);
  ~LeptonTimeLifeInfoTableProducer() override{};

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
  //--- private utility methods
  const reco::Track* getTrack(const T&);
  void produceAndFillIPInfo(const T&, const TransientTrackBuilder&, const reco::Vertex&, TrackTimeLifeInfo&);
  void produceAndFillSVInfo(const T&, const TransientTrackBuilder&, const reco::Vertex&, TrackTimeLifeInfo&);
  static bool fitVertex(const std::vector<reco::TransientTrack>& transTrk, TransientVertex& transVtx) {
    if (transTrk.size() < 2)
      return false;
    KalmanVertexFitter kvf(true);
    transVtx = kvf.vertex(transTrk);
    return transVtx.hasRefittedTracks() && transVtx.refittedTracks().size() == transTrk.size();
  }

  //--- configuration parameters
  edm::EDGetTokenT<std::vector<T>> leptonsToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transTrackBuilderToken_;
  std::string name_, doc_;
  bool extension_;
  const StringCutObjectSelector<T> selector_;
  int pvChoice_;
  std::vector<std::unique_ptr<Variable<TrackTimeLifeInfo>>> vars_;

  enum PVChoice { useFront = 0, useClosestInDz };
  template <typename ValType>
  using TrackTimeLifeInfoVar = FuncVariable<TrackTimeLifeInfo, StringObjectFunction<TrackTimeLifeInfo>, ValType>;
};

template <typename T>
LeptonTimeLifeInfoTableProducer<T>::LeptonTimeLifeInfoTableProducer(const edm::ParameterSet& cfg)
    : leptonsToken_(consumes<std::vector<T>>(cfg.getParameter<edm::InputTag>("src"))),
      pvToken_(consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("pvSource"))),
      transTrackBuilderToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
      name_(cfg.getParameter<std::string>("name")),
      doc_(cfg.getParameter<std::string>("doc")),
      extension_(cfg.getParameter<bool>("extension")),
      selector_(cfg.getParameter<std::string>("selection")),
      pvChoice_(cfg.getParameter<int>("pvChoice")) {
  //Time-life info variables
  edm::ParameterSet const& varsPSet = cfg.getParameter<edm::ParameterSet>("variables");
  for (const std::string& vname : varsPSet.getParameterNamesForType<edm::ParameterSet>()) {
    const auto& varPSet = varsPSet.getParameter<edm::ParameterSet>(vname);
    const std::string& type = varPSet.getParameter<std::string>("type");
    if (type == "int")
      vars_.push_back(std::make_unique<TrackTimeLifeInfoVar<int32_t>>(vname, varPSet));
    else if (type == "uint")
      vars_.push_back(std::make_unique<TrackTimeLifeInfoVar<uint32_t>>(vname, varPSet));
    else if (type == "float")
      vars_.push_back(std::make_unique<TrackTimeLifeInfoVar<float>>(vname, varPSet));
    else if (type == "double")
      vars_.push_back(std::make_unique<TrackTimeLifeInfoVar<double>>(vname, varPSet));
    else if (type == "int8")
      vars_.push_back(std::make_unique<TrackTimeLifeInfoVar<int8_t>>(vname, varPSet));
    else if (type == "uint8")
      vars_.push_back(std::make_unique<TrackTimeLifeInfoVar<uint8_t>>(vname, varPSet));
    else if (type == "int16")
      vars_.push_back(std::make_unique<TrackTimeLifeInfoVar<int16_t>>(vname, varPSet));
    else if (type == "uint16")
      vars_.push_back(std::make_unique<TrackTimeLifeInfoVar<uint16_t>>(vname, varPSet));
    else if (type == "bool")
      vars_.push_back(
          std::make_unique<FuncVariable<TrackTimeLifeInfo, StringCutObjectSelector<TrackTimeLifeInfo>, bool>>(vname,
                                                                                                              varPSet));
    else
      throw cms::Exception("Configuration", "unsupported type " + type + " for variable " + vname);
  }
  produces<nanoaod::FlatTable>();
}

template <typename T>
void LeptonTimeLifeInfoTableProducer<T>::produce(edm::Event& evt, const edm::EventSetup& es) {
  // Get leptons
  edm::Handle<std::vector<T>> leptons;
  evt.getByToken(leptonsToken_, leptons);

  // Get the vertices
  edm::Handle<reco::VertexCollection> vertices;
  evt.getByToken(pvToken_, vertices);

  // Get transient track builder
  const TransientTrackBuilder& transTrackBuilder = es.getData(transTrackBuilderToken_);

  std::vector<const TrackTimeLifeInfo*> infos;
  infos.reserve(leptons->size());

  for (const auto& lepton : *leptons) {
    TrackTimeLifeInfo* info = new TrackTimeLifeInfo();

    // Do nothing for lepton not passing selection
    if (!selector_(lepton)) {
      infos.push_back(info);
      continue;
    }
    size_t pv_idx = 0;
    if (pvChoice_ == useClosestInDz && getTrack(lepton) != nullptr) {
      float dz_min = 999;
      size_t vtx_idx = 0;
      for (const auto& vtx : *vertices) {
        float dz_tmp = std::abs(getTrack(lepton)->dz(vtx.position()));
        if (dz_tmp < dz_min) {
          dz_min = dz_tmp;
          pv_idx = vtx_idx;
        }
        vtx_idx++;
      }
    }
    const reco::Vertex& pv = !vertices->empty() ? (*vertices)[pv_idx] : reco::Vertex();

    // Obtain IP vector and set related info into lepton
    produceAndFillIPInfo(lepton, transTrackBuilder, pv, *info);

    // Fit SV and set related info for taus or do nothing for other lepton types
    produceAndFillSVInfo(lepton, transTrackBuilder, pv, *info);
    infos.push_back(info);
  }  // end of lepton loop

  // Define and fill table
  auto timelifeTable = std::make_unique<nanoaod::FlatTable>(leptons->size(), name_, false, extension_);
  timelifeTable->setDoc(doc_);
  for (const auto& var : vars_)
    var->fill(infos, *timelifeTable);
  for (auto info : infos)
    delete info;

  // Store table in the event
  evt.put(std::move(timelifeTable));
}

template <>
const reco::Track* LeptonTimeLifeInfoTableProducer<pat::Electron>::getTrack(const pat::Electron& electron) {
  return electron.gsfTrack().isNonnull() ? electron.gsfTrack().get() : nullptr;
}

template <>
const reco::Track* LeptonTimeLifeInfoTableProducer<pat::Muon>::getTrack(const pat::Muon& muon) {
  return muon.innerTrack().isNonnull() ? muon.innerTrack().get() : nullptr;
}

template <>
const reco::Track* LeptonTimeLifeInfoTableProducer<pat::Tau>::getTrack(const pat::Tau& tau) {
  const reco::Track* track = nullptr;
  if (tau.leadChargedHadrCand().isNonnull())
    track = tau.leadChargedHadrCand()->bestTrack();
  return track;
}

template <typename T>
void LeptonTimeLifeInfoTableProducer<T>::produceAndFillIPInfo(const T& lepton,
                                                              const TransientTrackBuilder& transTrackBuilder,
                                                              const reco::Vertex& pv,
                                                              TrackTimeLifeInfo& info) {
  const reco::Track* track = getTrack(lepton);
  if (track != nullptr) {
    info.setTrack(track);
    info.setBFiled_z(transTrackBuilder.field()->inInverseGeV(GlobalPoint(track->vx(), track->vy(), track->vz())).z());

    // Extrapolate track to the point closest to PV
    reco::TransientTrack transTrack = transTrackBuilder.build(track);
    AnalyticalImpactPointExtrapolator extrapolator(transTrack.field());
    TrajectoryStateOnSurface closestState =
        extrapolator.extrapolate(transTrack.impactPointState(), RecoVertex::convertPos(pv.position()));
    GlobalPoint pca = closestState.globalPosition();
    GlobalError pca_cov = closestState.cartesianError().position();
    GlobalVector ip_vec = GlobalVector(pca.x() - pv.x(), pca.y() - pv.y(), pca.z() - pv.z());
    GlobalError ip_cov = pca_cov + GlobalError(pv.covariance());
    VertexDistance3D pca_dist;
    Measurement1D ip_mes = pca_dist.distance(pv, VertexState(pca, pca_cov));
    if (ip_vec.dot(GlobalVector(lepton.px(), lepton.py(), lepton.pz())) < 0)
      ip_mes = Measurement1D(-1. * ip_mes.value(), ip_mes.error());

    // Store PCA info
    info.setPCA(pca, pca_cov);
    info.setIP(ip_vec, ip_cov);
    info.setIPLength(ip_mes);
  }
}

template <typename T>
void LeptonTimeLifeInfoTableProducer<T>::produceAndFillSVInfo(const T& lepton,
                                                              const TransientTrackBuilder& transTrackBuilder,
                                                              const reco::Vertex& pv,
                                                              TrackTimeLifeInfo& info) {}

template <>
void LeptonTimeLifeInfoTableProducer<pat::Tau>::produceAndFillSVInfo(const pat::Tau& tau,
                                                                     const TransientTrackBuilder& transTrackBuilder,
                                                                     const reco::Vertex& pv,
                                                                     TrackTimeLifeInfo& info) {
  // Fit SV with tracks of charged tau decay products
  int fitOK = 0;
  if (tau.signalChargedHadrCands().size() + tau.signalLostTracks().size() > 1) {
    // Get tracks from tau signal charged candidates
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
    fitOK = fitVertex(transTrks, transVtx) ? 1 : -1;
    if (fitOK > 0) {
      reco::Vertex sv = transVtx;
      // Get flight-length
      // Full PV->SV flight vector with its covariance
      GlobalVector flight_vec = GlobalVector(sv.x() - pv.x(), sv.y() - pv.y(), sv.z() - pv.z());
      GlobalError flight_cov = transVtx.positionError() + GlobalError(pv.covariance());
      //MB: can be taken from tau itself (but with different fit of PV and SV) as follows:
      //tau.flightLength().mag2());
      //tau.flightLengthSig();
      VertexDistance3D sv_dist;
      Measurement1D flightLength_mes = sv_dist.signedDistance(pv, sv, GlobalVector(tau.px(), tau.py(), tau.pz()));

      // Store SV info
      info.setSV(sv);
      info.setFlightVector(flight_vec, flight_cov);
      info.setFlightLength(flightLength_mes);
    }
  }
}

template <typename T>
void LeptonTimeLifeInfoTableProducer<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
  desc.add<std::string>("selection", "")->setComment("Selection required to produce and store time-life information");
  desc.add<std::string>("name")->setComment("Name of the branch in the time-life info table for " +
                                            ClassName<T>::name());
  desc.add<std::string>("doc", "")->setComment("Few words of self documentation");
  desc.add<bool>("extension", true)->setComment("Whether or not to extend an existing same table (default: true)");
  desc.add<int>("pvChoice", useFront)
      ->setComment(
          "Define PV to compute IP: 0: first PV, 1: PV with the smallest dz of the tau leading track (default: " +
          std::to_string(useFront) + ")");

  // variables
  edm::ParameterSetDescription variable;
  variable.add<std::string>("expr")->setComment("a function to define the content of the branch in the flat table");
  variable.add<std::string>("doc")->setComment("few words description of the branch content");
  variable.ifValue(
      edm::ParameterDescription<std::string>(
          "type", "int", true, edm::Comment("the c++ type of the branch in the flat table")),
      edm::allowedValues<std::string>("int", "uint", "float", "double", "int8", "uint8", "int16", "uint16", "bool"));
  variable.addOptionalNode(
      edm::ParameterDescription<int>(
          "precision", true, edm::Comment("the precision with which to store the value in the flat table")) xor
          edm::ParameterDescription<std::string>(
              "precision", true, edm::Comment("the precision with which to store the value in the flat table")),
      false);
  edm::ParameterSetDescription variables;
  variables.setComment("a parameters set to define additional variables describing main vertex");
  variables.addNode(edm::ParameterWildcard<edm::ParameterSetDescription>("*", edm::RequireZeroOrMore, true, variable));
  desc.add<edm::ParameterSetDescription>("variables", variables)->setComment("Time-life info variables");

  descriptions.addWithDefaultLabel(desc);
}

#include "FWCore/Framework/interface/MakerMacros.h"
typedef LeptonTimeLifeInfoTableProducer<pat::Electron> ElectronTimeLifeInfoTableProducer;
DEFINE_FWK_MODULE(ElectronTimeLifeInfoTableProducer);
typedef LeptonTimeLifeInfoTableProducer<pat::Muon> MuonTimeLifeInfoTableProducer;
DEFINE_FWK_MODULE(MuonTimeLifeInfoTableProducer);
typedef LeptonTimeLifeInfoTableProducer<pat::Tau> TauTimeLifeInfoTableProducer;
DEFINE_FWK_MODULE(TauTimeLifeInfoTableProducer);
