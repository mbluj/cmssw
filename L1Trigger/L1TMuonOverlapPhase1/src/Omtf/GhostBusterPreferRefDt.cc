#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GhostBusterPreferRefDt.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFConfiguration.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <sstream>

namespace {

  struct AlgoMuonEtaFix {
    AlgoMuonEtaFix(const AlgoMuonPtr& mu) : mu(mu), fixedEta(mu->getEtaHw()) {}
    const AlgoMuonPtr mu;
    unsigned int fixedEta;
    const AlgoMuon* operator -> () {
      return mu.get();
    }
  };

}  // namespace

AlgoMuons GhostBusterPreferRefDt::select(AlgoMuons muonsIN, int charge) {
  // sorting within GB.
  auto customLess = [&](const AlgoMuons::value_type& a, const AlgoMuons::value_type& b) -> bool {
    if (!a->isValid()) {
      return true;
    }
    if (!b->isValid()) {
      return false;
    }

    int aRefLayerLogicNum = omtfConfig->getRefToLogicNumber()[a->getRefLayer()];
    int bRefLayerLogicNum = omtfConfig->getRefToLogicNumber()[b->getRefLayer()];
    if (a->getQ() > b->getQ())
      return false;
    else if (a->getQ() == b->getQ() && aRefLayerLogicNum < bRefLayerLogicNum) {
      return false;
    } else if (a->getQ() == b->getQ() && aRefLayerLogicNum == bRefLayerLogicNum && a->getDisc() > b->getDisc())
      return false;
    else if (a->getQ() == b->getQ() && aRefLayerLogicNum == bRefLayerLogicNum && a->getDisc() == b->getDisc() &&
             a->getPatternNumber() > b->getPatternNumber())
      return false;
    else if (a->getQ() == b->getQ() && aRefLayerLogicNum == bRefLayerLogicNum && a->getDisc() == b->getDisc() &&
             a->getPatternNumber() == b->getPatternNumber() && a->getRefHitNumber() < b->getRefHitNumber())
      return false;
    else
      return true;
  };

  auto customLessByFPLLH = [&](const AlgoMuons::value_type& a, const AlgoMuons::value_type& b) -> bool {
    if (!a->isValid()) {
      return true;
    }
    if (!b->isValid()) {
      return false;
    }

    if (a->getQ() > b->getQ())
      return false;
    else if (a->getQ() == b->getQ()) {
      return false;
    } else if (a->getQ() == b->getQ() && a->getDisc() > b->getDisc())
      return false;
    else if (a->getQ() == b->getQ() && a->getDisc() == b->getDisc() && a->getPatternNumber() > b->getPatternNumber())
      return false;
    else if (a->getQ() == b->getQ() && a->getDisc() == b->getDisc() && a->getPatternNumber() == b->getPatternNumber() &&
             a->getRefHitNumber() < b->getRefHitNumber())
      return false;
    else
      return true;
  };

  auto customLessByLLH = [&](const AlgoMuons::value_type& a, const AlgoMuons::value_type& b) -> bool {
    if (!a->isValid()) {
      return true;
    }
    if (!b->isValid()) {
      return false;
    }

    if (a->getDisc() > b->getDisc())
      return false;
    else if (a->getDisc() == b->getDisc() && a->getPatternNumber() > b->getPatternNumber())
      return false;
    else if (a->getDisc() == b->getDisc() && a->getPatternNumber() == b->getPatternNumber() &&
             a->getRefHitNumber() < b->getRefHitNumber())
      return false;
    else
      return true;
  };

  auto customByRefLayer = [&](const AlgoMuons::value_type& a, const AlgoMuons::value_type& b)->bool {
    if(!a->isValid()) {
      return true;
    }
    if (!b->isValid()) {
      return false;
    }

    int aRefLayerLogicNum = omtfConfig->getRefToLogicNumber()[a->getRefLayer()];
    int bRefLayerLogicNum = omtfConfig->getRefToLogicNumber()[b->getRefLayer()];

    if (aRefLayerLogicNum < bRefLayerLogicNum) {
      return false;
    }
    // if(a->getQ() > b->getQ())
    //   return false;
    else if (aRefLayerLogicNum == bRefLayerLogicNum && a->getDisc() > b->getDisc()) //TODO how about getPdfSumUpt ????
      return false;
    else if (aRefLayerLogicNum == bRefLayerLogicNum && a->getDisc() == b->getDisc() &&
             a->getPatternNumber() > b->getPatternNumber())
      return false;
    else if (aRefLayerLogicNum == bRefLayerLogicNum && a->getDisc() == b->getDisc() &&
             a->getPatternNumber() == b->getPatternNumber() && a->getRefHitNumber() < b->getRefHitNumber())
      return false;
    else
      return true;
  };

  if (omtfConfig->getGhostBusterType() == "byLLH")
    std::sort(muonsIN.rbegin(), muonsIN.rend(), customLessByLLH);
  else if (omtfConfig->getGhostBusterType() == "byFPLLH")
    std::sort(muonsIN.rbegin(), muonsIN.rend(), customLessByFPLLH);
  else if (omtfConfig->getGhostBusterType() == "byRefLayer")
    std::sort(muonsIN.rbegin(), muonsIN.rend(), customByRefLayer);
  else
    std::sort(muonsIN.rbegin(), muonsIN.rend(), customLess);

  // actual GhostBusting. Overwrite eta in case of no DT info.
  std::vector<AlgoMuonEtaFix> refHitCleanCandsFixedEta;

  for (unsigned int iMu1 = 0; iMu1 < muonsIN.size(); iMu1++) {
    auto& muIN1 = muonsIN[iMu1];
    if (!muIN1->isValid() || muIN1->isKilled())
      continue;

    refHitCleanCandsFixedEta.push_back(muIN1);
    for (unsigned int iMu2 = iMu1+1; iMu2 < muonsIN.size(); iMu2++) {
      auto& muIN2 = muonsIN[iMu2];
      if (muIN2->isValid() &&
          std::abs(omtfConfig->procPhiToGmtPhi(muIN1->getPhi()) - omtfConfig->procPhiToGmtPhi(muIN2->getPhi())) < 8) {
        //the candidates are sorted, so only the  muIN2 can be killed, as it is "worse" than the muIN1
        muonsIN[iMu2]->kill();
        muonsIN[iMu1]->getKilledMuons().emplace_back(muIN2);

        if ((omtfConfig->fwVersion() >= 6) &&
            ((abs(muIN1->getEtaHw()) == 75 || abs(muIN1->getEtaHw()) == 79 || abs(muIN1->getEtaHw()) == 92)) &&
            ((abs(muIN2->getEtaHw()) != 75 && abs(muIN2->getEtaHw()) != 79 && abs(muIN2->getEtaHw()) != 92))) {

          refHitCleanCandsFixedEta.back().fixedEta = muIN2->getEtaHw();

        }
      }
    }
  }

  // fill outgoing collection
  AlgoMuons refHitCleanCands;
  for (const auto& mu : refHitCleanCandsFixedEta) {
    refHitCleanCands.emplace_back(new AlgoMuon( *(mu.mu) ));
    refHitCleanCands.back()->setEta(mu.fixedEta);
    if (refHitCleanCands.size() >= 3)
      break;
  }

  while (refHitCleanCands.size() < 3)
    refHitCleanCands.emplace_back(new AlgoMuon());

  return refHitCleanCands;
}
