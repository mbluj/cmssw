/*
 * DataROOTDumper2.cc
 *
 *  Created on: Dec 11, 2019
 *      Author: kbunkow
 */

#include "L1Trigger/L1TMuonOverlapPhase1/interface/Tools/DataROOTDumper2.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "TObjString.h"

/*
#include <boost/range/adaptor/reversed.hpp>
#include <boost/timer/timer.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
*/

DataROOTDumper2::DataROOTDumper2(const edm::ParameterSet& edmCfg,
    const OMTFConfiguration* omtfConfig,
    CandidateSimMuonMatcher* candidateSimMuonMatcher)
: EmulationObserverBase(edmCfg, omtfConfig), candidateSimMuonMatcher(candidateSimMuonMatcher) {
  edm::LogVerbatim("l1tOmtfEventPrint") << " omtfConfig->nTestRefHits() " << omtfConfig->nTestRefHits()
                                            << " event.omtfGpResultsPdfSum.num_elements() " << endl;
  initializeTTree("dump.root"); //TODO

  if (edmCfg.exists("dumpKilledOmtfCands"))
    if (edmCfg.getParameter<bool>("dumpKilledOmtfCands"))
      dumpKilledOmtfCands = true;

  edm::LogVerbatim("l1tOmtfEventPrint") << " DataROOTDumper2 created. dumpKilledOmtfCands " <<dumpKilledOmtfCands<< std::endl;
}

DataROOTDumper2::~DataROOTDumper2() { saveTTree(); }

void DataROOTDumper2::initializeTTree(std::string rootFileName) {
  edm::Service<TFileService> fs;

  //TFileDirectory subDir = fs->mkdir("OmtfDataDumper");

  rootTree = fs->make<TTree>("OMTFHitsTree", "");

  rootTree->Branch("eventNum", &omtfEvent.eventNum);
  rootTree->Branch("muonEvent", &omtfEvent.muonEvent);

  rootTree->Branch("muonPt", &omtfEvent.muonPt);
  rootTree->Branch("muonEta", &omtfEvent.muonEta);
  rootTree->Branch("muonPhi", &omtfEvent.muonPhi);
  rootTree->Branch("muonCharge", &omtfEvent.muonCharge);

  rootTree->Branch("muonDxy", &omtfEvent.muonDxy);
  rootTree->Branch("muonRho", &omtfEvent.muonRho);

  rootTree->Branch("omtfPt", &omtfEvent.omtfPt);
  rootTree->Branch("omtfEta", &omtfEvent.omtfEta);
  rootTree->Branch("omtfPhi", &omtfEvent.omtfPhi);
  rootTree->Branch("omtfCharge", &omtfEvent.omtfCharge);

  rootTree->Branch("omtfHwEta", &omtfEvent.omtfHwEta);

  rootTree->Branch("omtfProcessor", &omtfEvent.omtfProcessor);
  rootTree->Branch("omtfScore", &omtfEvent.omtfScore);
  rootTree->Branch("omtfQuality", &omtfEvent.omtfQuality);
  rootTree->Branch("omtfRefLayer", &omtfEvent.omtfRefLayer);
  rootTree->Branch("omtfRefHitNum", &omtfEvent.omtfRefHitNum);

  rootTree->Branch("omtfFiredLayers", &omtfEvent.omtfFiredLayers);  //<<<<<<<<<<<<<<<<<<<<<<!!!!TODOO

  rootTree->Branch("killed", &omtfEvent.killed);
  
  rootTree->Branch("hits", &omtfEvent.hits);

  ptGenPos = fs->make<TH1I>("ptGenPos", "ptGenPos, eta at vertex 0.8 - 1.24", 400, 0, 200); //TODO
  ptGenNeg = fs->make<TH1I>("ptGenNeg", "ptGenNeg, eta at vertex 0.8 - 1.24", 400, 0, 200);
}

void DataROOTDumper2::saveTTree() {

}

void DataROOTDumper2::observeProcesorEmulation(unsigned int iProcessor,
    l1t::tftype mtfType,
    const std::shared_ptr<OMTFinput>&,
    const AlgoMuons& algoCandidates,
    const AlgoMuons& gbCandidates,
    const std::vector<l1t::RegionalMuonCand>& candMuons) {

}

void DataROOTDumper2::observeEventEnd(const edm::Event& iEvent,
    std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates) {
  /*
  int muonCharge = 0;
  if (simMuon) {
    if (fabs(simMuon->momentum().eta()) < 0.8 || fabs(simMuon->momentum().eta()) > 1.24)
      return;

    muonCharge = (abs(simMuon->type()) == 13) ? simMuon->type() / -13 : 0;
    if (muonCharge > 0)
      ptGenPos->Fill(simMuon->momentum().pt());
    else
      ptGenNeg->Fill(simMuon->momentum().pt());
  }

  if (simMuon == nullptr || !omtfCand->isValid())  //no sim muon or empty candidate
    return;

  omtfEvent.muonPt = simMuon->momentum().pt();
  omtfEvent.muonEta = simMuon->momentum().eta();

  //TODO add cut on ete if needed
    if(fabs(event.muonEta) < 0.8 || fabs(event.muonEta) > 1.24)
    return;

  omtfEvent.muonPhi = simMuon->momentum().phi();
  omtfEvent.muonCharge = muonCharge;  //TODO
   */

  std::vector<MatchingResult> matchingResults = candidateSimMuonMatcher->getMatchingResults();
  LogTrace("l1tOmtfEventPrint") << "\nDataROOTDumper2::observeEventEnd matchingResults.size() " << matchingResults.size() << std::endl;

  //candidateSimMuonMatcher should use the  trackingParticles, because the simTracks are not stored for the pile-up events
  for (auto& matchingResult : matchingResults) {
    omtfEvent.eventNum = iEvent.id().event();

    if (matchingResult.trackingParticle) {
      auto trackingParticle = matchingResult.trackingParticle;

      omtfEvent.muonEvent = trackingParticle->eventId().event();

      omtfEvent.muonPt = trackingParticle->pt();
      omtfEvent.muonEta = trackingParticle->momentum().eta();

      omtfEvent.muonPhi = trackingParticle->momentum().phi();
      omtfEvent.muonCharge = (abs(trackingParticle->pdgId()) == 13) ? trackingParticle->pdgId() / -13 : 0;;  //TODO

      //omtfEvent.muonDxy = TODO
      if(trackingParticle->parentVertex().isNonnull())
        omtfEvent.muonRho = trackingParticle->parentVertex()->position().Rho();


      LogTrace("l1tOmtfEventPrint") << "DataROOTDumper2::observeEventEnd trackingParticle: eventId " << trackingParticle->eventId().event() << " pdgId " << std::setw(3)
                 << trackingParticle->pdgId() << " trackId " << trackingParticle->g4Tracks().at(0).trackId() << " pt "
                 << std::setw(9) << trackingParticle->pt()  //<<" Beta "<<simMuon->momentum().Beta()
                 << " eta " << std::setw(9) << trackingParticle->momentum().eta() << " phi " << std::setw(9)
                 << trackingParticle->momentum().phi() << std::endl;

      if(fabs(omtfEvent.muonEta) > 0.8 && fabs(omtfEvent.muonEta) < 1.24) {
        if (omtfEvent.muonCharge > 0)
          ptGenPos->Fill(omtfEvent.muonPt);
        else
          ptGenNeg->Fill(omtfEvent.muonPt);
      }
    }
    else if (matchingResult.simTrack) {
      auto simTrack = matchingResult.simTrack;

      omtfEvent.muonEvent = simTrack->eventId().event();

      omtfEvent.muonPt = simTrack->momentum().pt();
      omtfEvent.muonEta = simTrack->momentum().eta();

      omtfEvent.muonPhi = simTrack->momentum().phi();
      omtfEvent.muonCharge = (abs(simTrack->type()) == 13) ? simTrack->type() / -13 : 0;;  //TODO

      //omtfEvent.muonDxy = TODO
      //if(trackingParticle->parentVertex().isNonnull()) //TODO!!!!!!!!!!!!!!!!!!
      //  omtfEvent.muonRho = trackingParticle->parentVertex()->position().Rho();

      LogTrace("l1tOmtfEventPrint") << "DataROOTDumper2::observeEventEnd trackingParticle: eventId " << simTrack->eventId().event() << " pdgId " << std::setw(3)
                 << simTrack->type() //<< " trackId " << simTrack->g4Tracks().at(0).trackId()
                 << " pt "<< std::setw(9) << simTrack->momentum().pt()  //<<" Beta "<<simMuon->momentum().Beta()
                 << " eta " << std::setw(9) << simTrack->momentum().eta() << " phi " << std::setw(9)
                 << simTrack->momentum().phi() << std::endl;

      if(fabs(omtfEvent.muonEta) > 0.8 && fabs(omtfEvent.muonEta) <1.24) {
        if (omtfEvent.muonCharge > 0)
          ptGenPos->Fill(omtfEvent.muonPt);
        else
          ptGenNeg->Fill(omtfEvent.muonPt);
      }
    }
    else {
      omtfEvent.muonEvent = -1;

      omtfEvent.muonPt = 0;
      omtfEvent.muonEta = 0;

      omtfEvent.muonPhi = 0;
      omtfEvent.muonCharge = 0;  //TODO

      omtfEvent.muonDxy = 0;
      omtfEvent.muonRho = 0;
    }


    auto addOmtfCand = [&](AlgoMuonPtr& procMuon) {
      omtfEvent.omtfPt = omtfConfig->hwPtToGev(procMuon->getPt());
      omtfEvent.omtfEta = omtfConfig->hwEtaToEta(procMuon->getEtaHw());
      omtfEvent.omtfPhi = procMuon->getPhi();
      omtfEvent.omtfCharge = procMuon->getCharge();
      omtfEvent.omtfScore = procMuon->getPdfSum();

      omtfEvent.omtfHwEta = procMuon->getEtaHw();

      omtfEvent.omtfFiredLayers = procMuon->getFiredLayerBits();
      omtfEvent.omtfRefLayer = procMuon->getRefLayer();
      omtfEvent.omtfRefHitNum = procMuon->getRefHitNumber();

      omtfEvent.hits.clear();

      auto& gpResult = procMuon->getGpResult();
      /*
        edm::LogVerbatim("l1tOmtfEventPrint")<<"DataROOTDumper2:;observeEventEnd muonPt "<<event.muonPt<<" muonCharge "<<event.muonCharge
            <<" omtfPt "<<event.omtfPt<<" RefLayer "<<event.omtfRefLayer<<" omtfPtCont "<<event.omtfPtCont
            <<std::endl;  */

      for (unsigned int iLogicLayer = 0; iLogicLayer < gpResult.getStubResults().size(); ++iLogicLayer) {
        auto& stubResult = gpResult.getStubResults()[iLogicLayer];
        if (stubResult.getMuonStub()) {  //&& stubResult.getValid() //TODO!!!!!!!!!!!!!!!!1
          OmtfEvent::Hit hit;
          hit.layer = iLogicLayer;
          hit.quality = stubResult.getMuonStub()->qualityHw;
          hit.eta = stubResult.getMuonStub()->etaHw;  //in which scale?
          hit.valid = stubResult.getValid();

          int hitPhi = stubResult.getMuonStub()->phiHw;
          unsigned int refLayerLogicNum = omtfConfig->getRefToLogicNumber()[procMuon->getRefLayer()];
          int phiRefHit = gpResult.getStubResults()[refLayerLogicNum].getMuonStub()->phiHw;

          if (omtfConfig->isBendingLayer(iLogicLayer)) {
            hitPhi = stubResult.getMuonStub()->phiBHw;
            phiRefHit = 0;  //phi ref hit for the bending layer set to 0, since it should not be included in the phiDist
          }

          //phiDist = hitPhi - phiRefHit;
          hit.phiDist = hitPhi - phiRefHit;

          /* LogTrace("l1tOmtfEventPrint")<<" muonPt "<<event.muonPt<<" omtfPt "<<event.omtfPt<<" RefLayer "<<event.omtfRefLayer
                <<" layer "<<int(hit.layer)<<" PdfBin "<<stubResult.getPdfBin()<<" hit.phiDist "<<hit.phiDist<<" valid "<<stubResult.getValid()<<" " //<<" phiDist "<<phiDist
                <<" getDistPhiBitShift "<<procMuon->getGoldenPatern()->getDistPhiBitShift(iLogicLayer, procMuon->getRefLayer())
                <<" meanDistPhiValue   "<<procMuon->getGoldenPatern()->meanDistPhiValue(iLogicLayer, procMuon->getRefLayer())//<<(phiDist != hit.phiDist? "!!!!!!!<<<<<" : "")
                <<endl;*/

          if (hit.phiDist > 504 || hit.phiDist < -512) {
            edm::LogVerbatim("l1tOmtfEventPrint")
                      << " muonPt "<< omtfEvent.muonPt<< " omtfPt " << omtfEvent.omtfPt << " RefLayer "
                      <<(int)omtfEvent.omtfRefLayer << " layer " << int(hit.layer) << " hit.phiDist " << hit.phiDist << " valid "
                      << stubResult.getValid() << " !!!!!!!!!!!!!!!!!!!!!!!!" << endl;
          }

          DetId detId(stubResult.getMuonStub()->detId);
          if (detId.subdetId() == MuonSubdetId::CSC) {
            CSCDetId cscId(detId);
            hit.z = cscId.chamber() % 2;
          }

          omtfEvent.hits.push_back(hit.rawData);
        }
      }

      LogTrace("l1tOmtfEventPrint") << "DataROOTDumper2::observeEventEnd adding omtfCand : " << std::endl;
      auto finalCandidate = matchingResult.muonCand;
      LogTrace("l1tOmtfEventPrint") << " hwPt " << finalCandidate->hwPt() << " hwSign " << finalCandidate->hwSign() << " hwQual "
          << finalCandidate->hwQual() << " hwEta " << std::setw(4) << finalCandidate->hwEta() << std::setw(4)
          << " hwPhi " << finalCandidate->hwPhi() << "    eta " << std::setw(9)
          << (finalCandidate->hwEta() * 0.010875) << " isKilled " << procMuon->isKilled()
          <<" tRefLayer "<<procMuon->getRefLayer()<<" RefHitNumber "<<procMuon->getRefHitNumber()<<std::endl;


    };

    if (matchingResult.muonCand && matchingResult.procMuon->getPt() >= 0 && matchingResult.muonCand->hwQual() >= 1) //TODO set the quality
    {  //&& matchingResult.genPt < 20

      omtfEvent.omtfQuality = matchingResult.muonCand->hwQual();  //procMuon->getQ();
      omtfEvent.killed = false;
      omtfEvent.omtfProcessor = matchingResult.muonCand->processor();

      if( matchingResult.muonCand->trackFinderType() == l1t::omtf_neg) {
        omtfEvent.omtfProcessor *= -1;
      }

      addOmtfCand(matchingResult.procMuon);
      rootTree->Fill();

      if(dumpKilledOmtfCands) {
        for(auto& killedCand : matchingResult.procMuon->getKilledMuons()) {
          omtfEvent.omtfQuality = 0;
          omtfEvent.killed =  true;
          if(killedCand->isKilled() == false) {
            edm::LogVerbatim("l1tOmtfEventPrint")<<" killedCand->isKilled() == false !!!!!!!!";
          }
          addOmtfCand(killedCand);
          rootTree->Fill();
        }
      }
    }
    else {
      LogTrace("l1tOmtfEventPrint") << "DataROOTDumper2::observeEventEnd no matching omtfCand" << std::endl;

      omtfEvent.omtfPt = 0;
      omtfEvent.omtfEta = 0;
      omtfEvent.omtfPhi = 0;
      omtfEvent.omtfCharge = 0;
      omtfEvent.omtfScore = 0;

      omtfEvent.omtfHwEta = 0;

      omtfEvent.omtfFiredLayers = 0;
      omtfEvent.omtfRefLayer = 0;
      omtfEvent.omtfRefHitNum = 0;
      omtfEvent.omtfProcessor = 10;

      omtfEvent.omtfQuality = 0;
      omtfEvent.killed = false;

      omtfEvent.hits.clear();

      rootTree->Fill();
    }
  }
  evntCnt++;
}

void DataROOTDumper2::endJob() { edm::LogVerbatim("l1tOmtfEventPrint") << " evntCnt " << evntCnt << endl; }
