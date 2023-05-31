/*
 * XMLEventWriter.cc
 *
 *  Created on: Oct 12, 2017
 *      Author: kbunkow
 */

#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/Common/interface/EventBase.h"
#include "FWCore/Framework/interface/Event.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OmtfName.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFinput.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/XMLEventWriter.h"

#include <boost/property_tree/xml_parser.hpp>

#include <bitset>

XMLEventWriter::XMLEventWriter(const OMTFConfiguration* aOMTFConfig, std::string fName)
    : omtfConfig(aOMTFConfig), fName(fName) {
  //std::string fName = "OMTF";
  eventNum = 0;

  unsigned int version = aOMTFConfig->patternsVersion();
  unsigned int mask16bits = 0xFFFF;

  version &= mask16bits;

  std::ostringstream stringStr;
  stringStr.str("");
  stringStr << "0x" << std::hex << std::setfill('0') << std::setw(4) << version;

  tree.put("OMTF.<xmlattr>.version", stringStr.str());
};

XMLEventWriter::~XMLEventWriter() {}

void XMLEventWriter::observeProcesorBegin(unsigned int iProcessor, l1t::tftype mtfType) {
  if (eventNum > 5000)
    return;

  procTree.clear();

  int endcap = (mtfType == l1t::omtf_neg) ? -1 : ((mtfType == l1t::omtf_pos) ? +1 : 0);
  OmtfName board(iProcessor, endcap);
  procTree.add("<xmlattr>.board", board.name());
  procTree.add("<xmlattr>.iProcessor", iProcessor);

  std::ostringstream stringStr;
  stringStr << (board.position() == 1 ? "+" : "")<< board.position();
  procTree.add("<xmlattr>.position", stringStr.str());
}

void XMLEventWriter::observeProcesorEmulation(unsigned int iProcessor,
                                              l1t::tftype mtfType,
                                              const std::shared_ptr<OMTFinput>& input,
                                              const AlgoMuons& algoCandidates,
                                              const AlgoMuons& gbCandidates,
                                              const std::vector<l1t::RegionalMuonCand>& candMuons) {
  if (eventNum > 5000)
    return;

  int endcap = (mtfType == l1t::omtf_neg) ? -1 : ((mtfType == l1t::omtf_pos) ? +1 : 0);
  OmtfName board(iProcessor, endcap);

  if (candMuons.empty())
    return;

  for (unsigned int iLayer = 0; iLayer < omtfConfig->nLayers(); ++iLayer) {
    boost::property_tree::ptree layerTree;

    for (unsigned int iHit = 0; iHit < input->getMuonStubs()[iLayer].size(); ++iHit) {
      int hitPhi = input->getPhiHw(iLayer, iHit);
      if (hitPhi >= (int)omtfConfig->nPhiBins())
        continue;

      auto& hitTree = layerTree.add("Hit", "");

      hitTree.add("<xmlattr>.iEta", OMTFConfiguration::eta2Bits(abs(input->getHitEta(iLayer, iHit))));
      hitTree.add("<xmlattr>.iInput", iHit);
      hitTree.add("<xmlattr>.iPhi", hitPhi);
    }

    if(layerTree.size()) {
      layerTree.add("<xmlattr>.iLayer", iLayer);
      procTree.add_child("Layer", layerTree);
    }
  }

  for (auto& algoCand : algoCandidates) {
    ///Dump only regions, where a candidate was found
    if (algoCand->isValid()) {
      auto& algoMuonTree = procTree.add("AlgoMuon", "");
      algoMuonTree.add("<xmlattr>.charge", algoCand->getCharge());
      algoMuonTree.add("<xmlattr>.disc", algoCand->getDisc());
      algoMuonTree.add("<xmlattr>.pdfSum", algoCand->getGpResult().getPdfSum());
      algoMuonTree.add("<xmlattr>.pdfSumUpt", algoCand->getGpResultUpt().getPdfSumUpt());
      algoMuonTree.add("<xmlattr>.etaCode", OMTFConfiguration::eta2Bits(abs(algoCand->getEtaHw())));
      algoMuonTree.add("<xmlattr>.iRefHit", algoCand->getRefHitNumber());
      algoMuonTree.add("<xmlattr>.iRefLayer", algoCand->getRefLayer());
      algoMuonTree.add("<xmlattr>.layers", std::bitset<18>(algoCand->getFiredLayerBits()));
      algoMuonTree.add("<xmlattr>.nHits", algoCand->getQ());
      algoMuonTree.add("<xmlattr>.patNum", algoCand->getHwPatternNumber());
      algoMuonTree.add("<xmlattr>.phiCode", algoCand->getPhi());
      algoMuonTree.add("<xmlattr>.phiRHit", algoCand->getPhiRHit());
      algoMuonTree.add("<xmlattr>.ptCode", algoCand->getPt());
    }
  }

  for (auto& candMuon : candMuons) {
    auto& candMuonTree = procTree.add("CandMuon", "");
    candMuonTree.add("<xmlattr>.hwEta", candMuon.hwEta());
    candMuonTree.add("<xmlattr>.hwPhi", candMuon.hwPhi());
    candMuonTree.add("<xmlattr>.hwPt", candMuon.hwPt());
    candMuonTree.add("<xmlattr>.hwQual", candMuon.hwQual());
    candMuonTree.add("<xmlattr>.hwSign", candMuon.hwSign());
    candMuonTree.add("<xmlattr>.hwSignValid", candMuon.hwSignValid());
    candMuonTree.add("<xmlattr>.hwTrackAddress", std::bitset<29>(candMuon.trackAddress().at(0)));
    candMuonTree.add("<xmlattr>.link", candMuon.link());
    candMuonTree.add("<xmlattr>.processor", candMuon.processor());

    std::ostringstream stringStr;
    if (candMuon.trackFinderType() == l1t::omtf_neg)
      stringStr << "OMTF_NEG";
    else if (candMuon.trackFinderType() == l1t::omtf_pos)
      stringStr << "OMTF_POS";
    else
      stringStr << candMuon.trackFinderType();
    candMuonTree.add("<xmlattr>.trackFinderType", stringStr.str());
  }

  if(procTree.size())
    eventTree->add_child("Processor", procTree);
}

void XMLEventWriter::observeEventBegin(const edm::Event& iEvent) {
  eventNum++;
  if (eventNum > 5000)
    //due to some bug if more events is written the memory consumption s very big and program crashes
    return;
  //currentElement = xmlWriter.writeEventHeader(iEvent.id().event());
  eventId = iEvent.id().event();

  eventTree = &(tree.add("OMTF.Event", ""));
  eventTree->add("<xmlattr>.iEvent", eventId);

  eventTree = &(eventTree->add("bx", ""));
  eventTree->add("<xmlattr>.iBx", 2*eventId);
}

void XMLEventWriter::observeEventEnd(const edm::Event& iEvent,
                                     std::unique_ptr<l1t::RegionalMuonCandBxCollection>& finalCandidates) {
}

void XMLEventWriter::endJob() {
  boost::property_tree::write_xml(fName, tree, std::locale(), boost::property_tree::xml_parser::xml_writer_make_settings<std::string>(' ', 2));
}
