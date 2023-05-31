#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFSorter.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GoldenPattern.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GoldenPatternWithStat.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/OMTFConfiguration.h"

#include <cassert>
#include <iostream>
#include <strstream>
#include <algorithm>
#include <bitset>

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/RPCTrigger/interface/RPCConst.h"
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
template <class GoldenPatternType>
AlgoMuons::value_type OMTFSorter<GoldenPatternType>::sortRefHitResults(
    unsigned int procIndx, unsigned int iRefHit, const GoldenPatternVec<GoldenPatternType>& gPatterns, int charge) {
  GoldenPatternType* bestGP = nullptr;  //the GoldenPattern with the best result for this iRefHit
  //std::cout <<" ====== sortRefHitResults: " << std::endl;

  GoldenPatternType* bestGpUpt = nullptr;

  for (auto& itGP : gPatterns) {
    if (!itGP->getResults()[procIndx][iRefHit].isValid())
      continue;

    if (charge != 0 && itGP->key().theCharge != charge)
      continue;  //charge==0 means ignore charge

    ///Accept only candidates with >2 hits
    if (itGP->getResults()[procIndx][iRefHit].getFiredLayerCnt() < 3)  //TODO - move 3 to the configuration??
      continue;

    if (bestGP == nullptr) {
      bestGP = itGP.get();
    } else if (myType == 0 && itGP->getResults()[procIndx][iRefHit].getFiredLayerCnt() >
                                  bestGP->getResults()[procIndx][iRefHit].getFiredLayerCnt()) {
      bestGP = itGP.get();
    } else if (myType == 1 || (itGP->getResults()[procIndx][iRefHit].getFiredLayerCnt() ==
                               bestGP->getResults()[procIndx][iRefHit].getFiredLayerCnt())) {
      if (itGP->getResults()[procIndx][iRefHit].getPdfSum() > bestGP->getResults()[procIndx][iRefHit].getPdfSum()) {
        //if the PdfWeigtSum is equal, we take the GP with the lower number, i.e. lower pt = check if this is ok for physics FIXME (KB)
        bestGP = itGP.get();
      }
    }


    if (bestGpUpt == nullptr) {
      if(itGP->getResults()[procIndx][iRefHit].getPdfSumUpt() > 0)
        bestGpUpt = itGP.get();
    } else if (myType == 0 && itGP->getResults()[procIndx][iRefHit].getFiredLayerCnt() >
                              bestGpUpt->getResults()[procIndx][iRefHit].getFiredLayerCnt()) {
      bestGpUpt = itGP.get();
    } else if (myType == 1 || (itGP->getResults()[procIndx][iRefHit].getFiredLayerCnt() ==
                               bestGpUpt->getResults()[procIndx][iRefHit].getFiredLayerCnt())) {
      if (itGP->getResults()[procIndx][iRefHit].getPdfSumUpt() > bestGpUpt->getResults()[procIndx][iRefHit].getPdfSumUpt()) {
        //if the PdfWeigtSum is equal, we take the GP with the lower number, i.e. lower pt = check if this is ok for physics FIXME (KB)
        bestGpUpt = itGP.get();
      }
    }

//    LogTrace("OMTFReconstruction")<<"OMTFSorter::sortRefHitResults() " << itGP->key()<<" getFiredLayer "<< itGP->getResults()[procIndx][iRefHit].getFiredLayerCnt()
//    <<" PdfSum "<< itGP->getResults()[procIndx][iRefHit].getPdfSum()
//    <<" PdfSumUpt "<< itGP->getResults()[procIndx][iRefHit].getPdfSumUpt()<< std::endl;
  }
  if (bestGP) {
    AlgoMuons::value_type candidate(new AlgoMuon(bestGP->getResults()[procIndx][iRefHit], bestGP, iRefHit));

    if(bestGpUpt) {
      candidate->setGpResultUpt(bestGpUpt->getResults()[procIndx][iRefHit]);
      candidate->setGoldenPaternUpt(bestGpUpt);
    }

    //std::cout<<__FUNCTION__<<" line "<<__LINE__ <<" return: " << *candidate << std::endl;
    return candidate;
  } else {
    AlgoMuons::value_type candidate(new AlgoMuon());
    candidate->setRefHitNumber(iRefHit);
    return candidate;
  }
}

template class OMTFSorter<GoldenPattern>;
template class OMTFSorter<GoldenPatternWithStat>;
template class OMTFSorter<GoldenPatternWithThresh>;
