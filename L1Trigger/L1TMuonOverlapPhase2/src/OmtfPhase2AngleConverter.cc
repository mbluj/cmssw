#include "L1Trigger/L1TMuonOverlapPhase2/interface/OmtfPhase2AngleConverter.h"

int OmtfPhase2AngleConverter::getProcessorPhi(int phiZero, l1t::tftype part, int dtScNum, int dtPhi) const {
  int dtPhiBins = 65536;                    //65536. per 0.8 radians
  double hsPhiPitch = 2 * M_PI / nPhiBins;  // width of phi Pitch, related to halfStrip at CSC station 2

  int sector = dtScNum + 1;  //NOTE: there is a inconsistency in DT sector numb. Thus +1 needed to get detector numb.

  double scale = 0.8 / dtPhiBins / hsPhiPitch;
  int scale_coeff = lround(scale * pow(2, 11));

  int ichamber = sector - 1;
  if (ichamber > 6)
    ichamber = ichamber - 12;

  int offsetGlobal = (int)nPhiBins * ichamber / 12;

  int phiConverted = floor(dtPhi * scale_coeff / pow(2, 11)) + offsetGlobal - phiZero;

  return config->foldPhi(phiConverted);
}
