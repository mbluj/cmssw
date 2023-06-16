#ifndef L1T_OmtfP1_AlgoMuon_H
#define L1T_OmtfP1_AlgoMuon_H

#include "L1Trigger/L1TMuonOverlapPhase1/interface/AlgoMuonBase.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GoldenPatternBase.h"
#include "L1Trigger/L1TMuonOverlapPhase1/interface/Omtf/GoldenPatternResult.h"
#include <ostream>

class AlgoMuon : public AlgoMuonBase {
public:
  AlgoMuon() : AlgoMuonBase() {}

  AlgoMuon(const GoldenPatternResult& gpResult, GoldenPatternBase* gp, unsigned int refHitNumber, int bx = 0)
      : AlgoMuonBase(gp->getConfig()),
        gpResult(gpResult),
        goldenPatern(gp),
        //m_q(gpResult.getFiredLayerCnt()),  //initial value of quality, can be altered later
        m_bx(bx),
        m_rhitNumb(refHitNumber) {}

  GoldenPatternBase* getGoldenPatern() const { return goldenPatern; }

  ~AlgoMuon() override {}

  const GoldenPatternResult& getGpResult() const { return gpResult; }

  PdfValueType getDisc() const;
  int getPhi() const;
  int getEtaHw() const override { return gpResult.getEta(); }
  int getRefLayer() const { return gpResult.getRefLayer(); }
  unsigned int getFiredLayerBits() const;
  int getQ() const;// { return m_q; }
  int getBx() const { return m_bx; }

  int getPt() const {
    if (goldenPatern == nullptr)
      return -1; //TODO maybe it should be 0?
    return goldenPatern->key().thePt;
  }

  int getPtUnconstrained() const {
    if (goldenPaternUpt == nullptr)
      return -1; //TODO maybe it should be 0?
    return goldenPaternUpt->key().thePt;
  }

  //TODO how about charge Unconstrained?
  int getCharge() const {
    if (goldenPatern == nullptr)
      return 0;
    return goldenPatern->key().theCharge;
  }
  int getPhiRHit() const { return gpResult.getRefHitPhi(); }

  //TODO how about goldenPaternUpt?
  unsigned int getPatternNumber() const {
    if (goldenPatern == nullptr)
      return 0;
    return goldenPatern->key().theNumber;
  }

  //TODO how about goldenPaternUpt?
  unsigned int getHwPatternNumber() const {
    if (goldenPatern == nullptr)
      return 0;
    return goldenPatern->key().getHwPatternNumber();
  }

  unsigned int getRefHitNumber() const { return m_rhitNumb; }

  //void setQ(int q) { m_q = q; }
  void setEta(int eta) { gpResult.setEta(eta); }

  void setRefHitNumber(unsigned int aRefHitNum) { m_rhitNumb = aRefHitNum; }

  bool isValid() const override;

  unsigned int getFiredLayerCnt() const override { return gpResult.getFiredLayerCnt(); }

  double getPdfSum() const override { return gpResult.getPdfSum(); }

  const StubResult& getStubResult(unsigned int iLayer) const override { return gpResult.getStubResults().at(iLayer); }

  const StubResults& getStubResults() const override { return gpResult.getStubResults(); }

  const bool isKilled() const { return killed; }

  void kill() { killed = true; }

  bool operator<(const AlgoMuon& o) const;

  friend std::ostream& operator<<(std::ostream& out, const AlgoMuon& o);

  std::vector<std::shared_ptr<AlgoMuon> >& getKilledMuons() {
    return killedMuons;
  }
  
  const GoldenPatternResult& getGpResultUpt() const {
    return gpResultUpt;
  }

  void setGpResultUpt(const GoldenPatternResult& gpResultUpt) {
    this->gpResultUpt = gpResultUpt;
  }

  GoldenPatternBase* getGoldenPaternUpt() const {
    return goldenPaternUpt;
  }

  void setGoldenPaternUpt(GoldenPatternBase* goldenPaternUpt) {
    this->goldenPaternUpt = goldenPaternUpt;
  }

  int getChargeNN() const {
    return chargeNN;
  }

  void setChargeNN(int chargeNn = 0) {
    chargeNN = chargeNn;
  }

  int getPtNN() const {
    return ptNN;
  }

  void setPtNN(int ptNn = 0) {
    ptNN = ptNn;
  }

private:
  ///FIXME maybe the gpResult cannot be a reference or pointer, ad not a copy
  GoldenPatternResult gpResult;

  //GoldenPatternResult without vertex constraint (unconstrained pt)
  //TODO make it pointer
  GoldenPatternResult gpResultUpt;

  GoldenPatternBase* goldenPatern = nullptr;

  //GoldenPattern without vertex constraint (unconstrained pt)
  GoldenPatternBase* goldenPaternUpt = nullptr;

  //int m_q = -1;
  int m_bx = 0;

  unsigned int m_rhitNumb = 0;

  bool killed = false;

  unsigned int index = 0;

  std::vector<std::shared_ptr<AlgoMuon>> killedMuons;

  int ptNN = 0;
  int chargeNN = 0;
};

typedef std::shared_ptr<AlgoMuon> AlgoMuonPtr;
typedef std::vector<AlgoMuonPtr> AlgoMuons;

#endif  //L1T_OmtfP1_AlgoMuon_H
