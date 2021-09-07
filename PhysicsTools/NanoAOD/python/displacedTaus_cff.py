import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

##################### Displaced taus collection #################################
# Production of displaced taus to be added to the master nano_cff by the
# following dedicated function
def nanoAOD_addDisplacedTaus(process):
    print("Add displaced taus")
    postfix = 'Displaced'
    from RecoTauTag.Configuration.tools.adaptToRunAtMiniAOD import adaptToRunAtMiniAOD
    tauAtMiniTools = adaptToRunAtMiniAOD(process,postfix=postfix)
    tauAtMiniTools.addTauReReco()
    tauAtMiniTools.adaptTauToMiniAODReReco()
    #Remove PATTau MC info as it is anyway recalculated for NanoAOD
    from PhysicsTools.PatAlgos.tools.coreTools import runOnData
    runOnData(process, names = ['Taus'], outputModules = [],postfix=postfix)
    #modify tau reco to be displaced-friendly
    pvAlgo = 'highestPtInEvent' #PV[0]
    maxDeltaZ = 100.
    maxTIP = 100.
    maxTrkChi2 = 1000.
    minPxlForDM = 1 #default 1
    qCuts = [
        getattr(process,'combinatoricRecoTaus'+postfix).builders[0].qualityCuts,
        getattr(process,'ak4PFJetsRecoTauChargedHadrons'+postfix).builders[0].qualityCuts,
        getattr(process,'ak4PFJetsRecoTauChargedHadrons'+postfix).builders[1].qualityCuts,
        getattr(process,'ak4PFJetsRecoTauChargedHadrons'+postfix).builders[2].qualityCuts,
        getattr(process,'ak4PFJetsLegacyHPSPiZeros'+postfix).builders[0].qualityCuts, #needed?
        getattr(process,'hpsPFTauPrimaryVertexProducer'+postfix).qualityCuts #need to refit proper PV
    ]
    for qCut in qCuts:
        qCut.pvFindingAlgo = pvAlgo
        qCut.signalQualityCuts.maxDeltaZ = maxDeltaZ
        qCut.signalQualityCuts.maxTransverseImpactParameter = maxTIP
        qCut.signalQualityCuts.maxTrackChi2 = maxTrkChi2
    getattr(process,'hpsSelectionDiscriminator'+postfix).minPixelHits = minPxlForDM
    getattr(process,'hpsPFTauDiscriminationByDecayModeFindingNewDMs'+postfix).minPixelHits = minPxlForDM
    getattr(process,'hpsPFTauDiscriminationByDecayModeFindingOldDMs'+postfix).minPixelHits = minPxlForDM
    getattr(process,'hpsPFTauDiscriminationByDecayModeFinding'+postfix).minPixelHits = minPxlForDM
    #add it to sequences
    process.finalDisplacedTaus.src = 'selectedPatTaus'+postfix
    process.displacedTauTask.add(getattr(process,'miniAODTausTask'+postfix))
    process.nanoTableTaskCommon.add(process.displacedTauTask)
    process.nanoTableTaskCommon.add(process.displacedTauTablesTask)
    process.nanoTableTaskFS.add(process.displacedTauMCTask)
    #FIXME: to be removed after studies
    print("Relax standard taus")
    process.finalTaus.cut = process.finalDisplacedTaus.cut.value()

##################### Import reusable funtions and objects from std taus ######## 
from PhysicsTools.NanoAOD.taus_cff import _tauId2WPMask,_tauId4WPMask,_tauId5WPMask,_tauId6WPMask,_tauId7WPMask,_tauId8WPMask,tausMCMatchLepTauForTable,tausMCMatchHadTauForTable,tauMCTable

##################### User floats producers, selectors ##########################

finalDisplacedTaus = cms.EDFilter("PATTauRefSelector",
    src = cms.InputTag("selectedPatTaus"),
    cut = cms.string("pt > 18 && tauID('decayModeFindingNewDMs')")
)

##################### Tables for final output and docs ##########################

displacedTauTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("finalDisplacedTaus"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name= cms.string("DisplacedTau"),
    doc = cms.string("displacedTaus after basic selection (" + finalDisplacedTaus.cut.value()+")"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the taus
    variables = cms.PSet() # PSet defined below in era dependent way
)
_tauVarsBase = cms.PSet(P4Vars,
       charge = Var("charge", int, doc="electric charge"),
       jetIdx = Var("?hasUserCand('jet')?userCand('jet').key():-1", int, doc="index of the associated jet (-1 if none)"),#MB: doesn't work w/o cross-linking, can add crosslinker it usefull
       decayMode = Var("decayMode()",int),
       idDecayMode = Var("tauID('decayModeFinding')", bool),
       idDecayModeNewDMs = Var("tauID('decayModeFindingNewDMs')", bool),

       leadTkPtOverTauPt = Var("leadChargedHadrCand.pt/pt ",float, doc="pt of the leading track divided by tau pt",precision=10),
       leadTkDeltaEta = Var("leadChargedHadrCand.eta - eta ",float, doc="eta of the leading track, minus tau eta",precision=8),
       leadTkDeltaPhi = Var("deltaPhi(leadChargedHadrCand.phi, phi) ",float, doc="phi of the leading track, minus tau phi",precision=8),

       dxy = Var("leadChargedHadrCand().dxy()",float, doc="d_{xy} of lead track with respect to PV, in cm (with sign)",precision=10),
       dz = Var("leadChargedHadrCand().dz()",float, doc="d_{z} of lead track with respect to PV, in cm (with sign)",precision=14),

       # these are too many, we may have to suppress some
       rawIso = Var( "tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits')", float, doc = "combined isolation (deltaBeta corrections)", precision=10),
       rawIsodR03 = Var( "(tauID('chargedIsoPtSumdR03')+max(0.,tauID('neutralIsoPtSumdR03')-0.072*tauID('puCorrPtSum')))", float, doc = "combined isolation (deltaBeta corrections, dR=0.3)", precision=10),
       chargedIso = Var( "tauID('chargedIsoPtSum')", float, doc = "charged isolation", precision=10),
       neutralIso = Var( "tauID('neutralIsoPtSum')", float, doc = "neutral (photon) isolation", precision=10),
       puCorr = Var( "tauID('puCorrPtSum')", float, doc = "pileup correction", precision=10),
       photonsOutsideSignalCone = Var( "tauID('photonPtSumOutsideSignalCone')", float, doc = "sum of photons outside signal cone", precision=10),

       idAntiMu = _tauId2WPMask("againstMuon%sSimple", doc= "Anti-muon discriminator Simple: "),
       idAntiEleDeadECal = Var("tauID('againstElectronDeadECAL')", bool, doc = "Anti-electron dead-ECal discriminator"),

)

_mvaIsoVars = cms.PSet(#FIXME check names when switched to taus@mini
    rawMVAnewDM = Var( "tauID('byIsolationMVArun2v1DBnewDMwLTraw')",float, doc="byIsolationMVArun2v1DBoldDMwLT raw output discriminator",precision=10),
    rawMVAoldDM = Var( "tauID('byIsolationMVArun2v1DBoldDMwLTraw')",float, doc="byIsolationMVArun2v1DBoldDMwLT raw output discriminator",precision=10),
    rawMVAoldDMdR03 = Var( "tauID('byIsolationMVArun2v1DBdR03oldDMwLTraw')",float, doc="byIsolationMVArun2v1DBoldDMwLT raw output discriminator (2015)",precision=10),
    idMVAnewDM = _tauId6WPMask( "by%sIsolationMVArun2v1DBnewDMwLT", doc="IsolationMVArun2v1DBnewDMwLT ID working point"),
    idMVAoldDM = _tauId6WPMask( "by%sIsolationMVArun2v1DBoldDMwLT", doc="IsolationMVArun2v1DBoldDMwLT ID working point"),
    idMVAoldDMdR03 = _tauId6WPMask( "by%sIsolationMVArun2v1DBdR03oldDMwLT", doc="IsolationMVArun2v1DBoldDMdR0p3wLT ID working point")
)
_mvaAntiEVars = cms.PSet(#FIXME check names when switched to taus@mini
       rawAntiEle = Var("tauID('againstElectronMVA6Raw2018')", float, doc= "Anti-electron MVA discriminator V6 raw output discriminator", precision=10),
       rawAntiEleCat = Var("tauID('againstElectronMVA6category2018')", int, doc="Anti-electron MVA discriminator V6 category"),
       idAntiEle = _tauId5WPMask("againstElectron%sMVA62018", doc= "Anti-electron MVA discriminator V6"),
)
_deepTauVars2017v2p1 = cms.PSet(
    rawDeepTau2017v2p1VSe = Var("tauID('byDeepTau2017v2p1VSeraw')", float, doc="byDeepTau2017v2p1VSe raw output discriminator (deepTau2017v2p1)", precision=10),
    rawDeepTau2017v2p1VSmu = Var("tauID('byDeepTau2017v2p1VSmuraw')", float, doc="byDeepTau2017v2p1VSmu raw output discriminator (deepTau2017v2p1)", precision=10),
    rawDeepTau2017v2p1VSjet = Var("tauID('byDeepTau2017v2p1VSjetraw')", float, doc="byDeepTau2017v2p1VSjet raw output discriminator (deepTau2017v2p1)", precision=10),
    idDeepTau2017v2p1VSe = _tauId8WPMask("by%sDeepTau2017v2p1VSe", doc="byDeepTau2017v2p1VSe ID working points (deepTau2017v2p1)"),
    idDeepTau2017v2p1VSmu = _tauId4WPMask("by%sDeepTau2017v2p1VSmu", doc="byDeepTau2017v2p1VSmu ID working points (deepTau2017v2p1)"),
    idDeepTau2017v2p1VSjet = _tauId8WPMask("by%sDeepTau2017v2p1VSjet", doc="byDeepTau2017v2p1VSjet ID working points (deepTau2017v2p1)"),
)

_variablesMiniV2 = cms.PSet(
    _tauVarsBase,
    _mvaAntiEVars,
    _mvaIsoVars,
    #_deepTauVars2017v2p1 #FIXME
)

displacedTauTable.variables = _variablesMiniV2

displacedTausMCMatchLepTauForTable = tausMCMatchLepTauForTable.clone(
    src = displacedTauTable.src
)

displacedTausMCMatchHadTauForTable = tausMCMatchHadTauForTable.clone(
    src         = displacedTauTable.src
)

displacedTauMCTable = tauMCTable.clone(
    src = displacedTauTable.src,
    mcMap = cms.InputTag("displacedTausMCMatchLepTauForTable"),
    mcMapVisTau = cms.InputTag("displacedTausMCMatchHadTauForTable"),
    objName = displacedTauTable.name
)


displacedTauTask = cms.Task(finalDisplacedTaus)
displacedTauTablesTask = cms.Task(displacedTauTable)
displacedTauMCTask = cms.Task(displacedTausMCMatchLepTauForTable,displacedTausMCMatchHadTauForTable,displacedTauMCTable)
