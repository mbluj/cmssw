#
# Common definition of time-life variables for pat-leptons updated 
# with pat{Electron,Muon,Tau}TimeLifeInfoUpdater
#
import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.nano_eras_cff import *
from PhysicsTools.PatAlgos.miniAODRefitVertexProducer_cfi import miniAODRefitVertexProducer

# common settings of lepton life-time info producer
prod_common = cms.PSet(
    pvSource = cms.InputTag("offlineSlimmedPrimaryVerticesWithBS"),
    pvChoice = cms.int32(0) #0: PV[0], 1: smallest dz
)

# impact parameter
ipVars = cms.PSet(
    ipLength = Var("?hasUserFloat('IP')?userFloat('IP'):0", float, doc="lenght of impact parameter (3d)", precision=16),
    ipLengthSig = Var("?hasUserFloat('IP_sig')?userFloat('IP_sig'):0", float, doc="Significance of impact parameter", precision=16),
    IPx = Var("?hasUserFloat('IP_x')?userFloat('IP_x'):0", float, doc="x coordinate of impact parameter vector", precision=16),
    IPy = Var("?hasUserFloat('IP_x')?userFloat('IP_x'):0", float, doc="x coordinate of impact parameter vector", precision=16),
    IPz = Var("?hasUserFloat('IP_x')?userFloat('IP_x'):0", float, doc="x coordinate of impact parameter vector", precision=16)
)

# track parameters and covariance at ref. point
trackVars = cms.PSet(
    track_qoverp = Var("?hasUserFloat('track_qoverp')?userFloat('track_qoverp'):0", float, doc="track q/p", precision=16),
    track_lambda = Var("?hasUserFloat('track_lambda')?userFloat('track_lambda'):0", float, doc="track lambda", precision=16),
    track_phi = Var("?hasUserFloat('track_phi')?userFloat('track_phi'):0", float, doc="track phi", precision=16),
    track_dxy = Var("?hasUserFloat('track_dxy')?userFloat('track_dxy'):0", float, doc="track dxy", precision=16),
    track_dsz = Var("?hasUserFloat('track_dsz')?userFloat('track_dsz'):0", float, doc="track dsz", precision=16),
    bField_z = Var("?hasUserFloat('bField_z')?userFloat('bField_z'):0", float, doc="z coordinate of magnetic field at track ref. point", precision=16),
)
# track covariance elements (adding to trackVars)
for i in range(0,5):
    for j in range(i,5):
        jistr = str(j)+str(i)
        setattr(trackVars, 'track_cov'+jistr,Var("?hasUserFloat('trackCov_"+jistr+"')?userFloat('trackCov_"+jistr+"'):0", float, doc="track covariance element ("+str(j)+","+str(i)+")", precision=16))

# secondary vertex
svVars = cms.PSet(
    # SV
    SVstatus = Var("?hasUserInt('refitSV_status')?userInt('refitSV_status'):-1", int, doc="Status of SV refit: 0: not fitted or not present, 1: OK, -1: failed"),
    SVx = Var("?hasUserFloat('refitSV_x')?userFloat('refitSV_x'):0", float, doc="x coordinate of SV", precision=16),
    SVy = Var("?hasUserFloat('refitSV_y')?userFloat('refitSV_y'):0", float, doc="y coordinate of SV", precision=16),
    SVz = Var("?hasUserFloat('refitSV_z')?userFloat('refitSV_z'):0", float, doc="z coordinate of SV", precision=16),
    SVchi2 = Var("?hasUserFloat('refitSV_chi2')?userFloat('refitSV_chi2'):0", float, doc="chi2 of SV fit", precision=10),
    SVndof = Var("?hasUserFloat('refitSV_ndof')?userFloat('refitSV_ndof'):0", float, doc="ndof of SV fit", precision=10),
    SVcxx = Var("?hasUserFloat('refitSV_cxx')?userFloat('refitSV_cxx'):0", float, doc="Covariance of SV (0,0)", precision=16),
    SVcyx = Var("?hasUserFloat('refitSV_cyx')?userFloat('refitSV_cyx'):0", float, doc="Covariance of SV (1,0)", precision=16),
    SVczx = Var("?hasUserFloat('refitSV_czx')?userFloat('refitSV_czx'):0", float, doc="Covariance of SV (2,0)", precision=16),
    SVcyy = Var("?hasUserFloat('refitSV_cyy')?userFloat('refitSV_cyy'):0", float, doc="Covariance of SV (1,1)", precision=16),
    SVczy = Var("?hasUserFloat('refitSV_czy')?userFloat('refitSV_czy'):0", float, doc="Covariance of SV (2,1)", precision=16),
    SVczz = Var("?hasUserFloat('refitSV_czz')?userFloat('refitSV_czz'):0", float, doc="Covariance of SV (2,2)", precision=16),
    # flight-length
    flightLength = Var("?hasUserFloat('flightLength')?userFloat('flightLength'):0", float, doc="flight-length,i.e. the PV to SV distance", precision=16),
    flightLengthSig = Var("?hasUserFloat('flightLength_sig')?userFloat('flightLength_sig'):0", float, doc="Significance of flight-length", precision=16)
)

# Module to refit PV with beam-spot constraint that is not present in Run-2 samples
refittedPV = miniAODRefitVertexProducer.clone(
    srcVertices = "offlineSlimmedPrimaryVertices",
)
run2_nanoAOD_ANY.toModify(
    prod_common, pvSource = "refittedPV")
