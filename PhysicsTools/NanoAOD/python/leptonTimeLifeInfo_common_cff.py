#
# Common definition of time-life variables for pat-leptons produced
# with {Electron,Muon,Tau}TimeLifeInfoTableProducer
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
    ipLength = Var("ipLength().value()", float, doc="lenght of impact parameter (3d)", precision=10),
    ipLengthSig = Var("ipLength().significance()", float, doc="significance of impact parameter", precision=10),
    IPx = Var("ipVector().x()", float, doc="x coordinate of impact parameter vector", precision=10),
    IPy = Var("ipVector().y()", float, doc="y coordinate of impact parameter vector", precision=10),
    IPz = Var("ipVector().z()", float, doc="z coordinate of impact parameter vector", precision=10)
)

# track parameters and covariance at ref. point
trackVars = cms.PSet(
    track_qoverp = Var("?hasTrack()?track().parameter(0):0", float, doc="track q/p", precision=10),
    track_lambda = Var("?hasTrack()?track().parameter(1):0", float, doc="track lambda", precision=10),
    track_phi = Var("?hasTrack()?track().parameter(2):0", float, doc="track phi", precision=10),
    #track_deltaPhi = Var("?hasTrack()?deltaPhi(track().parameter(2), phi):0", float, doc="track phi minus lepton phi", precision=10),
    track_dxy = Var("?hasTrack()?track().parameter(3):0", float, doc="track dxy", precision=10),
    track_dsz = Var("?hasTrack()?track().parameter(4):0", float, doc="track dsz", precision=10),
    bField_z = Var("?hasTrack()?bField_z:0", float, doc="z coordinate of magnetic field at track ref. point", precision=10),
)
# track covariance elements (adding to trackVars)
for i in range(0,5):
    for j in range(i,5):
        jistr = str(j)+str(i)
        setattr(trackVars, 'track_cov'+jistr, Var("?hasTrack()?track().covariance("+str(j)+","+str(i)+"):0", float, doc="track covariance element ("+str(j)+","+str(i)+")", precision=10))

# secondary vertex
svVars = cms.PSet(
    # SV
    hasRefitSV = Var("hasSV()", bool, doc="has SV refit using miniAOD quantities"),
    refitSVx = Var("?hasSV()?sv().x():0", float, doc="x coordinate of SV", precision=10),
    refitSVy = Var("?hasSV()?sv().y():0", float, doc="y coordinate of SV", precision=10),
    refitSVz = Var("?hasSV()?sv().z():0", float, doc="z coordinate of SV", precision=10),
    refitSVchi2 = Var("?hasSV()?sv().chi2():0", float, doc="chi2 of SV fit", precision=8),
    refitSVndof = Var("?hasSV()?sv().ndof():0", float, doc="ndof of SV fit", precision=8),
    # flight-length
    #refitFlightLength = Var("?hasSV()?flightLength().value():0", float, doc="flight-length,i.e. the PV to SV distance", precision=10),
    #refitFlightLengthSig = Var("?hasSV()?flightLength().significance():0", float, doc="Significance of flight-length", precision=10)
)
# secondary vertex covariance elements (adding to svVars)
for i in range(0,3):
    for j in range(i,3):
        jistr = str(j)+str(i)
        setattr(svVars, 'refitSVcov'+jistr, Var("?hasSV()?sv().covariance("+str(j)+","+str(i)+"):0", float, doc="Covariance of SV ("+str(j)+","+str(i)+")", precision=10))

# Module to refit PV with beam-spot constraint that is not present in Run-2 samples
refittedPV = miniAODRefitVertexProducer.clone(
    srcVertices = "offlineSlimmedPrimaryVertices",
)
run2_nanoAOD_ANY.toModify(
    prod_common, pvSource = "refittedPV")
