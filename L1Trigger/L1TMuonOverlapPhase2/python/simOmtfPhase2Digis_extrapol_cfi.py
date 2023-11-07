import FWCore.ParameterSet.Config as cms

###OMTF emulator configuration
simOmtfPhase2Digis = cms.EDProducer("L1TMuonOverlapPhase2TrackProducer",
                              
  srcDTPh = cms.InputTag('simDtTriggerPrimitiveDigis'),
  srcDTTh = cms.InputTag('simDtTriggerPrimitiveDigis'),
  srcCSC = cms.InputTag('simCscTriggerPrimitiveDigis','MPCSORTED'),
  srcRPC = cms.InputTag('simMuonRPCDigis'), 
  srcDTPhPhase2 = cms.InputTag('dtTriggerPhase2PrimitiveDigis'),
  
  #g4SimTrackSrc = cms.InputTag('g4SimHits'),                             
  dumpResultToXML = cms.bool(False),
  dumpDetailedResultToXML = cms.bool(False),
  XMLDumpFileName = cms.string("TestEvents.xml"),                                     
  dumpGPToXML = cms.bool(False),  
  readEventsFromXML = cms.bool(False),
  eventsXMLFiles = cms.vstring("TestEvents.xml"),
  
  dropRPCPrimitives = cms.bool(False),                                    
  dropDTPrimitives = cms.bool(False),                                    
  dropCSCPrimitives = cms.bool(False),
  processorType = cms.string("OMTFProcessor"),
  
  #ghostBusterType = cms.string("GhostBusterPreferRefDt"),
  
  #patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/Patterns_0x00020007.xml")
  #patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/Patterns_0x0003.xml")                               

  #if commented the default values are 0-0
  #-3 to 4 is the range of the OMTF DAQ readout, so should be used e.g. in the DQM data to emulator comparison
  bxMin = cms.int32(0),
  bxMax = cms.int32(0),
  
  minDtPhiQuality = cms.int32(2),
  minDtPhiBQuality = cms.int32(4),
  
  dtRefHitMinQuality =  cms.int32(4),  
  
  usePhiBExtrapolationFromMB1 = cms.bool(True),
  usePhiBExtrapolationFromMB2 = cms.bool(True),
  useStubQualInExtr  = cms.bool(False),
  useEndcapStubsRInExtr  = cms.bool(False),
  useFloatingPointExtrapolation  = cms.bool(False),

  extrapolFactorsFilename = cms.string("ExtrapolationFactors_simple.xml"),
  sorterType = cms.string("byLLH"),
  ghostBusterType = cms.string("byRefLayer") # byLLH byRefLayer GhostBusterPreferRefDt
)

