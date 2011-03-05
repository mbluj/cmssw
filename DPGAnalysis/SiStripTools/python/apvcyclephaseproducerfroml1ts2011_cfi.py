import FWCore.ParameterSet.Config as cms

APVPhases = cms.EDProducer('APVCyclePhaseProducerFromL1TS',
                           defaultPartitionNames = cms.vstring("TI",
                                                               "TO",
                                                               "TP",
                                                               "TM"
                                                               ),
                           defaultPhases = cms.vint32(60,60,60,60),
                           l1TSCollection = cms.InputTag("scalersRawToDigi"),
                           )

