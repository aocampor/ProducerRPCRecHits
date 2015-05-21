import FWCore.ParameterSet.Config as cms

demo = cms.EDProducer('MyRPCProducer',
                      src = cms.InputTag('AhmedRecHitCollection')
)
