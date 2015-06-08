import FWCore.ParameterSet.Config as cms

process = cms.Process("Prod")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.Geometry.GeometryExtended2015Reco_cff")
#process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023XML_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
#process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(

   'root://eoscms//eos/cms/store/data/Run2015A/RPCMonitor/RAW/v1/000/247/081/00000/1A498EDF-460B-E511-8828-02163E0121C5.root'
    )
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('MyRPCDigiCollection.root')
)

process.demo = cms.EDProducer('MyRPCProducer',
                      src = cms.InputTag('hltMuonRPCDigis')
)


# process.MessageLogger = cms.Service("MessageLogger",
#     debugModules = cms.untracked.vstring('*'),
#     cout = cms.untracked.PSet(threshold = cms.untracked.string('INFO')),
#     # cout = cms.untracked.PSet(threshold = cms.untracked.string('DEBUG')),
#     destinations = cms.untracked.vstring('cout')
# ) 

process.load("RecoLocalMuon.RPCRecHit.rpcRecHits_cfi")
#process.p = cms.Path(process.rpcRecHits)
process.p = cms.Path(process.demo)
process.ep = cms.EndPath(process.out)