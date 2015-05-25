import FWCore.ParameterSet.Config as cms

process = cms.Process("Prod")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("Configuration.Geometry.GeometryExtended2015Reco_cff")
#process.load("Geometry.CMSCommonData.cmsExtendedGeometry2023XML_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
#process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")

process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'/store/data/Commissioning2015/Cosmics/RECO/PromptReco-v1/000/234/304/00000/0060FA00-11B7-E411-A028-02163E0123B6.root'
        #'/store/data/Commissioning2015/MinimumBias/RECO/PromptReco-v1/000/245/194/00000/00E5F0C3-4101-E511-9348-02163E01448D.root'
        'file:00E5F0C3-4101-E511-9348-02163E01448D.root'
        #'file:/afs/cern.ch/work/a/aocampor/public/CMSSW_7_3_5/src/CMSRPCDPGUserCode/RPCHitsFilter/test/FC76C5CE-9EFF-E411-8594-02163E01383E.root'
    )
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('MyRecHitsCollection.root')
)

process.demo = cms.EDProducer('MyRPCProducer',
                      src = cms.InputTag('rpcRecHits')
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
