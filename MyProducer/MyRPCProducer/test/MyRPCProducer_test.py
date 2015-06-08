import FWCore.ParameterSet.Config as cms

process = cms.Process("NewRPCDigis")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")

process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.MuonCommonData.muonIdealGeometryXML_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:myfile.root'
        #'/store/data/Run2015A/RPCMonitor/RAW/v1/000/247/444/00000/A291CEFA-EA0D-E511-96CF-02163E011E08.root',
        '/store/data/Run2015A/RPCMonitor/RAW/v1/000/247/398/00000/105F4AFC-870D-E511-854C-02163E011B28.root',
        '/store/data/Run2015A/RPCMonitor/RAW/v1/000/247/398/00000/2CBB45FC-870D-E511-831B-02163E011D69.root',
        '/store/data/Run2015A/RPCMonitor/RAW/v1/000/247/398/00000/387779C8-880D-E511-9FAF-02163E012298.root',
        '/store/data/Run2015A/RPCMonitor/RAW/v1/000/247/398/00000/468D7155-940D-E511-926A-02163E01456E.root',
        '/store/data/Run2015A/RPCMonitor/RAW/v1/000/247/398/00000/508F5E9C-900D-E511-9073-02163E01257B.root',
        '/store/data/Run2015A/RPCMonitor/RAW/v1/000/247/398/00000/603FEDF6-870D-E511-90D0-02163E0146EE.root',
        '/store/data/Run2015A/RPCMonitor/RAW/v1/000/247/398/00000/8AA6B78D-900D-E511-9880-02163E014565.root',
        
    )
)

process.myProducerLabel = cms.EDProducer('MyRPCProducer',
                                         src = cms.InputTag('hltMuonRPCDigis')
)

process.p = cms.Path( process.myProducerLabel )

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('MyDigiCollection.root')
    
)

  
process.p = cms.Path(process.myProducerLabel)

process.e = cms.EndPath(process.out)
