
#Ahmed
import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Geometry.CMSCommonData.cmsExtendedGeometryPostLS1XML_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.MuonCommonData.muonIdealGeometryXML_upscope_cfi")
process.load("Geometry.RPCGeometry.rpcGeometry_cfi")
process.load("Geometry.CSCGeometry.cscGeometry_cfi")
process.load("Alignment.CommonAlignmentProducer.FakeAlignmentSource_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
            fileNames = cms.untracked.vstring(
            
	        # 'file:/afs/cern.ch/work/a/aliah/private/work/extention/CMSSW_7_3_2_patch3/src/SingleMuPt100_cfi_RECO_gas_leak.root',
	        'file:/afs/cern.ch/work/a/aliah/private/work/extention/CMSSW_7_3_2_patch3/src/SingleMuPt100_cfi_RECO.root',

            )
)

process.demo = cms.EDAnalyzer('genAnalyzer',
                              
							 # genParticles  = cms.InputTag("genParticles"),
							  GTReadoutRcd     = cms.InputTag("gtDigis"),
                              GMTReadoutRcd    = cms.InputTag("gtDigis"),
                              RootFileName = cms.untracked.string("output.root"),
                              #Debug        = cms.untracked.bool(False),
							


)


process.p = cms.Path(process.demo)
GenParticles = cms.string("prunedGenParticles"),




