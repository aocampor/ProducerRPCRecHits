
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
            #'file:/afs/cern.ch/user/a/aliah/source/001FB128-1DD9-E111-BBB3-E41F1318165C.root'
			#'file:/afs/cern.ch/work/a/aliah/private/work/emulator/CMSSW_7_3_2_patch3/src/L1Trigger/RPCTrigger/test/l1_gas_leak.root',
			'file:/afs/cern.ch/work/a/aliah/private/work/producer/CMSSW_7_3_2_patch3/src/SingleMuPt100_cfi_RECO.root',

            )
)

process.demo = cms.EDAnalyzer('genAnalyzer',
                              RootFileName = cms.untracked.string("genmuons.root"),


)


process.p = cms.Path(process.demo)
GenParticles = cms.string("prunedGenParticles"),


