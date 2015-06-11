import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring( *(

	
	'file:/afs/cern.ch/work/a/aliah/private/emulator/CMSSW_7_4_2/src/L1Trigger/RPCTrigger/test/l1_leaky_filter.root',






    ) )
)
	

process.demo = cms.EDAnalyzer('MyTriggerAnalyzer',
                              GTReadoutRcd     = cms.InputTag("hltGtDigis"),
                              GMTReadoutRcd    = cms.InputTag("hltGtDigis"),
                              RootFileName = cms.untracked.string("work.root"),
                              Debug        = cms.untracked.bool(False),

)


process.p = cms.Path(process.demo)
