import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Services_cff")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
    '/store/data/Run2012D/MuEG/AOD/22Jan2013-v1/20000/FE5E0B33-0D82-E211-AB78-00237DA15C7C.root',
]

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("skim.root"),
    outputCommands = cms.untracked.vstring(
        "keep *",
        "drop *_selectedMuons_*_SKIM",
        "drop *_selectedElectrons_*_SKIM",
        "drop *_z*Cands_*_SKIM",
        "drop *_selectedJets_*_SKIM",
    ),
    SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
)
process.outPath = cms.EndPath(process.out)

process.load("TopAnalysis.Configuration.commonFilters_cff")
process.load("HLTrigger.HLTfilters.hltHighLevel_cfi")
process.hltHighLevel.throw = False
process.hltHighLevel.HLTPaths = ["HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*", "HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",]

process.p = cms.Path(
    process.hltHighLevel
  + process.commonSequenceForData
  + process.filterMuEGSequence
)

