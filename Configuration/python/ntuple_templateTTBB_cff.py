import FWCore.ParameterSet.Config as cms

from KrAFT.GeneratorTools.pileupWeight_cff import *
from KrAFT.GeneratorTools.pdfWeight_cff import *
from KrAFT.RecoSelectorTools.leptonSelector_cfi import *
from KrAFT.RecoSelectorTools.jetSelector_cfi import *
from KrAFT.RecoSelectorTools.jpsiToMuMu_cfi import *

TFileService = cms.Service("TFileService",
    fileName = cms.string("ntuple.root"),
)

from KrAFT.GenericNtuple.genericNtupleMakerTTBB_cfi import *
event.eventCounters = ["nEventsTotal", "nEventsClean", "nEventsPAT",]

MuMu = event.clone()
ElEl = event.clone()
MuEl = event.clone()
MuMu.muon.minNumber = 2
ElEl.electron.minNumber = 2
MuEl.muon.minNumber = 1
MuEl.electron.minNumber = 1

MuJets = event.clone()
MuJets.muon.minNumber = 1
MuJets.jetMET.minNumber = 3

ElJets = event.clone()
ElJets.electron.minNumber = 1
ElJets.jetMET.minNumber = 3

nEventsNtupleElEl = cms.EDProducer("EventCountProducer")
nEventsNtupleMuMu = cms.EDProducer("EventCountProducer")
nEventsNtupleMuEl = cms.EDProducer("EventCountProducer")
nEventsNtupleMuJets = cms.EDProducer("EventCountProducer")
nEventsNtupleElJets = cms.EDProducer("EventCountProducer")

ntupleSequenceElEl = cms.Sequence(
    pileupWeight + pdfWeight
  + nEventsNtupleElEl
  + jetUnc
  + goodMuons + goodElectrons * goodJets
  + jpsiToMuMu
  * ElEl
)

ntupleSequenceMuMu = cms.Sequence(
    pileupWeight + pdfWeight
  + nEventsNtupleMuMu
  + jetUnc
  + goodMuons + goodElectrons * goodJets
  + jpsiToMuMu
  * MuMu
)

ntupleSequenceMuEl = cms.Sequence(
    pileupWeight + pdfWeight
  + nEventsNtupleMuEl
  + jetUnc
  + goodMuons + goodElectrons * goodJets
  + jpsiToMuMu
  * MuEl
)

ntupleSequenceMuJets = cms.Sequence(
    pileupWeight + pdfWeight
  + nEventsNtupleMuJets
  + jetUnc
  + goodMuons + goodElectrons * goodJets
  + jpsiToMuMu
  * MuJets
)

ntupleSequenceElJets = cms.Sequence(
    pileupWeight + pdfWeight
  + nEventsNtupleElJets
  + jetUnc
  + goodMuons + goodElectrons * goodJets
  + jpsiToMuMu
  * ElJets
)

