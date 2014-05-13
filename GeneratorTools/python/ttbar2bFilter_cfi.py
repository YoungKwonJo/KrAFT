import FWCore.ParameterSet.Config as cms

ttbar2bFilter = cms.EDFilter('TTbar2bGenFilter',
        genParticlesLabel = cms.InputTag('genParticles'),
        #genParticlesLabel = cms.InputTag('genParticlesForJetsNoNu'),
        #genJetsLabel = cms.InputTag('ak5GenJets'),
        genJetsLabel = cms.InputTag('ak5GenJetsNoNu'),
)

