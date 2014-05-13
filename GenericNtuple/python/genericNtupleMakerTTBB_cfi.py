import FWCore.ParameterSet.Config as cms

event = cms.EDAnalyzer("KGenericNtupleMakerTTBB",
    isMC = cms.bool(False),

    genEventInfo = cms.InputTag("generator"),
    genParticle = cms.InputTag("genParticles"),
#    genParticle = cms.InputTag("genParticlesForJetsNoNu"),
    genJet = cms.InputTag('ak5GenJetsNoNu'),
#    genJet = cms.InputTag("ak5GenJets"),
    recoToGenJetMap = cms.InputTag("recoToGenJetMap"),
    genJetToPartonsMap = cms.InputTag("genJetToPartonsMap"),

    flavorHistoryFilter = cms.InputTag("flavorHistoryFilter"),
    bFlavorHistoryProducer = cms.InputTag("bFlavorHistoryProducer", "bPartonFlavorHistory"),
    cFlavorHistoryProducer = cms.InputTag("cFlavorHistoryProducer", "cPartonFlavorHistory"),

    pdfWeights = cms.InputTag("pdfWeight"),
    puWeight = cms.InputTag("pileupWeight"),
    vertex = cms.InputTag("goodOfflinePrimaryVertices"),
    eventCounters = cms.vstring(),

    electron = cms.PSet(
        src = cms.InputTag("goodElectrons"),
        minNumber = cms.uint32(0),
        maxNumber = cms.uint32(999),
    ),
    muon = cms.PSet(
        src = cms.InputTag("goodMuons"),
        minNumber = cms.uint32(0),
        maxNumber = cms.uint32(999),
    ),
    jetMET = cms.PSet(
        jet = cms.InputTag("goodJets"),
        met = cms.InputTag("patMETsPFlow"),
        unc = cms.InputTag("jetUnc"),
        minNumber = cms.uint32(0),
        leptonDeltaR = cms.double(0.5),
        bTagType = cms.string("combinedSecondaryVertexBJetTags"),
    ),
    jpsi = cms.PSet(
        src = cms.InputTag("jpsiToMuMu"),
    ),
)
