import FWCore.ParameterSet.Config as cms

def initialise(runOnMC, decayMode, doOutModule=False, doPAT=True):
    process = cms.Process("KrAFT")

    process.load("Configuration.StandardSequences.Services_cff")
    process.load("Configuration.Geometry.GeometryDB_cff")
    process.load("Configuration.StandardSequences.MagneticField_cff")
    process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

    process.load("FWCore.MessageLogger.MessageLogger_cfi")
    process.MessageLogger.cerr.FwkReport.reportEvery = 10000

    process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
    process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*")

    process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
    from Configuration.AlCa.autoCond import autoCond
    if runOnMC: process.GlobalTag.globaltag = autoCond['startup']
    else: process.GlobalTag.globaltag = autoCond['com10']

    outputModuleForTriggerMatch = ""
    outputModules = []
    if doOutModule:
        from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
        process.out = cms.OutputModule("PoolOutputModule",
            fileName = cms.untracked.string("out.root"),
            outputCommands = cms.untracked.vstring(
                'drop *',
                'keep recoPFCandidates_particleFlow_*_*',
                *patEventContentNoCleaning
            )
        )
        process.outPath = cms.EndPath(process.out)

        outputModuleForTriggerMatch = "out"
        outputModules.append(process.out)

    ## Load PAT
    process.load("PhysicsTools.PatAlgos.patSequences_cff")

    ## Apply MVA
    process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
    process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
    process.patElectrons.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
    process.patElectrons.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
    process.patDefaultSequence.replace( process.patElectrons, process.eidMVASequence * process.patElectrons )

    ## Load trigger matching
    process.load("KrAFT.Configuration.hltFilters_cff")
    #from PhysicsTools.PatAlgos.tools.trigTools import *
    #switchOnTriggerMatchEmbedding(process, outputModule="")

    ## Apply PF2PAT
    from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
    if runOnMC: jecLevels = ['L1FastJet','L2Relative','L3Absolute']
    else: jecLevels = ['L1FastJet','L2Relative', 'L3Absolute', 'L2L3Residual']

    postfix="PFlow"
    #usePFBRECO(process,runPFBRECO=True,
    usePF2PAT(process, runPF2PAT=True,
              runOnMC=runOnMC, outputModules = outputModules, postfix=postfix,
              jetAlgo="AK5", jetCorrections=("AK5PFchs", jecLevels),
              pvCollection=cms.InputTag('goodOfflinePrimaryVertices'),
              typeIMetCorrections=True)

    #
    process.pfPileUpPFlow.checkClosestZVertex = cms.bool(False)

    # top projections in PF2PAT:
    process.pfNoPileUpPFlow.enable = True
    process.pfNoMuonPFlow.enable = True
    process.pfNoElectronPFlow.enable = True
    process.pfNoTauPFlow.enable = False
    process.pfNoJetPFlow.enable = False

    # verbose flags for the PF2PAT modules
    process.pfNoMuonPFlow.verbose = False

    # Use non-isolated muons and electrons
    process.patMuonsPFlow.pfMuonSource = "pfMuonsPFlow"
    process.patElectronsPFlow.pfElectronSource = "pfElectronsPFlow"

    # And turn on delta-beta corrections while building pfIsolated*PFlow
    process.pfIsolatedMuonsPFlow.doDeltaBetaCorrection = True
    process.pfIsolatedElectronsPFlow.doDeltaBetaCorrection = True

    # Change DR cone size to 0.3
    process.pfIsolatedMuonsPFlow.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('muPFIsoValueCharged03PFlow'))
    process.pfIsolatedMuonsPFlow.deltaBetaIsolationValueMap = cms.InputTag('muPFIsoValuePU03PFlow')
    process.pfIsolatedMuonsPFlow.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('muPFIsoValueNeutral03PFlow'),
                                                                         cms.InputTag('muPFIsoValueGamma03PFlow'),)
    process.pfMuonsPFlow.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('muPFIsoValueCharged03PFlow') )
    process.pfMuonsPFlow.deltaBetaIsolationValueMap = cms.InputTag('muPFIsoValuePU03PFlow')
    process.pfMuonsPFlow.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('muPFIsoValueNeutral03PFlow'),
                                                                 cms.InputTag('muPFIsoValueGamma03PFlow'),)
    process.patMuonsPFlow.isolationValues.pfNeutralHadrons   = cms.InputTag('muPFIsoValueNeutral03PFlow')
    process.patMuonsPFlow.isolationValues.pfChargedAll       = cms.InputTag('muPFIsoValueChargedAll03PFlow')
    process.patMuonsPFlow.isolationValues.pfPUChargedHadrons = cms.InputTag('muPFIsoValuePU03PFlow')
    process.patMuonsPFlow.isolationValues.pfPhotons          = cms.InputTag('muPFIsoValueGamma03PFlow')
    process.patMuonsPFlow.isolationValues.pfChargedHadrons   = cms.InputTag('muPFIsoValueCharged03PFlow')

    process.pfIsolatedElectronsPFlow.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFlow'))
    process.pfIsolatedElectronsPFlow.deltaBetaIsolationValueMap = cms.InputTag('elPFIsoValuePU03PFIdPFlow')
    process.pfIsolatedElectronsPFlow.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('elPFIsoValueNeutral03PFIdPFlow'),
                                                                             cms.InputTag('elPFIsoValueGamma03PFIdPFlow'))
    process.pfElectronsPFlow.isolationValueMapsCharged  = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFlow'))
    process.pfElectronsPFlow.deltaBetaIsolationValueMap = cms.InputTag('elPFIsoValuePU03PFIdPFlow')
    process.pfElectronsPFlow.isolationValueMapsNeutral  = cms.VInputTag(cms.InputTag('elPFIsoValueNeutral03PFIdPFlow'),
                                                                     cms.InputTag('elPFIsoValueGamma03PFIdPFlow'))
    process.patElectronsPFlow.isolationValues.pfNeutralHadrons   = cms.InputTag('elPFIsoValueNeutral03PFIdPFlow')
    process.patElectronsPFlow.isolationValues.pfChargedAll       = cms.InputTag('elPFIsoValueChargedAll03PFIdPFlow')
    process.patElectronsPFlow.isolationValues.pfPUChargedHadrons = cms.InputTag('elPFIsoValuePU03PFIdPFlow')
    process.patElectronsPFlow.isolationValues.pfPhotons          = cms.InputTag('elPFIsoValueGamma03PFIdPFlow')
    process.patElectronsPFlow.isolationValues.pfChargedHadrons   = cms.InputTag('elPFIsoValueCharged03PFIdPFlow')

    ## Add common filters
    process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
        vertexCollection = cms.InputTag('offlinePrimaryVertices'),
        minimumNDOF = cms.uint32(4) ,
        maxAbsZ = cms.double(24), 
        maxd0 = cms.double(2) 
    )

    from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
    process.goodOfflinePrimaryVertices = cms.EDFilter(
        "PrimaryVertexObjectFilter",
        filterParams = pvSelector.clone( maxZ = cms.double(24.0) ),
        src=cms.InputTag('offlinePrimaryVertices')
    )

    #process.load( "TopQuarkAnalysis.Configuration.patRefSel_goodVertex_cfi" )
    #process.goodOfflinePrimaryVertices.filter = True

    process.load( 'TopQuarkAnalysis.Configuration.patRefSel_eventCleaning_cff' )
    process.trackingFailureFilter.VertexSource = cms.InputTag('goodOfflinePrimaryVertices')
    if runOnMC: process.eventCleaning += process.eventCleaningMC
    else: process.eventCleaning += process.eventCleaningData

    #reason : removing many events.
    process.eventCleaning.remove(process.CSCTightHaloFilter)

    ## Lepton veto filters for L+J channels
    process.muonVetoFilter = cms.EDFilter("PATCandViewCountFilter",
        src = cms.InputTag("selectedPatMuonsPFlow"),
        maxNumber = cms.uint32(0),
        minNumber = cms.uint32(0),
    )
    process.electronVetoFilter = cms.EDFilter("PATCandViewCountFilter",
        src = cms.InputTag("selectedPatElectronsPFlow"),
        maxNumber = cms.uint32(0),
        minNumber = cms.uint32(0),
    )

    # event counters
    process.nEventsTotal = cms.EDProducer("EventCountProducer")
    process.nEventsClean = cms.EDProducer("EventCountProducer")
    process.nEventsPAT   = cms.EDProducer("EventCountProducer")
    process.nEventsHLTElEl = cms.EDProducer("EventCountProducer")
    process.nEventsHLTMuMu = cms.EDProducer("EventCountProducer")
    process.nEventsHLTMuEl = cms.EDProducer("EventCountProducer")
    process.nEventsHLTMuJets = cms.EDProducer("EventCountProducer")
    process.nEventsHLTElJets = cms.EDProducer("EventCountProducer")


    ###############################
    ### Add AK5GenJetsNoNu ####
    ###############################
    #process.load("RecoJets.Configuration.GenJetParticles_cff")
    from RecoJets.Configuration.GenJetParticles_cff import genParticlesForJetsNoNu
    from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
    process.ak5GenJetsNoNu = ak5GenJets.clone( src = cms.InputTag("genParticlesForJetsNoNu") )

    process.ak5GenJetsSeq = cms.Sequence(process.genParticlesForJets*process.genParticlesForJetsNoNu*process.ak5GenJetsNoNu)

    # Flavor history stuff
    process.load("PhysicsTools.HepMCCandAlgos.flavorHistoryPaths_cfi")
    process.flavorHistoryFilter.pathToSelect = cms.int32(-1)
    process.cFlavorHistoryProducer.matchedSrc = cms.InputTag("ak5GenJetsNoNu")
    process.bFlavorHistoryProducer.matchedSrc = cms.InputTag("ak5GenJetsNoNu")
    #ptMinParticle = cms.double(0.0)
    #ptMinShower = cms.double(0.0)

    #process.load( "RecoJets.Configuration.GenJetParticles_cff")
    #process.load( "RecoJets.Configuration.RecoGenJets_cff")

    #process.ReGenJetSequence = cms.Sequence(
    #    process.genParticlesForJets
    #  + process.genParticlesForJetsNoNu
    #  + process.ak5GenJets
    #  + process.ak5GenJetsNoNu
    #)

    process.commonFilterSequence = cms.Sequence(
        process.goodOfflinePrimaryVertices
      * process.eventCleaning
      * process.ak5GenJetsSeq
      * process.primaryVertexFilter 
      * process.nEventsClean
      * getattr(process,"patPF2PATSequence"+postfix) # main PF2PAT 
      * process.nEventsPAT
      * process.flavorHistorySeq
    )

#    process.patSequenceComplete = cms.Sequence(
    #  + process.patDefaultSequence
    #  + process.patPFBRECOSequencePFlow
#        process.patPF2PATSequencePFlow
#      + process.nEventsPAT
#    )

    process.load("KrAFT.GeneratorTools.ttbar2bFilter_cfi")

    ## Defile paths
    if decayMode in ("all", "dilepton", "ElEl", "ee"):
        process.pElEl = cms.Path(
            process.nEventsTotal
          * process.commonFilterSequence
          * process.ttbar2bFilter
          * process.hltElEl * process.nEventsHLTElEl
        )
#        if doPAT: process.pElEl += process.patSequenceComplete
    if decayMode in ("all", "dilepton", "MuMu", "mumu"):
        process.pMuMu = cms.Path(
            process.nEventsTotal
          * process.commonFilterSequence
          * process.hltMuMu * process.nEventsHLTMuMu
        )
#        if doPAT: process.pMuMu += process.patSequenceComplete
    if decayMode in ("all", "dilepton", "MuEl", "emu"):
        process.pMuEl = cms.Path(
            process.nEventsTotal
          * process.commonFilterSequence
          * process.hltMuEl * process.nEventsHLTMuEl
        )
#        if doPAT: process.pMuEl += process.patSequenceComplete
    if decayMode in ("all", "MuJets"):
        process.selectedPatMuonsPFlow.cut = 'isPFMuon && (isGlobalMuon || isTrackerMuon) && pt > 10 && abs(eta) < 2.5 && (chargedHadronIso + max(0.,neutralHadronIso+photonIso-0.5*puChargedHadronIso))/pt < 0.15'
        process.selectedPatElectronsPFlow.cut = 'pt > 20 && abs(eta) < 2.5 && electronID("mvaTrigV0") > 0. && (chargedHadronIso + max(0.,neutralHadronIso+photonIso-0.5*puChargedHadronIso))/pt < 0.15'
        process.muonVetoFilter.maxNumber = 1

        process.pMuJets = cms.Path(
            process.nEventsTotal
          * process.commonFilterSequence
          * process.hltMuJets * process.nEventsHLTMuJets
        )
#        if doPAT: process.pMuJets += process.patSequenceComplete
        process.pMuJets *= process.muonVetoFilter
        process.pMuJets += process.electronVetoFilter

        if runOnMC: process.pMuJets.remove(process.hltMuJets)
    if decayMode in ("all", "ElJets"):
        process.selectedPatMuonsPFlow.cut = 'isPFMuon && (isGlobalMuon || isTrackerMuon) && pt > 10 && abs(eta) < 2.5 && (chargedHadronIso + max(0.,neutralHadronIso+photonIso-0.5*puChargedHadronIso))/pt < 0.15'
        process.selectedPatElectronsPFlow.cut = 'pt > 20 && abs(eta) < 2.5 && electronID("mvaTrigV0") > 0. && (chargedHadronIso + max(0.,neutralHadronIso+photonIso-0.5*puChargedHadronIso))/pt < 0.15'
        process.electronVetoFilter.maxNumber = 1

        process.pElJets = cms.Path(
            process.nEventsTotal
          * process.commonFilterSequence
          * process.hltElJets * process.nEventsHLTElJets
        )
#        if doPAT: process.pElJets += process.patSequenceComplete
        process.pElJets *= process.muonVetoFilter
        process.pElJets += process.electronVetoFilter

        if runOnMC: process.pElJets.remove(process.hltElJets)

    return process

def addNtupleStep(process, runOnMC):
    # Add ntuple production
    process.load("KrAFT.Configuration.ntuple_templateTTBB_cff")
    process.goodJets.isMC = runOnMC

    for mode in ('ElEl', 'MuMu', 'MuEl', 'ElJets', 'MuJets'):
        if not hasattr(process, 'p'+mode): continue

        getattr(process, mode).isMC = runOnMC
        p = getattr(process, 'p'+mode)
        ntupleStep = getattr(process, 'ntupleSequence'+mode)
        if not runOnMC:
            ntupleStep.remove(process.pdfWeight)
            ntupleStep.remove(process.pileupWeight)
        p += ntupleStep

        getattr(process, mode).eventCounters.extend([
            "nEventsHLT%s" % mode, "nEventsNtuple%s" % mode,
            "nEventsNtuple%s" % mode,
        ])

