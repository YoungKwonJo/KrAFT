import FWCore.ParameterSet.Config as cms
runOnMC = True

from KrAFT.Configuration.customiseTTBB_cff import *
process = initialise(decayMode="dilepton", runOnMC=runOnMC)
addNtupleStep(process, runOnMC=runOnMC)

import os
hostName = os.environ['HOSTNAME']
if 'cern.ch' in hostName:
    process.source.fileNames = [
'/store/temp/user/youngjo/ttbar_01.root',
'/store/temp/user/youngjo/ttbar_02.root',
#'/store/temp/user/youngjo/zjets_01.root',
#'/store/temp/user/youngjo/zjets_02.root',
#'/store/temp/user/youngjo/zjets_03.root',
#'/store/temp/user/youngjo/zjets_04.root',
#'/store/temp/user/youngjo/zjets_05.root',
    ]
elif 'uos.ac.kr' in hostName:
    process.source.fileNames = [
        '/store/mc/Summer12_DR53X/TTJets_FullLeptMGDecays_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7C-v2/10000/FECC62BD-3898-E211-82F9-003048FFD7BE.root'
    ]

process.maxEvents.input = 100

