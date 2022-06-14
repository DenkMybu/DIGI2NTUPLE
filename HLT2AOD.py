# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: --python_filename HLT2AOD.py --eventcontent AODSIM --customise Configuration/DataProcessing/Utils.addMonitoring --datatier AODSIM --fileout file:HSCP_Gluino_Mass1800_AOD.root --conditions 123X_mcRun3_2021_realistic_v4 --step RAW2DIGI,L1Reco,RECO,RECOSIM,EI --geometry DB:Extended --filein file:hltOutput.root --era Run3 --runUnscheduled --no_exec --mc -n 10
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('RECO',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/opt/sbg/cms/ui2_data1/rhaeberl/Prod/prodMarch2022_Gluino1800/HLT/CaloTowers/0000/HSCP_Gluino_Mass1800_HLT_135.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('--python_filename nevts:1000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)
# Custom Keep Drop
process.CustomKeep = cms.PSet(
    outputCommands = cms.untracked.vstring(
    'keep *_hltTowerMakerForAll_*_*',
    'keep GlobalAlgBlk_hltGtStage2Digis_*_HLT',
    'keep GlobalExtBlk_hltGtStage2Digis_*_HLT',
    'keep *_hltMet_*_HLT',
    'drop TriggerFilterObjectWithRefs_*_*_*'
    ),
    eventAutoFlushCompressedSize=cms.untracked.int32(30*1024*1024),
    compressionAlgorithm=cms.untracked.string("LZMA"),
    compressionLevel=cms.untracked.int32(4),
)



# Output definition

process.AODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('AODSIM'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    #fileName = cms.untracked.string('file:HSCP_Gluino_Mass1800_AOD_towers.root'),
    fileName = cms.untracked.string('file:HSCP_Gluino_Mass1800_AOD_CT_MET.root'),
    outputCommands = process.AODSIMEventContent.outputCommands 
)
process.AODSIMoutput.outputCommands.extend(process.CustomKeep.outputCommands)

'''
process.AODSIMoutput.outputCommands.extend(
    'keep GlobalAlbGlk_hltGtStage2Digis*_*_*',
    'keep GlobalExtBlk_hltGtStage2Digis*_*_*',
    'keep EventAux_*_*_*',
    'keep LumiSummary_*_*_*',
    'keep edmMergeableCounter_*_*_*',
    'keep GenRunInfoProduct_*_*_*',
    'keep GenEventInfoProduct_generator_*_*',
    'keep *_genParticlesSkimmed_*_*',
    'keep *_genParticlePlusGeant_*_*',
    'keep *_offlinePrimaryVertices_*_*',
    'keep recoTracks_generalTracks_*_*',
    'keep recoTracks_standAloneMuons_*_*',
    'me = cms.untracked.string('file:HSCP_Gluino_Mass1800_AOD_CT_MET.root'),keep recoTrackExtras_standAloneMuons_*_*',
    'keep TrackingRecHitsOwned_standAloneMuons_*_*',
    'keep recoTracks_globalMuons_*_*',
    'keep recoTrackExtras_globalMuons_*_*',
    'keep recoMuons_muons_*_*',
    'keep recoMuonTimeExtraedmValueMap_muons_*_*',
    'keep edmTriggerResults_TriggerResults_*_* ',
    'keep *_ak4PFJetsCHS__*',
    'keep recoPFMETs_pfMet__*',
    'keep *_HSCParticleProducer_*_*',
    'keep *_HSCPIsolation*_*_*',
    'keep *_dedxHitInfo*_*_*',
    'keep triggerTriggerEvent_hltTriggerSummaryAOD_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_MuonSegmentProducer_*_*',
    'keep *_g4SimHits_StoppedParticles*_*',
    'keep PileupSummaryInfos_addPileupInfo_*_*',
    'keep *_dt4DSegments__*',
    'keep *_cscSegments__*',
    'keep *_hltTowerMakerEcal_*_*',
    'keep *_hltGtStage2Digis_Muon_HLT',
    'keep *_*_EcalRecHitsEB_*',
    'keep *_*_EcalRecHitsEE_*',
    'keep *_rpcRecHits_*_*'
)
'''



# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '123X_mcRun3_2021_realistic_v4', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
process.eventinterpretaion_step = cms.Path(process.EIsequence)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.AODSIMoutput_step = cms.EndPath(process.AODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,process.eventinterpretaion_step,process.endjob_step,process.AODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring 

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

#'keep *_hltTowerMakerForAll_*_HLT',
#'keep *_genParticles'
#'drop *_*_QualityMask_HLT',
#'drop TriggerFilterObjectWithRefs_*_*_*',
#'drop vector<reco::Vertex>_*_*_*'
#'drop vector<reco::TrackExtra>_*_*_*',
#'drop vector<reco::PFTau>_*_*_*',
#'drop *_*_MVAValues_HLT'
#'drop OwnVector<TrackingRecHit,edm::ClonePolicy<TrackingRecHit>>_*_*_*',
#'drop ContainerMask<edmNew::DetSetVector<SiStripCluster>>_*_*_*',
