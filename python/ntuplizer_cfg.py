import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
params = VarParsing('analysis')

params.register(
  'inputFile',
  'file:',
  VarParsing.multiplicity.singleton,
  VarParsing.varType.string,
  'Input file'
)
params.register(
  'numEvents',
  1000,
  VarParsing.multiplicity.singleton,
  VarParsing.varType.int,
  'Number of events to process'
)
params.register(
  'globalTag',
  '150X_dataRun3_HLT_v1',
  VarParsing.multiplicity.singleton,
  VarParsing.varType.string,
  'GlobalTag'
)
params.register(
  'process',
  'HLT',
  VarParsing.multiplicity.singleton,
  VarParsing.varType.string,
  'Process name'
)

process = cms.Process('Ntuplizer')
params.parseArguments()

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet(
  allowUnscheduled = cms.untracked.bool(True),
  wantSummary      = cms.untracked.bool(True),
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(params.numEvents))

process.source = cms.Source(
  "PoolSource", 
  fileNames = cms.untracked.vstring(params.inputFile)
)

# Transient track
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.MagneticField_cff')

# Global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '150X_dataRun3_HLT_v1', '')
process.GlobalTag = GlobalTag(process.GlobalTag, params.globalTag, '')

# File service
process.TFileService = cms.Service("TFileService",
  fileName = cms.string("output.root")
)

# 2024 paths
HLTPaths = ["DST_PFScouting_DoubleMuon_v", "DST_PFScouting_DoubleEG_v", "DST_PFScouting_JetHT_v", "DST_PFScouting_AXONominal_v", "DST_PFScouting_AXOTight_v", "DST_PFScouting_AXOVTight_v", "DST_PFScouting_SingleMuon_v", "DST_PFScouting_SinglePhotonEB_v", "DST_PFScouting_ZeroBias_v"]
L1Seeds = ["L1_SingleMu5_BMTF", "L1_SingleMu11_SQ14_BMTF", "L1_SingleMu13_SQ14_BMTF"]

process.load("EventFilter.L1TRawToDigi.gtStage2Digis_cfi")
process.gtStage2Digis.InputLabel = cms.InputTag( "hltFEDSelectorL1" , "", params.process)

# Pruning gen particles
# MINI
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
from PhysicsTools.PatAlgos.slimming.prunedGenParticles_cfi import prunedGenParticles
process.prunedGenParticles = prunedGenParticles.clone()
process.prunedGenParticles.src = cms.InputTag("genParticles", "", params.process)
# NANO
from PhysicsTools.NanoAOD.genparticles_cff import finalGenParticles
process.finalGenParticles = finalGenParticles.clone()
process.finalGenParticles.src = "prunedGenParticles"

process.ntuplizer = cms.EDAnalyzer(
  'Ntuplizer',
  triggerResults = cms.InputTag("TriggerResults", "", params.process),
  muonsNoVtx = cms.InputTag("hltScoutingMuonPackerNoVtx", "", params.process),
  muonsVtx = cms.InputTag("hltScoutingMuonPackerVtx", "", params.process),
  PV = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx", params.process),
  SVNoVtx = cms.InputTag("hltScoutingMuonPackerNoVtx", "displacedVtx", params.process),
  SVVtx = cms.InputTag("hltScoutingMuonPackerVtx", "displacedVtx", params.process),
  genParticles = cms.InputTag("prunedGenParticles"),
  hltPaths = cms.vstring(HLTPaths),
  doL1 = cms.bool(True),
  l1Seeds = cms.vstring(L1Seeds),
  AlgInputTag = cms.InputTag("gtStage2Digis"),
  l1tAlgBlkInputTag = cms.InputTag("gtStage2Digis"),
  l1tExtBlkInputTag = cms.InputTag("gtStage2Digis"),
  ReadPrescalesFromFile = cms.bool( False ),
)

#process.p = cms.Path(process.gtStage2Digis + process.ntuplizer)
process.p = cms.Path(process.gtStage2Digis + process.prunedGenParticles + process.ntuplizer)
#process.p = cms.Path(process.gtStage2Digis + process.prunedGenParticles + process.finalGenParticles + process.ntuplizer)