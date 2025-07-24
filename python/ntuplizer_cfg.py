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
process.GlobalTag = GlobalTag(process.GlobalTag, '150X_dataRun3_HLT_v1', '')

# File service
process.TFileService = cms.Service("TFileService",
  fileName = cms.string("output.root")
)

process.ntuplizer = cms.EDAnalyzer(
  'Ntuplizer',
  triggerResults = cms.InputTag("TriggerResults", "", "HLTX"),
  muonsNoVtx = cms.InputTag("hltScoutingMuonPackerNoVtx", "", "HLTX"),
  muonsVtx = cms.InputTag("hltScoutingMuonPackerVtx", "", "HLTX"),
  PV = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx", "HLTX"),
  SVNoVtx = cms.InputTag("hltScoutingMuonPackerNoVtx", "displacedVtx", "HLTX"),
  SVVtx = cms.InputTag("hltScoutingMuonPackerVtx", "displacedVtx", "HLTX"),
)

process.p = cms.Path(process.ntuplizer)