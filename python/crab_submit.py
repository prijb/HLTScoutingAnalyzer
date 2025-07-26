
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'SingleMuonMonitor_v1'
config.General.workArea = 'ReHLTNtuples'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'analysis'
config.JobType.psetName = 'python/ntuplizer_cfg.py'
config.JobType.maxMemoryMB = 2500
config.JobType.pyCfgParams = ['numEvents=-1']

config.Data.inputDataset = '/EphemeralHLTPhysics3/ppradeep-HLTRerun_SingleMuonMonitor_v1-8e1bfffa23d6a4aff79676759cfd0adc/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/ppradeep/Data/'
config.Data.publication = True
config.Data.outputDatasetTag = 'SingleMuonMonitor_v1_ntuple'

config.Site.storageSite = 'T2_UK_London_IC'
