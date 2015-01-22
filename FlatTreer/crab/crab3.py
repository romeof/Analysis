from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'TTbHWW_37'
config.General.workArea = 'TTbHWW_37'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/f/fromeo/CMSSW_7_0_6_patch1/src/Test/MiniAnalyzer/python/Config_miniaod_flatree.py'

config.section_("Data")
config.Data.inputDataset = '/TTbarH_HToWWTo2L2Nu_M-125_13TeV_amcatnlo-pythia8-tauola/Spring14miniaod-141029_PU40bx50_PLS170_V6AN2-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = 100
config.Data.outLFN = '/store/user/fromeo/'

config.section_("Site")
config.Site.storageSite = 'T2_IT_Pisa' 
