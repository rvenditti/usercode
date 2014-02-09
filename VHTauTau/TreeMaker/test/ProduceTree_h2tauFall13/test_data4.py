import FWCore.ParameterSet.Config as cms
process = cms.Process("HTauTauTree")
#------------------------
# Message Logger Settings
#------------------------
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
#--------------------------------------
# Event Source & # of Events to process
#---------------------------------------
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring()
)
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )
#-----------------------------
# Geometry
#-----------------------------
process.load("Configuration.Geometry.GeometryIdeal_cff")
#-----------------------------
# Magnetic Field
#-----------------------------
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#-------------
# Global Tag
#-------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'FT_53_V21_AN4::All'
#-------------
# Output ROOT file
#-------------
process.TFileService = cms.Service("TFileService",
  fileName = cms.string('/lustre/cms/store/user/rosma/UpgradePh1_FullSim_16July/GluGluToHToTauTau_Fall13_Tree/tree4.root')
)
#--------------------------------------------------
# VHTauTau Tree Specific
#--------------------------------------------------
process.load("VHTauTau.TreeMaker.TreeCreator_cfi")
process.load("VHTauTau.TreeMaker.TreeWriter_cfi")
process.load("VHTauTau.TreeMaker.TreeContentConfig_cff")
#-------------------------------------------------------
# PAT 
#------------------------------------------------------
process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
process.load("PhysicsTools.PatAlgos.patSequences_cff")
import PhysicsTools.PatAlgos.tools.tauTools as tauTools
tauTools.switchToPFTauHPS(process) # For HPS Taus
## --
## Switch on PAT trigger
## --
import PhysicsTools.PatAlgos.tools.trigTools as trigTools
trigTools.switchOnTrigger( process, outputModule='' ) # This is optional and can be omitted.
process.p = cms.Path(
  process.treeCreator +
  process.treeContentSequence +
  process.treeWriter
)
# List File names here
#---------------------------------------
process.PoolSource.fileNames = [
'file:/lustre/cms/store/user/rosma/GluGluToHToTauTau_M-125_13TeV_Fall13dr_PU20bx25_POSTLS162/patTuple_28_1_eHa.root',
'file:/lustre/cms/store/user/rosma/GluGluToHToTauTau_M-125_13TeV_Fall13dr_PU20bx25_POSTLS162/patTuple_29_1_vjS.root',
'file:/lustre/cms/store/user/rosma/GluGluToHToTauTau_M-125_13TeV_Fall13dr_PU20bx25_POSTLS162/patTuple_30_1_LU4.root',
'file:/lustre/cms/store/user/rosma/GluGluToHToTauTau_M-125_13TeV_Fall13dr_PU20bx25_POSTLS162/patTuple_31_1_9QC.root',
'file:/lustre/cms/store/user/rosma/GluGluToHToTauTau_M-125_13TeV_Fall13dr_PU20bx25_POSTLS162/patTuple_3_1_85w.root',
]
