import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("FLASHggMicroAOD")

process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.throw = cms.bool(False)
process.hltHighLevel.HLTPaths = ["HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v*",
                                 "HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v*"]

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_74_V8::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10)

process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
"/store/relval/CMSSW_7_4_0_pre9/DoubleElectron/MINIAOD/GR_R_74_V8A_RelVal_zEl2012D-v1/00000/5A04EF0A-29D4-E411-BB12-003048FFCC2C.root",
"/store/relval/CMSSW_7_4_0_pre9/DoubleElectron/MINIAOD/GR_R_74_V8A_RelVal_zEl2012D-v1/00000/A6D4F50A-29D4-E411-98A2-002618943838.root",
"/store/relval/CMSSW_7_4_0_pre9/DoubleElectron/MINIAOD/GR_R_74_V8A_RelVal_zEl2012D-v1/00000/DE947AB8-FED3-E411-995F-0025905AA9F0.root",
        )
)

process.MessageLogger.cerr.threshold = 'ERROR'

process.load("flashgg/MicroAOD/flashggMicroAODZeeValidationSequence_cff")
process.load("flashgg/Taggers/flashggDiPhotonMVA_cfi")
process.flashggDiPhotonMVA.DiPhotonTag = cms.untracked.InputTag('flashggZeeDiPhotons')

#from flashgg.MicroAOD.flashggMicroAODOutputCommands_cff import microAODDebugOutputCommand
from flashgg.MicroAOD.flashggMicroAODZeeValidationOutputCommands_cff import microAODZeeValidationOutputCommand
process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('ZeeMicroAODOutputFile2012.root'),
                               SelectEvents = cms.untracked.PSet(
                                SelectEvents = cms.vstring('p1')
                               ),
                               outputCommands = microAODZeeValidationOutputCommand
                               )
#process.out.outputCommands += microAODDebugOutputCommand # extra items for debugging, CURRENTLY REQUIRED

process.p1 = cms.Path(process.hltHighLevel*process.flashggMicroAODZeeValidationSequence*process.flashggDiPhotonMVA)
process.e = cms.EndPath(process.out)

from flashgg.MicroAOD.MicroAODCustomize import customize
customize(process)
