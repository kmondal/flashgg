import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils

process = cms.Process("FLASHggMicroAOD")

process.load('HLTrigger.HLTfilters.hltHighLevel_cfi')
process.hltHighLevel.throw = cms.bool(False)
process.hltHighLevel.HLTPaths = ["HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon18_AND_HE10_R9Id65_Mass95_v*",
                                 "HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_v*"]

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_74_V8A'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10)

process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring(
"/store/relval/CMSSW_7_4_0_pre9/DoubleElectron/MINIAOD/GR_R_74_V8A_RelVal_zEl2011A-v1/00000/24BE70F0-FCD3-E411-8A77-0025905A6104.root",
"/store/relval/CMSSW_7_4_0_pre9/DoubleElectron/MINIAOD/GR_R_74_V8A_RelVal_zEl2011A-v1/00000/D64563F1-FCD3-E411-A7FA-0025905A6092.root",
        )
)

process.MessageLogger.cerr.threshold = 'ERROR'

process.load("flashgg/MicroAOD/flashggMicroAODZeeValidationSequence_cff")
process.load("flashgg/Taggers/flashggDiPhotonMVA_cfi")
process.flashggDiPhotonMVA.DiPhotonTag = cms.untracked.InputTag('flashggZeeDiPhotons')

#from flashgg.MicroAOD.flashggMicroAODOutputCommands_cff import microAODDebugOutputCommand
from flashgg.MicroAOD.flashggMicroAODZeeValidationOutputCommands_cff import microAODZeeValidationOutputCommand
process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string('ZeeMicroAODOutputFile2011.root'),
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
