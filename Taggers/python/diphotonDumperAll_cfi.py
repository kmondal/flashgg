import FWCore.ParameterSet.Config as cms

from flashgg.Taggers.diphotonDumpConfig_cff import diphotonDumpConfig

diphotonDumperAll = cms.EDAnalyzer('CutBasedDiPhotonDumper',
                                **diphotonDumpConfig.parameters_()
                                )
