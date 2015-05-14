import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing

from flashgg.MetaData.samples_utils import SamplesManager

process = cms.Process("zeeValidationDumper")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR_R_74_V8::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring("file:/afs/cern.ch/work/k/kmondal/public/FLASHgg/CMSSW_7_4_0_pre9/src/flashgg/MicroAOD/test/2012/ZeeMicroAODOutputFile2012.root",)
)

import flashgg.Taggers.dumperConfigTools as cfgTools
process.load("flashgg.Taggers.diphotoMVAWithZeeDumper_cff")
process.load("flashgg.Taggers.photonDumper_cfi")

process.load("flashgg.Taggers.diphotonDumper_cfi")
process.diphotonDumper.src = cms.InputTag("flashggZeeDiPhotons")
process.diphotonDumperAll.src = cms.InputTag("flashggZeeDiPhotons")

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("zeeValidation.root"),
                                   closeFileFast = cms.untracked.bool(True)
)

#DIPHOTON MVA
cfgTools.addCategories(process.DiPhotonWithZeeMVADumper,
                       [("All","1", 0),],
                       variables=["dipho_mva:=mvaValue"],
                       histograms=["dipho_mva>>dipho_mva(100,-1,1)",]
)
# split tree, histogram and datasets by process
process.DiPhotonWithZeeMVADumper.nameTemplate ="zeevalidation_$SQRTS_$LABEL_$SUBCAT"

cfgTools.addCategories(process.diphotonDumperAll,
                       ## categories definition
                       ## cuts are applied in cascade. Events getting to these categories have already failed the "Reject" selection
                       [("All","1", 0)],
                       variables=["leadPhIso := leadingPhoton.pfPhoIso03",
                                  "subleadPhIso := subLeadingPhoton.pfPhoIso03",
                                  "sieie := leadingPhoton.full5x5_sigmaIetaIeta",
                                  "covieip := leadingPhoton.sieip",
                                  "etawidth := leadingPhoton.superCluster.etaWidth",
                                  "phiwidth := leadingPhoton.superCluster.phiWidth",
                                  "s4ratio := leadingPhoton.s4",
                                  "leadR9 := leadingPhoton.r9",
                                  "subleadR9 := subLeadingPhoton.r9",
                                  "scEta := leadingPhoton.superCluster.eta",
                                  "CMS_hgg_mass[120,60,120]:=mass", 
                                  "leadPt      :=leadingPhoton.pt",
                                  "subleadPt   :=subLeadingPhoton.pt",
                                  "minR9       :=min(leadingPhoton.r9,subLeadingPhoton.r9)",
                                  "maxEta      :=max(abs(leadingPhoton.superCluster.eta),abs(leadingPhoton.superCluster.eta))",
                                  "subIDMVA    := subLeadPhotonId",
                                  "leadIDMVA   := leadPhotonId",
                                  "vtxProb     := vtxProbMVA",
                                  "cosdphi     := cos(leadingPhoton.phi-subLeadingPhoton.phi)",
                                  "sigmaEoE1   := leadingPhoton.sigEOverE",
                                  "sigmaEoE2   := subLeadingPhoton.sigEOverE",
                                  "sigmaMoM    := .5*sqrt(leadingPhoton.sigEOverE*leadingPhoton.sigEOverE+subLeadingPhoton.sigEOverE*subLeadingPhoton.sigEOverE)",
                                  "leadChIsoRv := leadingPhoton.pfChgIsoWrtChosenVtx03",
                                  "subleadChIsoRv := subLeadingPhoton.pfChgIsoWrtChosenVtx03",
                                  ],
                       histograms=["leadPhIso>>leadPhIso(100, 0, 10)",
                                   "subleadPhIso>>subleadPhIso(100, 0, 10)",
                                   "sieie>>sieie(200, 0, 0.04)",
                                   "covieip>>covieip(100, -0.001, 0.001)",
                                   "etawidth>>etawidth(100, 0, 0.04)",
                                   "phiwidth>>phiwidth(120, 0, 0.12)",
                                   "s4ratio>>s4ratio(100, 0, 1)",
                                   "leadR9>>leadR9(110,0,1.1)",
                                   "subleadR9>>subleadR9(110,0,1.1)",
                                   "scEta>>sceta(50, -2.5, 2.5)",
                                   "CMS_hgg_mass>>mass(120,60,120)",
                                   "leadPt>>leadPt(100, 0, 100)",
                                   "subleadPt>>subleadPt(100, 0, 100)",
                                   "minR9>>minR9(110, 0, 1.1)",
                                   "maxEta>>maxEta(50, -2.5, 2.5)",
                                   "subIDMVA>>subidmva(90,-0.3, 0.6)",
                                   "leadIDMVA>>leadidmva(90, -0.3, 0.6)",
                                   "vtxProb>>vtxProb(200, -1, 1)",
                                   "cosdphi>>cosdphi(100, -1, 1)",
                                   "sigmaEoE1>>sigmaEoE1(100, 0, 0.1)",
                                   "sigmaEoE2>>sigmaEoE2(100, 0, 0.1)",
                                   "sigmaMoM>>sigmaMoM(100, 0, 0.1)",
                                   "leadChIsoRv>>leadChIsoRv(120, 0, 12)",
                                   "subleadChIsoRv>>subleadChIsoRv(120, 0, 12)",
                                   ]
                       )

process.diphotonDumperAll.nameTemplate ="zeevalidation_$SQRTS_$LABEL_$SUBCAT"

cfgTools.addCategories(process.diphotonDumper,
                       ## categories definition
                       ## cuts are applied in cascade. Events getting to these categories have already failed the "Reject" selection
                       [
                         ("EE", "abs(leadingPhoton.superCluster.eta)>1.566",0),
                         ("EB", "abs(leadingPhoton.superCluster.eta)<1.4442", 0),
                        ],
                       variables=["leadPhIso := leadingPhoton.pfPhoIso03",
                                  "subleadPhIso := subLeadingPhoton.pfPhoIso03",
                                  "sieie := leadingPhoton.full5x5_sigmaIetaIeta",
                                  "covieip := leadingPhoton.sieip",
                                  "etawidth := leadingPhoton.superCluster.etaWidth",
                                  "phiwidth := leadingPhoton.superCluster.phiWidth",
                                  "s4ratio := leadingPhoton.s4",
                                  "leadR9 := leadingPhoton.r9",
                                  "subleadR9 := subLeadingPhoton.r9",
                                  "scEta := leadingPhoton.superCluster.eta",
                                  "CMS_hgg_mass[120,60,120]:=mass",
                                  "leadPt      :=leadingPhoton.pt",
                                  "subleadPt   :=subLeadingPhoton.pt",
                                  "minR9       :=min(leadingPhoton.r9,subLeadingPhoton.r9)",
                                  "maxEta      :=max(abs(leadingPhoton.superCluster.eta),abs(leadingPhoton.superCluster.eta))",
                                  "subIDMVA    := subLeadPhotonId",
                                  "leadIDMVA   := leadPhotonId",
                                  "vtxProb     := vtxProbMVA",
                                  "cosdphi     := cos(leadingPhoton.phi-subLeadingPhoton.phi)",
                                  "sigmaEoE1   := leadingPhoton.sigEOverE",
                                  "sigmaEoE2   := subLeadingPhoton.sigEOverE",
                                  "sigmaMoM    := .5*sqrt(leadingPhoton.sigEOverE*leadingPhoton.sigEOverE+subLeadingPhoton.sigEOverE*subLeadingPhoton.sigEOverE)",
                                  "leadChIsoRv := leadingPhoton.pfChgIsoWrtChosenVtx03",
                                  "subleadChIsoRv := subLeadingPhoton.pfChgIsoWrtChosenVtx03",
                                  ],
                       histograms=["leadPhIso>>leadPhIso(100, 0, 10)",
                                   "subleadPhIso>>subleadPhIso(100, 0, 10)",
                                   "sieie>>sieie(200, 0, 0.04)",
                                   "covieip>>covieip(100, -0.001, 0.001)",
                                   "etawidth>>etawidth(100, 0, 0.04)",
                                   "phiwidth>>phiwidth(120, 0, 0.12)",
                                   "s4ratio>>s4ratio(100, 0, 1)",
                                   "leadR9>>leadR9(110,0,1.1)",
                                   "subleadR9>>subleadR9(110,0,1.1)",
                                   "scEta>>sceta(50, -2.5, 2.5)",
                                   "CMS_hgg_mass>>mass(120,60,120)",  
                                   "leadPt>>leadPt(100, 0, 100)",
                                   "subleadPt>>subleadPt(100, 0, 100)",
                                   "minR9>>minR9(110, 0, 1.1)",
                                   "maxEta>>maxEta(50, -2.5, 2.5)",
                                   "subIDMVA>>subidmva(90,-0.3, 0.6)",
                                   "leadIDMVA>>leadidmva(90, -0.3, 0.6)",
                                   "vtxProb>>vtxProb(200, -1, 1)",
                                   "cosdphi>>cosdphi(100, -1, 1)",
                                   "sigmaEoE1>>sigmaEoE1(100, 0, 0.1)",
                                   "sigmaEoE2>>sigmaEoE2(100, 0, 0.1)",
                                   "sigmaMoM>>sigmaMoM(100, 0, 0.1)",
                                   "leadChIsoRv>>leadChIsoRv(120, 0, 12)",
                                   "subleadChIsoRv>>subleadChIsoRv(120, 0, 12)",
                                  ]
                       )

process.diphotonDumper.nameTemplate ="zeevalidation_$SQRTS_$LABEL_$SUBCAT"

from flashgg.MetaData.JobConfig import customize
customize.setDefault("maxEvents",-1)
customize.setDefault("targetLumi",1.e+4)
customize(process)

process.p = cms.Path(process.DiPhotonWithZeeMVADumper*process.diphotonDumperAll*process.diphotonDumper)
