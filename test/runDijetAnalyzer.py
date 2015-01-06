import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

###############################
####### Parameters ############
###############################

options = VarParsing ('python')

options.register('reportEvery', 1000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=1000)"
)
options.register('outputFilename', 'histos.root',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Output file name"
)
options.register('dijetMassCut', 1118.,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.float,
    "DIjet mass cut (default is 1118 GeV)"
)
options.register('process', 'QCD',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "MC-simulated event type"
)
options.register('wantSummary', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Print out trigger and timing summary"
)

## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 50000)

options.parseArguments()

process = cms.Process("ANA")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # /QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/Phys14DR-PU20bx25_trkalmb_PHYS14_25_V1-v1/MINIAODSIM
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/5AE5B0FC-986B-E411-ACED-20CF3027A57B.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/1020E374-B26B-E411-8F91-E0CB4E29C513.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/1EF51024-986B-E411-A6F6-20CF300E9EAF.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/6C6C572A-9A6B-E411-93F1-00259073E364.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/7A60863A-9A6B-E411-A19D-002590D0AF6E.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/7EADDD13-9B6B-E411-9AC4-E0CB4E29C4F7.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/984A6622-9D6B-E411-849F-E0CB4E1A1149.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/C0587C5C-996B-E411-8388-20CF305B04D2.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/DA3D4205-9A6B-E411-AA7F-20CF3027A57B.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/EE72ECF8-996B-E411-B541-20CF305B057C.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/10000/307D472E-9A6B-E411-BC68-20CF3027A629.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/10000/389F3552-9A6B-E411-A8A3-20CF3019DEE8.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/10000/40BAFDD0-996B-E411-B886-002590D0AFEE.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/10000/50C8205D-9A6B-E411-885D-E0CB4EA0A917.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/10000/702436D6-996B-E411-8762-E0CB4E29C4D8.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/10000/82DD6BBC-996B-E411-ADD5-0025907B4FB6.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/10000/9CB0BB2E-9A6B-E411-8F3E-0025907B503A.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/10000/AEE02092-996B-E411-9A66-00259073E3C8.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/10000/C4A542EF-996B-E411-9F22-002590D0B098.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/10000/EEB3F427-9A6B-E411-8C4A-20CF305616F4.root',
        '/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/10000/F267CC51-9A6B-E411-9089-00259073E4E8.root'
    )
)

if options.process == "Qstar":
    process.source.fileNames = [
        # /QstarToJJ_M_4000_Tune4C_13TeV_pythia8/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM
        '/store/mc/Phys14DR/QstarToJJ_M_4000_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/8C78B451-576C-E411-89E0-002481E15204.root',
        '/store/mc/Phys14DR/QstarToJJ_M_4000_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/A496733E-576C-E411-ABE0-002590DB9214.root',
        '/store/mc/Phys14DR/QstarToJJ_M_4000_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/26D3D566-546C-E411-A4C7-00237DDCBEA4.root',
        '/store/mc/Phys14DR/QstarToJJ_M_4000_Tune4C_13TeV_pythia8/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/10000/BA78AEA4-546C-E411-8E2B-002590DB9214.root'
    ]

## Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string(options.outputFilename)
)

## Options and Output Report
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
    allowUnscheduled = cms.untracked.bool(True)
)

## Initialize the analyzer
process.dijetAna = cms.EDAnalyzer('DijetAnalyzer',
   jetSrc       = cms.InputTag("slimmedJets"),
   vertexSrc    = cms.InputTag("offlineSlimmedPrimaryVertices"),
   rhoSrc       = cms.InputTag("fixedGridRhoFastjetAll"),
   JESbias      = cms.double(1.0),
   dijetMassCut = cms.double(options.dijetMassCut),
   L1corr       = cms.FileInPath('CMSDAS2015/EXODijetResonances/data/PHYS14_25_V2_L1FastJet_AK4PFchs.txt'),
   L2corr       = cms.FileInPath('CMSDAS2015/EXODijetResonances/data/PHYS14_25_V2_L2Relative_AK4PFchs.txt'),
   L3corr       = cms.FileInPath('CMSDAS2015/EXODijetResonances/data/PHYS14_25_V2_L3Absolute_AK4PFchs.txt')
)

## Let it run
process.p = cms.Path(process.dijetAna)
