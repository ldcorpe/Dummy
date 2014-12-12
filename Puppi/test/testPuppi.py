import FWCore.ParameterSet.Config as cms

process = cms.Process('TestPuppi')
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.GlobalTag.globaltag = 'START53_V7G::All'

process.load('Dummy/Puppi/Puppi_cff')   

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.source = cms.Source("PoolSource",
  #fileNames  = cms.untracked.vstring('file:/tmp/pharris/RSGravitonToWW_kMpl01_M_3000_Tune4C_13TeV_pythia8_PU_S14_PAT.root')
  fileNames  = cms.untracked.vstring('/store/mc/Phys14DR/QCD_Pt-15to3000_Tune4C_Flat_13TeV_pythia8/MINIAODSIM/PU20bx25_trkalmb_PHYS14_25_V1-v1/00000/1EF51024-986B-E411-A6F6-20CF300E9EAF.root')
)
process.source.inputCommands = cms.untracked.vstring("keep *",
                                                     "drop *_MEtoEDMConverter_*_*")

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  Rethrow     = cms.untracked.vstring('ProductNotFound'),
  fileMode    = cms.untracked.string('NOMERGE')
)


process.puppiSequence = cms.Sequence(process.puppi)
process.p = cms.Path(process.puppiSequence)
process.output = cms.OutputModule("PoolOutputModule",                                                                                                                                                     
                                  outputCommands = cms.untracked.vstring('keep *'),
                                  #outputCommands = cms.untracked.vstring('drop *','keep *_*_*_RECO','drop *_*_Cleaned_*','keep *_puppi_*_*'),
                                  fileName       = cms.untracked.string ("Output.root")                                                                                                                   
)
# schedule definition                                                                                                       
process.outpath  = cms.EndPath(process.output) 
