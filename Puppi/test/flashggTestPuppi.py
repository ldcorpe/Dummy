import FWCore.ParameterSet.Config as cms

process = cms.Process('TestPuppi')
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.GlobalTag.globaltag = 'START53_V7G::All'

#process.load('Dummy/Puppi/Puppi_cff')   
#process.load('RecoJets/JetProducers/ak4PFJets_cfi')   

from RecoJets.JetProducers.AnomalousCellParameters_cfi import *

puppiCentral = cms.VPSet(
		cms.PSet(
			algoId           = cms.untracked.int32(5),  #0 is default Puppi
			useCharged       = cms.untracked.bool(True),
			applyLowPUCorr   = cms.untracked.bool(True),
			combOpt          = cms.untracked.int32(0),
			cone             = cms.untracked.double(0.3),
			rmsPtMin         = cms.untracked.double(0.1),
			rmsScaleFactor   = cms.untracked.double(1.0)
			)
		)

puppiForward = cms.VPSet(
		cms.PSet(
			algoId         = cms.untracked.int32(5),  #0 is default Puppi
			useCharged     = cms.untracked.bool(False),
			applyLowPUCorr = cms.untracked.bool(True),
			combOpt        = cms.untracked.int32(0),
			cone           = cms.untracked.double(0.3),
			rmsPtMin       = cms.untracked.double(0.5),
			rmsScaleFactor = cms.untracked.double(1.0)
			)
		)

process.puppi = cms.EDProducer("PuppiProducer",
		PuppiName      = cms.untracked.string("Puppi"),
		UseDeltaZCut   = cms.untracked.bool  (True),
		DeltaZCut      = cms.untracked.double(0.2),
		#candName       = cms.untracked.string('particleFlow'),
		#vertexName     = cms.untracked.string('offlinePrimaryVertices'),
		candName      = cms.untracked.string('packedPFCandidates'),
		vertexName     = cms.untracked.string('offlineSlimmedPrimaryVertices'),
		applyCHS       = cms.untracked.bool  (True),
		useExp         = cms.untracked.bool  (False),
		MinPuppiWeight = cms.untracked.double(0.01),
		algos          = cms.VPSet( 
			cms.PSet( 
				etaMin = cms.untracked.double(0.),
				etaMax = cms.untracked.double( 2.5),
				ptMin  = cms.untracked.double(0.),
				MinNeutralPt   = cms.untracked.double(0.2),
				MinNeutralPtSlope   = cms.untracked.double(0.02),
				puppiAlgos = puppiCentral
				),
			cms.PSet( 
				etaMin = cms.untracked.double(2.5),
				etaMax = cms.untracked.double(3.0),
				ptMin  = cms.untracked.double(0.0),
				MinNeutralPt        = cms.untracked.double(1.0),
				MinNeutralPtSlope   = cms.untracked.double(0.005),
				puppiAlgos = puppiForward
				),
			cms.PSet( 
				etaMin = cms.untracked.double(3.0),
				etaMax = cms.untracked.double(10.0),
				ptMin  = cms.untracked.double(0.0),
				MinNeutralPt        = cms.untracked.double(1.5),
				MinNeutralPtSlope   = cms.untracked.double(0.005),
				puppiAlgos = puppiForward
				)

			)
			)

PFJetParameters = cms.PSet(
					srcPVs = cms.InputTag(''),
					jetType = cms.string('PFJet'),
					doOutputJets = cms.bool(True),
					jetPtMin = cms.double(3.0),
					inputEMin = cms.double(0.0),
					inputEtMin = cms.double(0.0),
					doPVCorrection = cms.bool(False),
# pileup with offset correction
					doPUOffsetCorr = cms.bool(False),
# if pileup is false, these are not read:
					nSigmaPU = cms.double(1.0),
					radiusPU = cms.double(0.5),
# fastjet-style pileup
					doAreaFastjet = cms.bool( False),
					doRhoFastjet = cms.bool( False),
					doAreaDiskApprox = cms.bool( False),
					Active_Area_Repeats = cms.int32( 1),
					GhostArea = cms.double(0.01),
					Ghost_EtaMax = cms.double( 5.0),
					Rho_EtaMax = cms.double( 4.4),
					voronoiRfact = cms.double(-0.9),
					useDeterministicSeed= cms.bool( True ),
					minSeed = cms.uint32( 14327 )
	)

process.ak4PFJets = cms.EDProducer(
			"FastjetJetProducer",
			PFJetParameters,
			AnomalousCellParameters,
			src = cms.InputTag('puppi','Puppi'),
			jetAlgorithm = cms.string("AntiKt"),
			rParam = cms.double(0.4)
			)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
#process.source = cms.Source("PoolSource",
		#fileNames  = cms.untracked.vstring('file:/tmp/pharris/RSGravitonToWW_kMpl01_M_3000_Tune4C_13TeV_pythia8_PU_S14_PAT.root')
		#fileNames  = cms.untracked.vstring('root://cmsxrootd-site.fnal.gov//store/results/top/StoreResults/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/USER/Spring14dr_PU_S14_POSTLS170_V6AN1_miniAOD706p1_814812ec83fce2f620905d2bb30e9100-v2/00000/0012F41F-FA17-E411-A1FF-0025905A48B2.root')
#)
process.source = cms.Source("PoolSource",fileNames=cms.untracked.vstring("/store/cmst3/user/gpetrucc/miniAOD/v1/GluGluToHToGG_M-125_13TeV-powheg-pythia6_Flat20to50_PAT.root"))
process.source.inputCommands = cms.untracked.vstring("keep *",
			"drop *_MEtoEDMConverter_*_*")

process.options = cms.untracked.PSet(
			wantSummary = cms.untracked.bool(True),
			Rethrow     = cms.untracked.vstring('ProductNotFound'),
			fileMode    = cms.untracked.string('NOMERGE')
			)


#process.puppiSequence = cms.Sequence(process.puppi)
#process.jetSequence = cms.Sequence(process.ak4PFJets)
process.p = cms.Path(process.puppi*
			process.ak4PFJets	
			)
process.output = cms.OutputModule("PoolOutputModule",                                                                                                                                                     
			#outputCommands = cms.untracked.vstring('drop *','keep *_*_*_RECO','drop *_*_Cleaned_*','keep *_puppi_*_*'),
			outputCommands = cms.untracked.vstring('keep *'),
			fileName       = cms.untracked.string ("Output.root")                                                                                                                   
			)
# schedule definition                                                                                                       
process.outpath  = cms.EndPath(process.output) 
