import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer16_ID_ISO_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v6', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v10', '')#mc

process.MessageLogger.cerr.FwkReport.reportEvery = 200
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1),SkipEvent = cms.untracked.vstring('ProductNotFound'))
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(       
 #'/store/data/Run2017B/Charmonium/MINIAOD/PromptReco-v1/000/297/050/00000/183B4680-3356-E711-B33A-02163E014487.root',
 #'/store/data/Run2017F/Charmonium/MINIAOD/17Nov2017-v1/710000/F2F520DE-A8F4-E711-AB19-0025905D1D50.root',#char working
 #'/store/data/Run2017C/SingleElectron/MINIAOD/17Nov2017-v1/80000/1A108E38-C2FE-E711-8A04-F01FAFD67D07.root' # single Ele
 #'/store/mc/RunIIFall17MiniAODv2/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext2-v1/70000/EC03D597-4A9F-E911-B70C-00259021A4A2.root'#ZZmc
 #'/store/mc/RunIIFall17MiniAODv2/GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUGenV7011_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/60000/F80D89CA-1C72-E811-963C-FA163EB264B5.root'#HZZ MC
 'file:HZJPsi_MINIAOD-prod_PAT_4.root'
 )
)
#Prueba de fuego.
#process.load("myAnalyzers.JPsiKsPAT.miniAODmuonsRootupler_cfi")
#process.rootuple.dimuons = cms.InputTag('miniaodPATMuonsWithTrigger') 

process.rootuple = cms.EDAnalyzer('miniAODmuons',
                          dimuons = cms.InputTag("slimmedMuons"),
                          #dimuons = cms.InputTag("selectedPatMuons"),
                          #dielectron = cms.InputTag("selectedPatElectrons"),
                          dielectron = cms.InputTag("slimmedElectrons"),
                          Trak = cms.InputTag("packedPFCandidates"),
                          primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                          bits = cms.InputTag("TriggerResults::HLT"),
                          objects = cms.InputTag("slimmedPatTrigger"),
                          #packed = cms.InputTag("packedGenParticles"),
                          pruned = cms.InputTag("prunedGenParticles"),
                          #genParticles=cms.InputTag("genParticles")
                          #objects = cms.InputTag("selectedPatTrigger"),
                          #prescales = cms.InputTag("patTrigger"),
                          #eleLooseIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
                          isMC = cms.bool(True),
                          )

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('Rootuple_Jpsi_2017-MiniAOD_HZJPsi_MCAll_4.root'),
)

process.p = cms.Path(process.rootuple)


