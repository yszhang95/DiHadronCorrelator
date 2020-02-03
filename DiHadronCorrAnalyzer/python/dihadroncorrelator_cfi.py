import FWCore.ParameterSet.Config as cms

dihadroncorrelator = cms.EDAnalyzer('DiHadronCorrAnalyzer',
  tracks = cms.untracked.InputTag('generalTracks'),
  primaryVertices = cms.untracked.InputTag('offlinePrimaryVertices'),
  nTrkOfflineBinEdges = cms.vuint32(
    0,
    20,
    40,
    80,
    120,
    10000000
  ),
  ptBinEdges = cms.vdouble(
    0.3,
    0.6,
    0.9,
    1.2,
    1.5,
    1.8,
    2.1,
    2.4,
    2.7,
    3
  ),
  zVtxBinEdges = cms.vdouble(
    -15,
    -13,
    -11,
    -9,
    -7,
    -5,
    -3,
    -1,
    1,
    3,
    5,
    7,
    9,
    11,
    13,
    15
  ),
  nMixEvts = cms.uint32(10)
)
