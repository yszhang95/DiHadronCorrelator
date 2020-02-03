import FWCore.ParameterSet.Config as cms
from DiHadronCorrelator.DiHadronCorrAnalyzer.dihadroncorrelator_cfi import *
dihadroncorrelator_new = dihadroncorrelator.clone()
#dihadroncorrelator_new.ptBinEdges = cms.vdouble(0., 100000.)
#dihadroncorrelator_new.nTrkOfflineBinEdges= cms.vuint32(0, 100000)
