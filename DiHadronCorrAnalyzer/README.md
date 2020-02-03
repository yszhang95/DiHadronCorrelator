# DiHadronCorrAnalyzer
## General Descriptions
It is an analyzer for two-particle correlations in small systems. 
It will correlate the tracks associated to each primary vertex.
By default all the vertices provided by config file would be read.
The output of the analzyer are a collections of histograms, which
are organized by the pT of trigger particles and Ntrkoffline of the
corresponding vertices. 

## Several variables used in the analysis
- Primary tracks: tracks fulfill the following criteria
    * zDCA and xyDCA w.r.t. the corresponding primary vertex < 3cm
    * pTerr/pT < 0.1 
    * high purity
    * pT >= 0.3 GeV and |eta| < 2.4
- Ntrkoffline: it is calculated by counting the
  number of primary tracks with pT larger than 0.4 GeV. 
  It is classfied by the parameter `nTrkOfflineBinEdges`.
- Trigger particles, primary tracks within some pT ranges. They are clustered 
  by the parameter `ptBinEdges`.
- Associated particles, primary tracks with pT between 0.3-3.0 GeV.
- Signal, pairs correlated in the same primary vertex
- Background, pairs correlated by mixing the events. Each trigger particle
  is correlated with the associated particles in other `nMixEvts` track collectioins,
  whose corresponding vertices have similar positions (`dz<2cm`) and Ntrkoffline.

## Default configuration
Default configuration can be found in 
`DiHadronCorrelator/DiHadronCorrAnalyzer/python/dihadroncorrelator_cfi.py`.

Configuration fragments can be found in
`DiHadronCorrelator/DiHadronCorrAnalyzer/python/dihadroncorrelator_cff.py`.

Parameters in `dihadroncorrelator_cfi.py` are listed below.

|name               |Default parameters                                                |
|-------------------|------------------------------------------------------------------|
|tracks             |generalTracks                                                     |
|primaryVertices    |offlinePrimaryVertices                                            |
|nTrkOfflineBinEdges|{0, 20, 40, 80, 100, 10000000}                                    |
|ptBinEdges         |{0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3} #GeV             |
|zVtxBinEdges       |{-15, -13, -11, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9, 11, 13, 15} #cm|
|nMixEvts           |10                                                                |


# myVtxObj
## General Descriptions
It is in 
`DiHadronCorrelator/DiHadronCorrAnalyzer/plugins/myVtxObj.*`. It defines serverl 
user-defined classes and template functions. One can check the source code for more details.
