# v00-05-03

* 2022-06-29 Andre Sailer ([PR#21](https://github.com/iLCSoft/MarlinFastJet/pull/21))
  - FastJetUtil: fix memory leak in clusterJets function. Change signature of this function to include the clusterSequence, fixes #20 
  - FastJetProcessor: fix issue for kt algorithm ,fixes #15

* 2021-12-03 Thomas Madlener ([PR#18](https://github.com/iLCSoft/MarlinFastJet/pull/18))
  - Migrate CI to github actions.

* 2021-12-03 Frank Gaede ([PR#17](https://github.com/iLCSoft/MarlinFastJet/pull/17))
  - fix order of fastjet libraries at linking step - needed on Ubuntu systems
  -  fixes https://github.com/iLCSoft/iLCInstall/issues/128

* 2021-12-03 Frank Gaede ([PR#16](https://github.com/iLCSoft/MarlinFastJet/pull/16))
  - minor fix of source paths in cmake file (for newer cmake versions, eg. 3.17)

# v00-05-02

* 2019-02-19 Andre Sailer ([PR#13](https://github.com/iLCSoft/MarlinFastJet/pull/13))
  - FastJetUtil: fix memory leak when creating jetDefinition with plugin (SiSCone, Valencia), needs #12 
  - MarlinFastJet: remove support for FastJet version 2

* 2018-05-18 Lars Rickard Strom ([PR#10](https://github.com/iLCSoft/MarlinFastJet/pull/10))
  - Implemented substructure parameters commonly used for top tagging. The definition of these new variables can be adjusted from the steering file (added options to steer how to calculate the energy correlation function (energyCorrelator) and how to calculate NSubjettiness (axesMode and measureMode). Only recommended options are implemented. The beta parameter would typically be the same as used for jet clustering but can also be varied separately (beta weights the angular distances between the jet constituents compared to their pt in the calculation of the energy correlation function and subjettiness). 
  - The substructure variables are added to the file as a collection "TopTaggerSubStructure" of seven elements (order: C2, D2, C3, D3, tau1, tau2, tau3). 
  - This update is backwards compatible, with the only difference being the addition of this collection.

# v00-05-01

# v00-05

* 2017-05-08 Lars Rickard Strom ([PR#9](https://github.com/iLCSoft/MarlinFastJet/pull/9))
  - FastJetUtil/FastJetProcessor: Print error message if value of Inclusive minPt parameter is not specified in the steering file.
    - set default value minPt="0.0"

* 2017-05-08 Lars Rickard Strom ([PR#8](https://github.com/iLCSoft/MarlinFastJet/pull/8))
  - FastJetUtil/FastJetProcessor: Added support for running ee_genkt_algorithm with non-default exponent and inclusive clustering, using parameters ee_genkt_algorithm R p, Inclusive minPt, reproduces the results of jet-based trimming

* 2017-05-05 Lars Rickard Strom ([PR#7](https://github.com/iLCSoft/MarlinFastJet/pull/7))
  - Updated instructions on how to use FastJetUtil.h in stand-alone project.
  - Moved getRecPars to FastJetUtil to enable broader use of conversion between PseudoJet objects and ReconstructedParticle (renamed convertFromPseudoJet).
  - Added FastJet JHTopTagger implementation and made changes in other files to accommodate such tool implementations.

# v00-04

* 2017-03-29 Andre Sailer ([PR#5](https://github.com/iLCSoft/MarlinFastJet/pull/5))
  - Update CI to not allow new warnings to be merged into MarlinFastJet

* 2017-03-24 Andre Sailer ([PR#6](https://github.com/iLCSoft/MarlinFastJet/pull/6))
  - Merged the FastJetClustering processor into the MarlinFastJet package
  - Fix Warnings in FastJetClustering processor

# v00-03

* 2017-03-29 Andre Sailer ([PR#5](https://github.com/iLCSoft/MarlinFastJet/pull/5))
  - Update CI to not allow new warnings to be merged into MarlinFastJet

* 2017-03-24 Andre Sailer ([PR#6](https://github.com/iLCSoft/MarlinFastJet/pull/6))
  - Merged the FastJetClustering processor into the MarlinFastJet package
  - Fix Warnings in FastJetClustering processor


