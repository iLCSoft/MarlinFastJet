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


