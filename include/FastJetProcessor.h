/*
 * FastJetProcessor.h
 *
 *  Created on: 25.05.2010
 *      Author: Lars Weuste (MPP Munich) - weuste@mpp.mpg.de
 *			iterative inclusive algorithm based on design by Marco Battaglia (CERN) - Marco.Battaglia@cern.ch
 */

#ifndef FASTJETPROCESSOR_H_
#define FASTJETPROCESSOR_H_

#include "EClusterMode.h"

#include "marlin/Processor.h"
#include "marlin/VerbosityLevels.h"
#include "EVENT/LCCollection.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/LCCollectionVec.h"
#include "LCIOSTLTypes.h"

#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>

#include <vector>
#include <string>

using namespace marlin;
using namespace lcio;
using namespace std;
using namespace fastjet;

//Forward Declaration
class FastJetUtil;

typedef vector< fastjet::PseudoJet > PseudoJetList;


class FastJetProcessor : Processor {
public:
  FastJetProcessor();
  virtual ~FastJetProcessor();


  virtual Processor* newProcessor() {
    return new FastJetProcessor();
  }

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader* ){};

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent * evt);

  virtual void check(LCEvent * ) {};

  /** Called after data processing for clean up.
   */
  virtual void end();

  friend class FastJetUtil;

private:

  // the LC Collection names for input/output
  string	_lcParticleInName;
  string	_lcParticleOutName;
  string	_lcJetOutName;

  EVENT::ReconstructedParticle* getRecPar(fastjet::PseudoJet& fj, const vector< fastjet::PseudoJet >& constituents);

  LCCollection*	_reconstructedPars;

  int	_statsFoundJets;
  int _statsNrEvents;
  int _statsNrSkippedEmptyEvents;
  int _statsNrSkippedFixedNrJets;
  int _statsNrSkippedMaxIterations;
  bool _storeParticlesInJets;

  FastJetUtil* _fju;

private:
  FastJetProcessor(const FastJetProcessor& rhs);
  FastJetProcessor& operator=(const FastJetProcessor&) { return *this; }

};

ostream& operator<<(ostream& out, EClusterMode& m);

#endif /* FASTJETPROCESSOR_H_ */
