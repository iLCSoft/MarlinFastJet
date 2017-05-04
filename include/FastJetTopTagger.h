/*
 * FastJetTopTagger.h
 *
 *  Created on: 26.09.2016
 *      Author: Rickard Stroem (CERN) - lars.rickard.stroem@cern.ch
 *          Marlin implementation of the JHTagger from FastJet
 *          incl. support for all jet algortihms, etc. implemented in fastjet 
 *          through new fastjet util header file
 */

#ifndef FASTJETTOPTAGGER_H_
#define FASTJETTOPTAGGER_H_

#include <marlin/Processor.h>
#include <marlin/VerbosityLevels.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include "IMPL/LCCollectionVec.h"
#include <LCIOSTLTypes.h>

#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/tools/JHTopTagger.hh>

//Forward declaration
class FastJetUtil;
typedef std::vector< fastjet::PseudoJet > PseudoJetList;

class FastJetTopTagger : marlin::Processor {  
 public:
  FastJetTopTagger();
  virtual ~FastJetTopTagger();
  
  virtual Processor* newProcessor(){ 
    return new FastJetTopTagger(); 
  }
  
  /** Called at the begin of the job before anything is read.
      Use to initialize the processor, e.g. book histograms.
  */
  virtual void init();
  
  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader*){};
  
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
  std::string _lcParticleInName;
  std::string _lcParticleOutName;
  std::string _lcJetOutName;
  std::string _lcTopTaggerOutName;

  int _statsFoundJets;
  int _statsNrEvents;
  int _statsNrSkippedEmptyEvents;
  int _statsNrSkippedFixedNrJets;
  int _statsNrSkippedMaxIterations;
  bool _storeParticlesInJets;  

  FastJetUtil* _fju;

  std::string _R;
  double _deltaP;
  double _deltaR;
  double _cos_theta_W_max;

  // list of pseudo particles / jets
  fastjet::JHTopTagger _jhtoptagger;

 private:
  FastJetTopTagger(const FastJetTopTagger& rhs) = delete;
  FastJetTopTagger & operator = (const FastJetTopTagger&) = delete;

};

std::ostream& operator<<(std::ostream&, const fastjet::PseudoJet&);

#endif /* FASTJETTOPTAGGER_H_ */
