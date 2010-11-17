#ifndef FastJetClustering_h
#define FastJetClustering_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include "UTIL/LCRelationNavigator.h"
#include <string>
#include <TFile.h>
#include <TTree.h>
#include "fastjet/JetDefinition.hh"

using namespace lcio ;
using namespace marlin ;

class FastJetClustering : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new FastJetClustering ; }
    
  FastJetClustering() ;
  
  virtual void init() ;

  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
    
  virtual void check( LCEvent * evt ) ; 

  virtual void end() ;

  fastjet::JetAlgorithm  GetAlgorithm()  const {return fAlgorithm;}

  void SetAlgorithm(fastjet::JetAlgorithm f)  {fAlgorithm = f;}

 protected:

  // root file and tree objects
  TFile * _rootfile;
  TTree * _Etree;

  std::string _inputCollection;
  std::string _outputCollection;

  LCCollectionVec* _jetsCol;

  fastjet::JetAlgorithm fAlgorithm; 

  std::string sAlgorithm; 

  // px, py, pz, E, nptc
  float _jetVector[5][50];

  float _eCMS;

  int _nRun, _nEvt, _nJets, _nJetsHE;

  int _print, _nJetMax, _fillTree;

  double _RPar, _rp, _eJet;
} ;

#endif
