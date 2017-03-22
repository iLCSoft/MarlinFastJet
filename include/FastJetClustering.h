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

  FastJetClustering(const FastJetClustering&) = delete;
  FastJetClustering& operator=(const FastJetClustering&) = delete;
  
  virtual void init() ;

  virtual void processRunHeader( LCRunHeader* run ) ;
  
  virtual void processEvent( LCEvent * evt ) ; 
    
  virtual void check( LCEvent * evt ) ; 

  virtual void end() ;

  fastjet::JetAlgorithm  GetAlgorithm()  const {return fAlgorithm;}

  void SetAlgorithm(fastjet::JetAlgorithm f)  {fAlgorithm = f;}

 protected:

  // root file and tree objects
  TFile * _rootfile=NULL;
  TTree * _Etree=NULL;

  std::string _inputCollection{};
  std::string _outputCollection{};

  LCCollectionVec* _jetsCol=NULL;

  fastjet::JetAlgorithm fAlgorithm{};

  std::string sAlgorithm{};

  // px, py, pz, E, nptc
  float _jetVector[5][50];

  float _eCMS=0.0;

  int _nRun=0, _nEvt=0, _nJets=0, _nJetsHE=0;

  int _print=0, _nJetMax=0, _fillTree=0;

  double _RPar=0.0, _rp=0.0, _eJet=0.0;
} ;

#endif
