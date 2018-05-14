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
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCGenericObjectImpl.h>

#include <fastjet/PseudoJet.hh>
#include <fastjet/SharedPtr.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/JHTopTagger.hh>
#include <fastjet/contrib/Njettiness.hh>
#include <fastjet/contrib/Nsubjettiness.hh>
#include <fastjet/contrib/EnergyCorrelator.hh>

//Forward declarations
class FastJetUtil;

typedef std::vector<fastjet::PseudoJet> PseudoJetList;

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
  std::string _lcSubStructureOutName;

  int _statsFoundJets;
  int _statsNrEvents;
  int _statsNrSkippedEmptyEvents;
  int _statsNrSkippedFixedNrJets;
  int _statsNrSkippedMaxIterations;
  bool _storeParticlesInJets;  

  FastJetUtil* _fju;

  bool _doSubstructure;
  double _beta;
  std::string _energyCorrelator;
  std::string _axesMode;
  std::string _measureMode;

  double _deltaP;
  double _deltaR;
  double _cos_theta_W_max;

  fastjet::JHTopTagger _jhtoptagger;

  FastJetTopTagger(const FastJetTopTagger& rhs) = delete;
  FastJetTopTagger & operator = (const FastJetTopTagger&) = delete;
  
  double getECF(fastjet::PseudoJet& jet, int whichECF, const std::string& energyCorr);

  // Simple class to store Axes along with a name
  class AxesStruct {
  private:
    // Shared Ptr so it handles memory management
    fastjet::SharedPtr<fastjet::contrib::AxesDefinition> _axes_def;
  public:
  AxesStruct(const fastjet::contrib::AxesDefinition & axes_def)
    : _axes_def(axes_def.create()) {}
    
    // Need special copy constructor to make it possible to put in a std::vector
  AxesStruct(const AxesStruct& myStruct)
    : _axes_def(myStruct._axes_def->create()) {}
    
    const fastjet::contrib::AxesDefinition & def() const {return *_axes_def;}
    std::string description() const {return _axes_def->description();}
    std::string short_description() const {return _axes_def->short_description();}
  };
  
  // Simple class to store Measures to make it easier to put in std::vector
  class MeasureStruct { 
  private:
    // Shared Ptr so it handles memory management
    fastjet::SharedPtr<fastjet::contrib::MeasureDefinition> _measure_def;
  public:
  MeasureStruct(const fastjet::contrib::MeasureDefinition& measure_def)
    : _measure_def(measure_def.create()) {}
    
    // Need special copy constructor to make it possible to put in a std::vector
    MeasureStruct(const MeasureStruct& myStruct)
    : _measure_def(myStruct._measure_def->create()) {}
    
    const fastjet::contrib::MeasureDefinition & def() const {return *_measure_def;}
    std::string description() const {return _measure_def->description();}  
  };

  std::map<std::string, fastjet::contrib::EnergyCorrelator::Measure> _energyCorrMeasureMap;
  int _maxECF;
  std::map<std::string, std::vector<fastjet::contrib::EnergyCorrelator> > _energyCorrMap;
  std::map<std::string, AxesStruct > _axesModeMap;
  std::map<std::string, MeasureStruct > _measureModeMap;
  std::map<int, fastjet::contrib::Nsubjettiness> _nsubjettinessMap;
    
};

std::ostream& operator<<(std::ostream&, const fastjet::PseudoJet&);

#endif /* FASTJETTOPTAGGER_H_ */

