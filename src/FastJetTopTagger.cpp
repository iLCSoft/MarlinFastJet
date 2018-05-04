/*
 * FastJetTopTagger.cpp
 *
 *  Created on: 26.09.2016 by Rickard Stroem (RS) (CERN) - lars.rickard.stroem@cern.ch
 *  Modified on 04.05.2018 by RS - adding substructure variables
 *      Marlin implementation of the JHTagger from FastJet
 *      incl. support for all jet algortihms, etc. implemented in fastjet 
 *      through new fastjet util header file
 */

#include <FastJetTopTagger.h>
#include <FastJetUtil.h>

#include <fastjet/Selector.hh>
#include "fastjet/contrib/Nsubjettiness.hh"

FastJetTopTagger aFastJetTopTagger;

using namespace EVENT;
using namespace IMPL;
using namespace fastjet;
using namespace fastjet::contrib;

FastJetTopTagger::FastJetTopTagger() : Processor("FastJetTopTagger"),
				       _lcParticleInName(""),
				       _lcParticleOutName(""),
				       _lcJetOutName(""),
				       _lcTopTaggerOutName(""),
				       _statsFoundJets(0),
				       _statsNrEvents(0),
				       _statsNrSkippedEmptyEvents(0),
				       _statsNrSkippedFixedNrJets(0),
				       _statsNrSkippedMaxIterations(0),
				       _storeParticlesInJets( false ),
				       _fju(new FastJetUtil()),
				       _beta(0.),
				       _energyCorrelator(""),
				       _axesMode(""),
				       _measureMode(""),
				       _deltaP(0.),
				       _deltaR(0.),
				       _cos_theta_W_max(0.),
				       _jhtoptagger(JHTopTagger()),
				       _energyCorrMap(),
				       _axesModeMap(),
				       _measureModeMap()

{
  _description = "Using the FastJet tool JHTagger to identify top jets";
  
  // the input & output collections
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "recParticleIn", 
			  "a list of all reconstructed particles we are searching for jets in.",
			  _lcParticleInName, "MCParticle");
  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE, "jetOut", 
			   "The identified jets", 
			   _lcJetOutName, "JetOut");
  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE, "recParticleOut", 
			   "a list of all reconstructed particles used to make jets. If no value specified collection is not created", 
			   _lcParticleOutName, "");
  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE, "topTaggerOut", 
			   "The top tagger output for each jet", _lcTopTaggerOutName, "TopTaggerOut");
  
  registerProcessorParameter("storeParticlesInJets",
			     "Store the list of particles that were clustered into jets in the recParticleOut collection",
			     _storeParticlesInJets,
			     false);
  
  //Fastjet parameters
  _fju->registerFastJetParameters( this );
  
  //Substructure parameters
  registerProcessorParameter("beta",
			     "Beta is called angular exponent and weights the angular distances between the jet constituents compated to their pt in the calculation of the energy correlation function and subjettiness.",
			     _beta,
			     1.);
  registerProcessorParameter("energyCorrelator",
			     "Options for contrib::EnergyCorrelator: pt_R (transverse momenta and boost-invariant angles), E_theta (energy and angle as dot product of vectors), E_inv (energy and angle as (2p_i*p_j/E_i E_j) (default: E_theta).",
			     _energyCorrelator,
			     std::string("E_theta"));
  registerProcessorParameter("axesMode",
			     "Option for NSubjettiness calculation. See fastjet NSubjettiness module for all options.",
			     _axesMode,
			     std::string("OnePass_WTA_KT_Axes"));
  registerProcessorParameter("measureMode",
			     "Option for NSubjettiness calculation. See fastjet NSubjettiness module for all options. Normalised measure recommended only for advanced users.",
			     _measureMode,
			     std::string("UnnormalizedMeasure"));
  
  //Top Tagger parameter
  registerProcessorParameter("deltaP",
			     "Subjets must carry at least this fraction of the original jet's p_t",
			     _deltaP,
			     0.05);
  registerProcessorParameter("deltaR",
			     "Subjets must be separated by at least this Manhattan distance",
			     _deltaR,
			     0.05);
  registerProcessorParameter("cos_theta_W_max",
			     "The maximal allowed value of the W helicity angle",
			     _cos_theta_W_max,
			     1.0);
}

FastJetTopTagger::~FastJetTopTagger(){
  delete _fju;
}

/** Called at the begin of the job before anything is read.
 * Use to initialize the processor, e.g. book histograms.
 */
void FastJetTopTagger::init()
{
  // its always a good idea to ..
  printParameters();
 
  // parse the given steering parameters
  _fju->init();
  streamlog_out(MESSAGE) << "Jet Algorithm: " << _fju->_jetAlgo->description() << std::endl << std::endl;
  
  // initate the top tagger
  _jhtoptagger = JHTopTagger(_deltaP, _deltaR, _cos_theta_W_max); //mW=80.4
  streamlog_out(MESSAGE) << "Top tagger implementation: " << _jhtoptagger.description() << std::endl;
  //_jhtoptagger.set_top_selector(SelectorMassRange(145, 205)); //<--Not used here, can be set in analysis
  //_jhtoptagger.set_W_selector(SelectorMassRange(65, 95)); //<--Not used here, can be set in analysis

  // initate substructure variables
  _energyCorrMap.insert(std::make_pair("pt_R", EnergyCorrelator::pt_R));
  _energyCorrMap.insert(std::make_pair("E_theta", EnergyCorrelator::E_theta));
  _energyCorrMap.insert(std::make_pair("E_inv", EnergyCorrelator::E_inv));

  _axesModeMap.insert(std::make_pair("KT_Axes", KT_Axes()));
  _axesModeMap.insert(std::make_pair("WTA_KT_Axes", WTA_KT_Axes()));
  _axesModeMap.insert(std::make_pair("OnePass_KT_Axes", OnePass_KT_Axes()));
  _axesModeMap.insert(std::make_pair("OnePass_WTA_KT_Axes", OnePass_WTA_KT_Axes()));

  _measureModeMap.insert(std::make_pair("UnnormalizedMeasure", UnnormalizedMeasure(_beta)));

  // counters
  _statsFoundJets = 0;
  _statsNrEvents = 0;
  _statsNrSkippedEmptyEvents = 0;
  _statsNrSkippedFixedNrJets = 0;
  _statsNrSkippedMaxIterations = 0;

} // end init

/** Called for every event - the working horse.
 */
void FastJetTopTagger::processEvent(LCEvent * evt){
	
  LCCollection* particleIn(NULL);
  try {
    // get the input collection if existent
    particleIn = evt->getCollection(_lcParticleInName);
    if (particleIn->getNumberOfElements() < 1){
      _statsNrSkippedEmptyEvents++;
      throw DataNotAvailableException("Collection is there, but its empty!");
    }
    
  } catch (DataNotAvailableException& e) {
    streamlog_out(WARNING) << e.what() << std::endl << "Skipping" << std::endl;

    //create dummy empty collection only in case there are processor that need the presence of them in later stages

    //create output collection and save every jet with its particles in it 
    LCCollectionVec* lccJetsOut = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

    //create output collection and save every particle which contributes to a jet
    LCCollectionVec* lccParticlesOut(NULL);
    if (_storeParticlesInJets){
      lccParticlesOut= new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
      lccParticlesOut->setSubset(true);
    }
    
    evt->addCollection(lccJetsOut, _lcJetOutName);
    if (_storeParticlesInJets) evt->addCollection(lccParticlesOut, _lcParticleOutName);

    //create output collection for the top jets
    LCCollectionVec* lccTopTaggerOut = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    LCCollectionVec* lccTopTaggerWOut = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    LCCollectionVec* lccTopTaggerW1Out = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    LCCollectionVec* lccTopTaggerW2Out = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    LCCollectionVec* lccTopTaggernonWOut = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    LCCollectionVec* lccTopTaggerCosThetaW = new LCCollectionVec(LCIO::LCGENERICOBJECT);
    
    evt->addCollection(lccTopTaggerOut, _lcTopTaggerOutName);
    evt->addCollection(lccTopTaggerWOut, _lcTopTaggerOutName+"_W");
    evt->addCollection(lccTopTaggernonWOut, _lcTopTaggerOutName+"_nonW");
    evt->addCollection(lccTopTaggerW1Out, _lcTopTaggerOutName+"_W1");
    evt->addCollection(lccTopTaggerW2Out, _lcTopTaggerOutName+"_W2");
    evt->addCollection(lccTopTaggerCosThetaW, _lcTopTaggerOutName+"_cos_theta_W");
    
    return;
  }
  
  // convert to pseudojet list
  PseudoJetList pjList = _fju->convertFromRecParticle(particleIn);
  
  //Jet finding
  PseudoJetList jets;
  try {
    // sort jets according to pt
    jets = sorted_by_pt(_fju->clusterJets(pjList, particleIn));
  } catch(SkippedFixedNrJetException& e){
    _statsNrSkippedFixedNrJets++; 
  } catch( SkippedMaxIterationException& e ) {
    jets = e._jets;
    _statsNrSkippedMaxIterations++; 
  }

  _statsNrEvents++;
  _statsFoundJets += jets.size();
  const unsigned nrJets = jets.size();
  
  // create output collection and save every jet with its particles in it
  LCCollectionVec* lccJetsOut = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  
  // create output collection and save every particle which contributes to a jet
  LCCollectionVec* lccParticlesOut(NULL);
  if (_storeParticlesInJets){
    lccParticlesOut= new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    lccParticlesOut->setSubset(true);
  }

  //Save TopTagger info of the jets into the lcio stream
  LCCollectionVec* lccTopTaggerOut = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec* lccTopTaggerWOut = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec* lccTopTaggerW1Out = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec* lccTopTaggerW2Out = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec* lccTopTaggernonWOut = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  
  //Save helicity information
  LCCollectionVec* lccTopTaggerCosThetaW = new LCCollectionVec(LCIO::LCGENERICOBJECT);
  LCGenericObjectImpl* lcgTopTaggerCosThetaW = new LCGenericObjectImpl(0, 0, 2);
  
  //Defining substructure functions
  Nsubjettiness nSubJettiness1(1, (_axesModeMap.find(_axesMode)->second).def(), (_measureModeMap.find(_measureMode)->second).def());
  Nsubjettiness nSubJettiness2(2, (_axesModeMap.find(_axesMode)->second).def(), (_measureModeMap.find(_measureMode)->second).def());
  Nsubjettiness nSubJettiness3(3, (_axesModeMap.find(_axesMode)->second).def(), (_measureModeMap.find(_measureMode)->second).def());
  
  //Save substructure functions
  LCCollectionVec* lccSubStructure = new LCCollectionVec(LCIO::LCGENERICOBJECT);
  LCGenericObjectImpl* lcgSubStructureC2 = new LCGenericObjectImpl(0, 0, 2);
  LCGenericObjectImpl* lcgSubStructureD2 = new LCGenericObjectImpl(0, 0, 2);
  LCGenericObjectImpl* lcgSubStructureC3 = new LCGenericObjectImpl(0, 0, 2);
  LCGenericObjectImpl* lcgSubStructureD3 = new LCGenericObjectImpl(0, 0, 2);
  LCGenericObjectImpl* lcgSubStructureTau1 = new LCGenericObjectImpl(0, 0, 2);
  LCGenericObjectImpl* lcgSubStructureTau2 = new LCGenericObjectImpl(0, 0, 2);
  LCGenericObjectImpl* lcgSubStructureTau3 = new LCGenericObjectImpl(0, 0, 2);
  
  //Loop over jets
  int index = 0;  
  PseudoJetList::iterator it;
  for(it=jets.begin(); it != jets.end(); it++, index++) {
    
    // create a reconstructed particle for this jet, and add all the containing particles to it
    ReconstructedParticle* rec = _fju->convertFromPseudoJet((*it), _fju->_cs->constituents(*it), particleIn);
    lccJetsOut->addElement( rec );
    
    if (_storeParticlesInJets) {
      for (unsigned int n = 0; n < _fju->_cs->constituents(*it).size(); ++n){
	ReconstructedParticle* p = static_cast<ReconstructedParticle*>(particleIn->getElementAt((_fju->_cs->constituents(*it))[n].user_index()));
        lccParticlesOut->addElement(p); 
      }
    }

    //Save substructure information

    //Substructure - energy correlation 
    double ECF1 = getECF(*it, 1, _beta, _energyCorrelator);
    double ECF2 = getECF(*it, 2, _beta, _energyCorrelator);   
    double ECF3 = getECF(*it, 3, _beta, _energyCorrelator);   
    double ECF4 = getECF(*it, 4, _beta, _energyCorrelator);   
    double C2 = ECF3*pow(ECF1,1)/pow(ECF2, 2); lcgSubStructureC2->setDoubleVal(index, C2); 
    double D2 = ECF3*pow(ECF1,3)/pow(ECF2, 3); lcgSubStructureD2->setDoubleVal(index, D2);
    double C3 = ECF4*pow(ECF2,1)/pow(ECF3, 2); lcgSubStructureC3->setDoubleVal(index, C3);
    double D3 = ECF4*pow(ECF2,3)/pow(ECF3, 3); lcgSubStructureD3->setDoubleVal(index, D3);
    
    //Substructure - NSubjettiness
    double tau1 = nSubJettiness1(*it); lcgSubStructureTau1->setDoubleVal(index, tau1);
    double tau2 = nSubJettiness2(*it); lcgSubStructureTau2->setDoubleVal(index, tau2);
    double tau3 = nSubJettiness3(*it); lcgSubStructureTau3->setDoubleVal(index, tau3);

    //John-Hopkins top tagger
    
    // search for top quark like structure in jet
    PseudoJet top_candidate = _jhtoptagger(*it);
    
    if (top_candidate == 0){ 
      
      lccTopTaggerOut->addElement(new ReconstructedParticleImpl());
      lccTopTaggerWOut->addElement(new ReconstructedParticleImpl());
      lccTopTaggernonWOut->addElement(new ReconstructedParticleImpl());
      lccTopTaggerW1Out->addElement(new ReconstructedParticleImpl());
      lccTopTaggerW2Out->addElement(new ReconstructedParticleImpl());
      lcgTopTaggerCosThetaW->setDoubleVal(index, 0.);

    } else {      
      // save top candidate
      ReconstructedParticle* t = _fju->convertFromPseudoJet(top_candidate, top_candidate.constituents(), particleIn);	    
      lccTopTaggerOut->addElement(t);

      // save W candidate
      PseudoJet top_candidate_W = top_candidate.structure_of<JHTopTagger>().W();
      ReconstructedParticle* W = _fju->convertFromPseudoJet(top_candidate_W, top_candidate_W.constituents(), particleIn);	    
      lccTopTaggerWOut->addElement(W);
      
      // save part 1 of W candidate
      PseudoJet top_candidate_W1 = top_candidate.structure_of<JHTopTagger>().W1();
      ReconstructedParticle* W1 = _fju->convertFromPseudoJet(top_candidate_W1, top_candidate_W1.constituents(), particleIn);	    
      lccTopTaggerW1Out->addElement(W1);
      
      // save part 2 of W candidate
      PseudoJet top_candidate_W2 = top_candidate.structure_of<JHTopTagger>().W2();
      ReconstructedParticle* W2 = _fju->convertFromPseudoJet(top_candidate_W2, top_candidate_W2.constituents(), particleIn);	    
      lccTopTaggerW2Out->addElement(W2);    

      // save non-W subjet of top candidate
      PseudoJet top_candidate_nonW = top_candidate.structure_of<JHTopTagger>().non_W();
      ReconstructedParticle* nonW = _fju->convertFromPseudoJet(top_candidate_nonW, top_candidate_nonW.constituents(), particleIn);	    
      lccTopTaggernonWOut->addElement(nonW);
         
      // save the polarisation angle of W
      double top_candidate_cos_theta_W = top_candidate.structure_of<JHTopTagger>().cos_theta_W();
      lcgTopTaggerCosThetaW->setDoubleVal(index, top_candidate_cos_theta_W);

    }
  }

  evt->addCollection(lccJetsOut, _lcJetOutName);
  if (_storeParticlesInJets) evt->addCollection(lccParticlesOut, _lcParticleOutName);
  
  lccSubStructure->addElement(lcgSubStructureC2);
  lccSubStructure->addElement(lcgSubStructureD2);
  lccSubStructure->addElement(lcgSubStructureC3);
  lccSubStructure->addElement(lcgSubStructureD3);
  lccSubStructure->addElement(lcgSubStructureTau1);
  lccSubStructure->addElement(lcgSubStructureTau2);
  lccSubStructure->addElement(lcgSubStructureTau3);
  evt->addCollection(lccSubStructure, "TopTaggerSubStructure");

  evt->addCollection(lccTopTaggerOut, _lcTopTaggerOutName);
  evt->addCollection(lccTopTaggerWOut, _lcTopTaggerOutName+"_W");
  evt->addCollection(lccTopTaggernonWOut, _lcTopTaggerOutName+"_nonW");
  evt->addCollection(lccTopTaggerW1Out, _lcTopTaggerOutName+"_W1");
  evt->addCollection(lccTopTaggerW2Out, _lcTopTaggerOutName+"_W2");

  lccTopTaggerCosThetaW->addElement(lcgTopTaggerCosThetaW);  
  evt->addCollection(lccTopTaggerCosThetaW, _lcTopTaggerOutName+"_cos_theta_W");
  
  // special case for the exclusive jet mode: we can save the transition y_cut value
  if (_fju->_clusterMode == FJ_exclusive_nJets && jets.size() == _fju->_requestedNumberOfJets) {
    // save the dcut value for this algorithm (although it might not be meaningful)
    LCParametersImpl &lccJetParams((LCParametersImpl &)lccJetsOut->parameters());
    
    lccJetParams.setValue(std::string("d_{n-1,n}"), (float)_fju->_cs->exclusive_dmerge(nrJets-1));
    lccJetParams.setValue(std::string("d_{n,n+1}"), (float)_fju->_cs->exclusive_dmerge(nrJets));
    lccJetParams.setValue(std::string("y_{n-1,n}"), (float)_fju->_cs->exclusive_ymerge(nrJets-1));
    lccJetParams.setValue(std::string("y_{n,n+1}"), (float)_fju->_cs->exclusive_ymerge(nrJets));
  }
  
} //end processEvent

double FastJetTopTagger::getECF(PseudoJet& jet, int whichECF, double beta, std::string energyCorrelator){
  /// ********** Energy Correlation Function *** /////////

  // options for EnergyCorrelator:
  //   pt_R (transverse momenta and boost-invariant angles)
  //   E_theta (energy and angle as dot product of vectors)
  //   E_inv (energy and angle as (2p_i \cdot p_j/E_i E_j)

  if(!jet.has_constituents())
    return -1.0;
    
  EnergyCorrelator ECF(whichECF, beta, _energyCorrMap[energyCorrelator]);

  return ECF(jet);
} //end

void FastJetTopTagger::end()
{
  streamlog_out(MESSAGE)
    << "Found jets: " << _statsFoundJets
    << " (" << (double)_statsFoundJets/_statsNrEvents << " per event) "
    << " - Skipped Empty events:" << _statsNrSkippedEmptyEvents
    << " - Skipped Events after max nr of iterations reached: " << _statsNrSkippedMaxIterations
    << " - Skipped Search for Fixed Nr Jets (due to insufficient nr of particles):" << _statsNrSkippedFixedNrJets
    << std::endl;
} //end end

std::ostream& operator<<(std::ostream& ostr, const PseudoJet& jet){
  ostr << "pt, y, phi =" << std::setprecision(6)
       << " " << std::setw(9) << jet.perp()
       << " " << std::setw(9) << jet.rap()
       << " " << std::setw(9) << jet.phi()
       << ", mass = " << std::setw(9) << jet.m();
  return ostr;
} // end ostream def

