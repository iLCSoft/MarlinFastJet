/*
 * FastJetTopTagger.cpp
 *
 *  Created on: 26.09.2016
 *      Author: Rickard Stroem (CERN) - lars.rickard.stroem@cern.ch
 *          Marlin implementation of the JHTagger from FastJet
 *          incl. support for all jet algortihms, etc. implemented in fastjet 
 *          through new fastjet util header file
 */

#include <FastJetTopTagger.h>
#include <FastJetUtil.h>

#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCGenericObjectImpl.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>

#include <fastjet/Selector.hh>
#include <fastjet/tools/JHTopTagger.hh>

#include <sstream>
#include <iomanip>

FastJetTopTagger aFastJetTopTagger;

using namespace EVENT;

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
				       _R(""),
				       _deltaP(0.),
				       _deltaR(0.),
				       _cos_theta_W_max(0.),
				       _jhtoptagger(fastjet::JHTopTagger())
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
  _jhtoptagger = fastjet::JHTopTagger(_deltaP, _deltaR, _cos_theta_W_max); //mW=80.4
  streamlog_out(MESSAGE) << "Top tagger implementation: " << _jhtoptagger.description() << std::endl;
  //_jhtoptagger.set_top_selector(fastjet::SelectorMassRange(145, 205));
  //_jhtoptagger.set_W_selector(fastjet::SelectorMassRange(65, 95));
  
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
    IMPL::LCCollectionVec* lccJetsOut = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

    //create output collection and save every particle which contributes to a jet
    IMPL::LCCollectionVec* lccParticlesOut(NULL);
    if (_storeParticlesInJets){
      lccParticlesOut= new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
      lccParticlesOut->setSubset(true);
    }
    
    evt->addCollection(lccJetsOut, _lcJetOutName);
    if (_storeParticlesInJets) evt->addCollection(lccParticlesOut, _lcParticleOutName);

    //create output collection for the top jets
    IMPL::LCCollectionVec* lccTopTaggerOut = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    IMPL::LCCollectionVec* lccTopTaggerWOut = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    IMPL::LCCollectionVec* lccTopTaggerW1Out = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    IMPL::LCCollectionVec* lccTopTaggerW2Out = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    IMPL::LCCollectionVec* lccTopTaggernonWOut = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    IMPL::LCCollectionVec* lccTopTaggerCosThetaW = new IMPL::LCCollectionVec(LCIO::LCGENERICOBJECT);
    
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
  IMPL::LCCollectionVec* lccJetsOut = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  
  // create output collection and save every particle which contributes to a jet
  IMPL::LCCollectionVec* lccParticlesOut(NULL);
  if (_storeParticlesInJets){
    lccParticlesOut= new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    lccParticlesOut->setSubset(true);
  }

  //Save TopTagger info of the jets into the lcio stream
  IMPL::LCCollectionVec* lccTopTaggerOut = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  IMPL::LCCollectionVec* lccTopTaggerWOut = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  IMPL::LCCollectionVec* lccTopTaggerW1Out = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  IMPL::LCCollectionVec* lccTopTaggerW2Out = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  IMPL::LCCollectionVec* lccTopTaggernonWOut = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  
  //Save helicity information
  IMPL::LCCollectionVec* lccTopTaggerCosThetaW = new IMPL::LCCollectionVec(LCIO::LCGENERICOBJECT);
  IMPL::LCGenericObjectImpl* lcgTopTaggerCosThetaW = new IMPL::LCGenericObjectImpl(0, 0, 2);
  
  int nTops = 0;  
  PseudoJetList::iterator it;
  for(it=jets.begin(); it != jets.end(); it++, nTops++) {
    
    // create a reconstructed particle for this jet, and add all the containing particles to it
    ReconstructedParticle* rec = _fju->convertFromPseudoJet((*it), _fju->_cs->constituents(*it), particleIn);
    lccJetsOut->addElement( rec );
    
    if (_storeParticlesInJets) {
      for (unsigned int n = 0; n < _fju->_cs->constituents(*it).size(); ++n){
	ReconstructedParticle* p = static_cast<ReconstructedParticle*>(particleIn->getElementAt((_fju->_cs->constituents(*it))[n].user_index()));
        lccParticlesOut->addElement(p); 
      }
    }
    
    // search for top quark like structure in jet
    fastjet::PseudoJet top_candidate = _jhtoptagger(*it);
    
    if (top_candidate == 0){ 
      
      lccTopTaggerOut->addElement(new ReconstructedParticleImpl());
      lccTopTaggerWOut->addElement(new ReconstructedParticleImpl());
      lccTopTaggernonWOut->addElement(new ReconstructedParticleImpl());
      lccTopTaggerW1Out->addElement(new ReconstructedParticleImpl());
      lccTopTaggerW2Out->addElement(new ReconstructedParticleImpl());
      lcgTopTaggerCosThetaW->setDoubleVal(nTops, 0.);

    } else {      
      // save top candidate
      ReconstructedParticle* t = _fju->convertFromPseudoJet(top_candidate, top_candidate.constituents(), particleIn);	    
      lccTopTaggerOut->addElement(t);
      
      // save W candidate
      fastjet::PseudoJet top_candidate_W = top_candidate.structure_of<fastjet::JHTopTagger>().W();
      ReconstructedParticle* W = _fju->convertFromPseudoJet(top_candidate_W, top_candidate_W.constituents(), particleIn);	    
      lccTopTaggerWOut->addElement(W);
      
      // save part 1 of W candidate
      fastjet::PseudoJet top_candidate_W1 = top_candidate.structure_of<fastjet::JHTopTagger>().W1();
      ReconstructedParticle* W1 = _fju->convertFromPseudoJet(top_candidate_W1, top_candidate_W1.constituents(), particleIn);	    
      lccTopTaggerW1Out->addElement(W1);
      
      // save part 2 of W candidate
      fastjet::PseudoJet top_candidate_W2 = top_candidate.structure_of<fastjet::JHTopTagger>().W2();
      ReconstructedParticle* W2 = _fju->convertFromPseudoJet(top_candidate_W2, top_candidate_W2.constituents(), particleIn);	    
      lccTopTaggerW2Out->addElement(W2);    

      // save non-W subjet of top candidate
      fastjet::PseudoJet top_candidate_nonW = top_candidate.structure_of<fastjet::JHTopTagger>().non_W();
      ReconstructedParticle* nonW = _fju->convertFromPseudoJet(top_candidate_nonW, top_candidate_nonW.constituents(), particleIn);	    
      lccTopTaggernonWOut->addElement(nonW);
         
      // save the polarisation angle of W
      double top_candidate_cos_theta_W = top_candidate.structure_of<fastjet::JHTopTagger>().cos_theta_W();
      lcgTopTaggerCosThetaW->setDoubleVal(nTops, top_candidate_cos_theta_W);

    }
  }

  evt->addCollection(lccJetsOut, _lcJetOutName);
  if (_storeParticlesInJets) evt->addCollection(lccParticlesOut, _lcParticleOutName);
  
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

std::ostream& operator<<(std::ostream& ostr, const fastjet::PseudoJet& jet){
  ostr << "pt, y, phi =" << std::setprecision(6)
       << " " << std::setw(9) << jet.perp()
       << " " << std::setw(9) << jet.rap()
       << " " << std::setw(9) << jet.phi()
       << ", mass = " << std::setw(9) << jet.m();
  return ostr;
} // end ostream def

