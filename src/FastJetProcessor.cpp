/*
 * FastJetProcessor.cpp
 *
 *  Created on: 25.05.2010
 *      Author: Lars Weuste (MPP Munich) - weuste@mpp.mpg.de
 *			iterative inclusive algorithm based on design by Marco Battaglia (CERN) - Marco.Battaglia@cern.ch
 */

#include "FastJetProcessor.h"
#include "FastJetUtil.h"

#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ParticleIDImpl.h>
#include <IMPL/VertexImpl.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>

#include <sstream>

FastJetProcessor aFastJetProcessor;

using namespace EVENT;

FastJetProcessor::FastJetProcessor() : Processor("FastJetProcessor"),
				       _lcParticleInName(""),
				       _lcParticleOutName(""),
				       _lcJetOutName(""),
				       _statsFoundJets(0),
				       _statsNrEvents(0),
				       _statsNrSkippedEmptyEvents(0),
				       _statsNrSkippedFixedNrJets(0),
				       _statsNrSkippedMaxIterations(0),
				       _storeParticlesInJets(false),
				       _fju(new FastJetUtil())
{
  _description = "Using the FastJet library to identify jets";

  // the input & output collections
  registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "recParticleIn", "a list of all reconstructed particles we are searching for jets in.", _lcParticleInName, "MCParticle");
  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE, "jetOut", "The identified jets", _lcJetOutName, "JetOut");

  registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE, "recParticleOut", "a list of all reconstructed particles used to make jets. If no value specified collection is not created", _lcParticleOutName, "");

  registerProcessorParameter(
			     "storeParticlesInJets",
			     "Store the list of particles that were clustered into jets in the recParticleOut collection",
			     _storeParticlesInJets,
			     false);


  _fju->registerFastJetParameters( this );


}

FastJetProcessor::FastJetProcessor(const FastJetProcessor& rhs):
  Processor( rhs.name() ),
  _lcParticleInName(rhs._lcParticleInName),
  _lcParticleOutName(rhs._lcParticleOutName),
  _lcJetOutName(rhs._lcJetOutName),
  _statsFoundJets(rhs._statsFoundJets),
  _statsNrEvents(rhs._statsNrEvents),
  _statsNrSkippedEmptyEvents(rhs._statsNrSkippedEmptyEvents),
  _statsNrSkippedFixedNrJets(rhs._statsNrSkippedFixedNrJets),
  _statsNrSkippedMaxIterations(rhs._statsNrSkippedMaxIterations),
  _storeParticlesInJets(rhs._storeParticlesInJets),
  _fju( new FastJetUtil(*rhs._fju) )
 {}
//	FastJetProcessor& operator=(const FastJetProcessor&) {}


FastJetProcessor::~FastJetProcessor() {
  delete _fju;
}


/** Called at the begin of the job before anything is read.
 * Use to initialize the processor, e.g. book histograms.
 */
void FastJetProcessor::init()
{
  // its always a good idea to ..
  printParameters();

  // parse the given steering parameters
  _fju->init();
  streamlog_out(MESSAGE) << "Jet Algorithm: " << _fju->_jetAlgo->description() << std::endl << std::endl;

  _statsFoundJets = 0;
  _statsNrEvents = 0;
  _statsNrSkippedEmptyEvents = 0;
  _statsNrSkippedFixedNrJets = 0;
  _statsNrSkippedMaxIterations = 0;

}

/** Called for every event - the working horse.
 */
void FastJetProcessor::processEvent(LCEvent * evt)
{
  // create output collection and save every jet with its particles in it
	auto lccJetsOut = std::make_unique<IMPL::LCCollectionVec>(LCIO::RECONSTRUCTEDPARTICLE);
  lccJetsOut->setSubset(true);
    // create output collection and save every particle which contributes to a jet
  auto lccParticlesOut = std::make_unique<IMPL::LCCollectionVec>(LCIO::RECONSTRUCTEDPARTICLE);
  lccParticlesOut->setSubset(true);

  LCCollection* particleIn(NULL);
  try
    {
      // get the input collection if existent
      particleIn = evt->getCollection(_lcParticleInName);
      if (particleIn->getNumberOfElements() < 1) {
	_statsNrSkippedEmptyEvents++;
	throw DataNotAvailableException("Collection is there, but its empty!");
      }

    } catch (const DataNotAvailableException& e) {
    streamlog_out(WARNING) << e.what() << std::endl << "Skipping" << std::endl;


    evt->addCollection(lccJetsOut.release(), _lcJetOutName);
    if (_storeParticlesInJets) evt->addCollection(lccParticlesOut.release(), _lcParticleOutName);

    return ;
  }

  // convert to pseudojet list
  PseudoJetList pjList = _fju->convertFromRecParticle(particleIn);

  PseudoJetList jets;
  fastjet::ClusterSequence cs = fastjet::ClusterSequence(pjList, *_fju->_jetAlgo);
  try {
    jets = _fju->clusterJets(pjList, cs, particleIn);
  } catch(const SkippedFixedNrJetException& e ) {
    _statsNrSkippedFixedNrJets++;
  } catch(const SkippedMaxIterationException& e ) {
    jets = e._jets;
    _statsNrSkippedMaxIterations++;
  }


  _statsNrEvents++;
  _statsFoundJets += jets.size();
  const unsigned nrJets = jets.size();

  PseudoJetList::iterator it;

  for (it=jets.begin(); it != jets.end(); it++) {
    // create a reconstructed particle for this jet, and add all the containing particles to it
    ReconstructedParticle* rec = _fju->convertFromPseudoJet( (*it), cs.constituents(*it), particleIn);
    lccJetsOut->addElement( rec );

    if (_storeParticlesInJets) {
      for (unsigned int n = 0; n < cs.constituents(*it).size(); ++n) {
	ReconstructedParticle* p =
	  static_cast<ReconstructedParticle*>(particleIn->getElementAt((cs.constituents(*it))[n].user_index()));
	lccParticlesOut->addElement( p );
      }
    }
  }

  // special case for the exclusive jet mode: we can save the transition y_cut value
  if (_fju->_clusterMode == FJ_exclusive_nJets && jets.size() == _fju->_requestedNumberOfJets) {
    // save the dcut value for this algorithm (although it might not be meaningful)
    LCParametersImpl &lccJetParams((LCParametersImpl &)lccJetsOut->parameters());

    lccJetParams.setValue(std::string("d_{n-1,n}"), (float)cs.exclusive_dmerge(nrJets-1));
    lccJetParams.setValue(std::string("d_{n,n+1}"), (float)cs.exclusive_dmerge(nrJets));
    lccJetParams.setValue(std::string("y_{n-1,n}"), (float)cs.exclusive_ymerge(nrJets-1));
    lccJetParams.setValue(std::string("y_{n,n+1}"), (float)cs.exclusive_ymerge(nrJets));
  }

  evt->addCollection(lccJetsOut.release(), _lcJetOutName);
  if (_storeParticlesInJets) evt->addCollection(lccParticlesOut.release(), _lcParticleOutName);

}

/** Called after data processing for clean up.
 */
void FastJetProcessor::end()
{
  streamlog_out(MESSAGE)
    << "Found jets: " << _statsFoundJets
    << " (" << (double)_statsFoundJets/_statsNrEvents << " per event) "
    << " - Skipped Empty events:" << _statsNrSkippedEmptyEvents
    << " - Skipped Events after max nr of iterations reached: " << _statsNrSkippedMaxIterations
    << " - Skipped Search for Fixed Nr Jets (due to insufficient nr of particles):" << _statsNrSkippedFixedNrJets
    << std::endl;

}
