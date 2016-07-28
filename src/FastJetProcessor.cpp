/*
 * FastJetProcessor.cpp
 *
 *  Created on: 25.05.2010
 *      Author: Lars Weuste (MPP Munich) - weuste@mpp.mpg.de
 *      		iterative inclusive algorithm based on design by Marco Battaglia (CERN) - Marco.Battaglia@cern.ch
 */

#include "FastJetProcessor.h"

#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ParticleIDImpl.h>
#include <IMPL/VertexImpl.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/MCParticle.h>


#include <fastjet/ClusterSequence.hh>
// plugins: (maybe some #ifdefs here?)

#include <fastjet/SISConePlugin.hh>
#include <fastjet/SISConeSphericalPlugin.hh>
#include <fastjet/CDFMidPointPlugin.hh>
#include <fastjet/CDFJetCluPlugin.hh>
#include <fastjet/NestedDefsPlugin.hh>
#include <fastjet/EECambridgePlugin.hh>
#include <fastjet/JadePlugin.hh>

#include <fastjet/contrib/ValenciaPlugin.hh>

//#include <fastjet/D0RunIIConePlugin.hh>
//#include <fastjet/TrackJetPlugin.hh>
//#include <fastjet/ATLASConePlugin.hh>
//#include <fastjet/CMSIterativeConePlugin.hh>


#include <sstream>

#define ITERATIVE_INCLUSIVE_MAX_ITERATIONS 20


FastJetProcessor aFastJetProcessor;

using namespace EVENT;

FastJetProcessor::FastJetProcessor() : Processor("FastJetProcessor")
{
	_description = "Using the FastJet library to identify jets";

	// the input & output collections
	registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE, "recParticleIn", "a list of all reconstructed particles we are searching for jets in.", _lcParticleInName, "MCParticle");
	registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE, "jetOut", "The identified jets", _lcJetOutName, "JetOut");

	registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE, "recParticleOut", "a list of all reconstructed particles used to make jets. If no value specified collection is not created", _lcParticleOutName, "");

	// the parameters. See description for details
	StringVec defAlgoAndParam;
	defAlgoAndParam.push_back("kt_algorithm");
	defAlgoAndParam.push_back("0.7");
	registerProcessorParameter(
			"algorithm",
			"Selects the algorithm and its parameters. E.g. 'kt_algorithm 0.7' or 'ee_kt_algorithm'. For a full list of supported algorithms, see the logfile after execution.",
			_jetAlgoNameAndParams,
			defAlgoAndParam);

	registerProcessorParameter(
			"recombinationScheme",
			"The recombination scheme used when merging 2 particles. Usually there is no need to use anything else than 4-Vector addition: E_scheme",
			_jetRecoSchemeName,
			string("E_scheme"));

	registerProcessorParameter(
			"storeParticlesInJets",
			"Store the list of particles that were clustered into jets in the recParticleOut collection",
			_storeParticlesInJets,
			false);

	StringVec defClusterMode;
	defClusterMode.push_back("Inclusive");
	registerProcessorParameter(
			"clusteringMode",
			"One of 'Inclusive <minPt>', 'InclusiveIterativeNJets <nrJets> <minE>', 'ExclusiveNJets <nrJets>', 'ExclusiveYCut <yCut>'. Note: not all modes are available for all algorithms.",
			_clusterModeNameAndParam,
			defClusterMode);


}

FastJetProcessor::~FastJetProcessor()
{
}


/** Called at the begin of the job before anything is read.
 * Use to initialize the processor, e.g. book histograms.
 */
void FastJetProcessor::init()
{
	// its always a good idea to ..
	printParameters();

	// parse the given steering parameters
	this->initStrategy();
	this->initRecoScheme();
	this->initClusterMode();
	this->initJetAlgo();

	_statsFoundJets = 0;
	_statsNrEvents = 0;
	_statsNrSkippedEmptyEvents = 0;
	_statsNrSkippedFixedNrJets = 0;
	_statsNrSkippedMaxIterations = 0;

}


void FastJetProcessor::initClusterMode()
{
	// parse the clustermode string

	// at least a name has to be given
	if (_clusterModeNameAndParam.size() == 0)
		throw Exception("Cluster mode not specified");

	// save the name of the cluster mode
	_clusterModeName = _clusterModeNameAndParam[0];

	_clusterMode = NONE;

	// check the different cluster mode possibilities, and check if the number of parameters are correct
	if (_clusterModeName.compare("Inclusive") == 0)
	{
		if (_clusterModeNameAndParam.size() != 2)
			throw Exception("Wrong Parameter(s) for Clustering Mode. Expected:\n <parameter name=\"clusteringMode\" type=\"StringVec\"> Inclusive <minPt> </parameter>");
		_minPt = atof(_clusterModeNameAndParam[1].c_str());
		_clusterMode = FJ_inclusive;
	}
	else if (_clusterModeName.compare("InclusiveIterativeNJets") == 0)
	{
		if (_clusterModeNameAndParam.size() != 3)
			throw Exception("Wrong Parameter(s) for Clustering Mode. Expected:\n <parameter name=\"clusteringMode\" type=\"StringVec\"> InclusiveIterativeNJets <NJets> <minE> </parameter>");
		_nrJets = atoi(_clusterModeNameAndParam[1].c_str());
		_minE = atoi(_clusterModeNameAndParam[2].c_str());
		_clusterMode = OWN_inclusiveIteration;
	}
	else if (_clusterModeName.compare("ExclusiveNJets") == 0)
	{
		if (_clusterModeNameAndParam.size() != 2)
			throw Exception("Wrong Parameter(s) for Clustering Mode. Expected:\n <parameter name=\"clusteringMode\" type=\"StringVec\"> ExclusiveNJets <NJets> </parameter>");
		_nrJets = atoi(_clusterModeNameAndParam[1].c_str());
		_clusterMode = FJ_exclusive_nJets;
	}
	else if (_clusterModeName.compare("ExclusiveYCut") == 0)
	{
		if (_clusterModeNameAndParam.size() != 2)
			throw Exception("Wrong Parameter(s) for Clustering Mode. Expected:\n <parameter name=\"clusteringMode\" type=\"StringVec\"> ExclusiveYCut <yCut> </parameter>");
		_yCut = atof(_clusterModeNameAndParam[1].c_str());
		_clusterMode = FJ_exclusive_yCut;
	}
	else
	{
		throw Exception("Unknown cluster mode.");
	}

	streamlog_out(MESSAGE) << "cluster mode: " << _clusterMode << endl;
}

ostream& operator<<(ostream& out, EClusterMode& m)
{
	switch (m)
	{
	case OWN_inclusiveIteration:	out << "InclusiveIterativeNJets"; break;
	case FJ_inclusive:				out << "Inclusive"; break;
	case FJ_exclusive_nJets:		out << "ExclusiveNJets"; break;
	case FJ_exclusive_yCut:			out << "ExclusiveYCut"; break;
	default:						out << "unknown"; break;
	}

	return out;
}

void FastJetProcessor::initJetAlgo()
{

	// parse the given jet algorithm string

	// sanity check
	if (_jetAlgoNameAndParams.size() < 1)
		throw Exception("No Jet algorithm provided!");

	// save the name
	this->_jetAlgoName = _jetAlgoNameAndParams[0];

	// check all supported algorithms and create the appropriate FJ instance

	_jetAlgo = NULL;
	streamlog_out(MESSAGE) << "Algorithms: ";	// the isJetAlgo function will write to streamlog_out(MESSAGE), so that we get a list of available algorithms in the log

	// example: kt_algorithm, needs 1 parameter, supports inclusive, inclusiveIterative, exlusiveNJets and exlusiveYCut clustering
	if (isJetAlgo("kt_algorithm", 1, FJ_inclusive | FJ_exclusive_nJets | FJ_exclusive_yCut | OWN_inclusiveIteration))
	{
		_jetAlgoType = kt_algorithm;
		_jetAlgo = new fastjet::JetDefinition(
				_jetAlgoType, atof(_jetAlgoNameAndParams[1].c_str()), _jetRecoScheme, _strategy);

	}

	if (isJetAlgo("cambridge_algorithm", 1, FJ_inclusive | FJ_exclusive_nJets | FJ_exclusive_yCut | OWN_inclusiveIteration))
	{
		_jetAlgoType = cambridge_algorithm;
		_jetAlgo = new fastjet::JetDefinition(
				_jetAlgoType, atof(_jetAlgoNameAndParams[1].c_str()), _jetRecoScheme, _strategy);
	}

	if (isJetAlgo("antikt_algorithm", 1, FJ_inclusive | OWN_inclusiveIteration))
	{
		_jetAlgoType = antikt_algorithm;
		_jetAlgo = new fastjet::JetDefinition(
				_jetAlgoType, atof(_jetAlgoNameAndParams[1].c_str()), _jetRecoScheme, _strategy);
	}

	if (isJetAlgo("genkt_algorithm", 2, FJ_inclusive | OWN_inclusiveIteration | FJ_exclusive_nJets | FJ_exclusive_yCut))
	{
		_jetAlgoType = genkt_algorithm;
		_jetAlgo = new fastjet::JetDefinition(
				_jetAlgoType, atof(_jetAlgoNameAndParams[1].c_str()), atof(_jetAlgoNameAndParams[2].c_str()), _jetRecoScheme, _strategy);
	}

	if (isJetAlgo("cambridge_for_passive_algorithm", 1, FJ_inclusive | OWN_inclusiveIteration | FJ_exclusive_nJets | FJ_exclusive_yCut))
	{
		_jetAlgoType = cambridge_for_passive_algorithm;
		_jetAlgo = new fastjet::JetDefinition(
				_jetAlgoType, atof(_jetAlgoNameAndParams[1].c_str()), _jetRecoScheme, _strategy);
	}

	if (isJetAlgo("genkt_for_passive_algorithm", 1, FJ_inclusive | OWN_inclusiveIteration))
	{
		_jetAlgoType = genkt_for_passive_algorithm;
		_jetAlgo = new fastjet::JetDefinition(
				_jetAlgoType, atof(_jetAlgoNameAndParams[1].c_str()), _jetRecoScheme, _strategy);
	}


	if (isJetAlgo("ee_kt_algorithm", 0, FJ_exclusive_nJets | FJ_exclusive_yCut))
	{
		_jetAlgoType = ee_kt_algorithm;
		_jetAlgo = new fastjet::JetDefinition(
				_jetAlgoType, _jetRecoScheme, _strategy);
	}

	if (isJetAlgo("ee_genkt_algorithm", 1, FJ_exclusive_nJets | FJ_exclusive_yCut))
	{
		_jetAlgoType = ee_genkt_algorithm;
		_jetAlgo = new fastjet::JetDefinition(
				_jetAlgoType, atof(_jetAlgoNameAndParams[1].c_str()), _jetRecoScheme, _strategy);
	}

	if (isJetAlgo("SISConePlugin", 2, FJ_inclusive | OWN_inclusiveIteration))
	{
		fastjet::SISConePlugin* pl;

		pl = new fastjet::SISConePlugin(
				atof(_jetAlgoNameAndParams[1].c_str()),
				atof(_jetAlgoNameAndParams[2].c_str())
				);


		_jetAlgo = new fastjet::JetDefinition(pl);
	}

	if (isJetAlgo("SISConeSphericalPlugin", 2, FJ_inclusive | OWN_inclusiveIteration))
	{

		fastjet::SISConeSphericalPlugin* pl;

		pl = new fastjet::SISConeSphericalPlugin(
				atof(_jetAlgoNameAndParams[1].c_str()),
				atof(_jetAlgoNameAndParams[2].c_str())
				);

		_jetAlgo = new fastjet::JetDefinition(pl);
	}

	if (isJetAlgo("ValenciaPlugin", 3, FJ_exclusive_nJets | FJ_exclusive_yCut))
	{

	  fastjet::contrib::ValenciaPlugin* pl;
	  pl = new fastjet::contrib::ValenciaPlugin(
						    atof(_jetAlgoNameAndParams[1].c_str()),  // R value
						    atof(_jetAlgoNameAndParams[2].c_str()),  // beta value 
						    atof(_jetAlgoNameAndParams[3].c_str())   // gamma value 
						    ); 

	  _jetAlgo = new fastjet::JetDefinition(pl);
	}

// TODO: Maybe we should complete the user defined (ATLAS, CMS, ..) algorithms
//	if (isJetAlgo("ATLASConePlugin", 1))
//	{
//		fastjet::ATLASConePlugin* pl;
//		pl = new fastjet::ATLASConePlugin(
//				atof(_jetAlgoAndParamsString[1].c_str())
//				);
//
//		_jetAlgo = new fastjet::JetDefinition(pl);
//	}
//
//	if (isJetAlgo("CMSIterativeConePlugin", 1))
//	{
//		fastjet::CMSIterativeConePlugin* pl;
//		pl = new fastjet::CMSIterativeConePlugin(
//				atof(_jetAlgoAndParamsString[1].c_str())
//				);
//
//		_jetAlgo = new fastjet::JetDefinition(pl);
//	}

	streamlog_out(MESSAGE) << endl;	// end of list of available algorithms

	if (!_jetAlgo)
	{
		streamlog_out(ERROR) << "The given algorithm \"" << _jetAlgoName
					<< "\" is unknown to me!" << endl;
		throw Exception("Unknown FastJet algorithm.");
	}

	streamlog_out(MESSAGE) << "jet algorithm: " << _jetAlgo->description() << endl;
}


bool FastJetProcessor::isJetAlgo(string algo, int nrParams, int supportedModes)
{
	streamlog_out(MESSAGE) << " " << algo;

	// check if the chosen algorithm is the same as it was passed
	if (_jetAlgoName.compare(algo) != 0)
		return false;

	streamlog_out(MESSAGE) << "*";	// mark the current algorithm as the selected one, even before we did our checks on nr of parameters etc.

	// check if we have enough number of parameters
	if ((int)_jetAlgoNameAndParams.size() - 1 != nrParams)
	{
		streamlog_out(ERROR) << endl
				<< "Wrong numbers of parameters for algorithm: " << algo << endl
				<< "We need " << nrParams << " params, but we got " << _jetAlgoNameAndParams.size() - 1 << endl;
		throw Exception("You have insufficient number of parameters for this algorithm! See log for more details.");
	}

	// check if the mode is supported via a binary AND
	if ((supportedModes & _clusterMode) != _clusterMode)
	{
		streamlog_out(ERROR) << endl
				<< "This algorithm is not capable of running in this clustering mode (" << _clusterMode << "). Sorry!" << endl;
		throw Exception("This algorithm is not capable of running in this mode");

	}

	return true;
}

void FastJetProcessor::initRecoScheme()
{
	// parse the steering parameter for the recombination scheme

	if (_jetRecoSchemeName.compare("E_scheme") == 0)
		_jetRecoScheme = fastjet::E_scheme;
	else if (_jetRecoSchemeName.compare("pt_scheme") == 0)
		_jetRecoScheme = fastjet::pt_scheme;
	else if (_jetRecoSchemeName.compare("pt2_scheme") == 0)
		_jetRecoScheme = fastjet::pt2_scheme;
	else if (_jetRecoSchemeName.compare("Et_scheme") == 0)
		_jetRecoScheme = fastjet::Et_scheme;
	else if (_jetRecoSchemeName.compare("Et2_scheme") == 0)
		_jetRecoScheme = fastjet::Et2_scheme;
	else if (_jetRecoSchemeName.compare("BIpt_scheme") == 0)
		_jetRecoScheme = fastjet::BIpt_scheme;
	else if (_jetRecoSchemeName.compare("BIpt2_scheme") == 0)
		_jetRecoScheme = fastjet::BIpt2_scheme;
	else
	{
		streamlog_out(ERROR)
				<< "Unknown recombination scheme: " << _jetRecoSchemeName << endl;
		throw Exception("Unknown FastJet recombination scheme! See log for more details.");
	}

	streamlog_out(MESSAGE) << "recombination scheme: " << _jetRecoSchemeName << endl;
}

void FastJetProcessor::initStrategy()
{
	// we could provide a steering parameter for this and it would be parsed here.
	// however, when using the 'Best' clustering strategy FJ will chose automatically from a big list
	// changing this would (most likely) only change the speed of calculation, not the outcome
	_strategy = fastjet::Best;
	_strategyName = "Best";
	streamlog_out(MESSAGE) << "Strategy: " << _strategyName << endl;
}



/** Called for every event - the working horse.
 */
void FastJetProcessor::processEvent(LCEvent * evt)
{

	// clear the list containing the (pseudo)-jets and the FJ-index-to-PFOReconstructedParticle map
	_pjList.clear();

	try
	{
		// get the input collection if existent
		LCCollection* particleIn = evt->getCollection(_lcParticleInName);

		if (particleIn->getNumberOfElements() < 1)
		{
			_statsNrSkippedEmptyEvents++;
			throw DataNotAvailableException("Collection is there, but its empty!");
		}

		// convert to pseudojet list

		_reconstructedPars = particleIn;
		convertFromRecParticle(particleIn);


	} catch (DataNotAvailableException e) {
		streamlog_out(WARNING) << e.what() << endl << "Skipping" << endl;

		//create dummy empty collection only in case there are processor that need the presence of them in later stages

		// create output collection and save every jet with its particles in it
		IMPL::LCCollectionVec* lccJetsOut = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
		// create output collection and save every particle which contributes to a jet

		IMPL::LCCollectionVec* lccParticlesOut(NULL);
		if (_storeParticlesInJets){
		  lccParticlesOut= new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
		  lccParticlesOut->setSubset(true);
		}

		evt->addCollection(lccJetsOut, _lcJetOutName);
		if (_storeParticlesInJets) evt->addCollection(lccParticlesOut, _lcParticleOutName);
		

		return ;
	}

	///////////////////////////////
	// do the jet finding for the user defined parameter jet finder
	ClusterSequence cs(_pjList, *_jetAlgo);

	vector< fastjet::PseudoJet > jets;		// will contain the found jets

	if (_clusterMode == FJ_inclusive)
	{
		jets = cs.inclusive_jets(_minPt);
	}
	else if (_clusterMode == FJ_exclusive_yCut)
	{
		jets = cs.exclusive_jets_ycut(_yCut);
	}
	else if (_clusterMode == FJ_exclusive_nJets)
	{
		// sanity check: if we have not enough particles, FJ will cause an assert
		if (_reconstructedPars->getNumberOfElements() < (int)_nrJets)
		{
			streamlog_out(WARNING) << "Not enough elements in the input collection to create " << _nrJets << " jets." << endl;
			_statsNrSkippedFixedNrJets++;
		}
		else
		{
			jets = cs.exclusive_jets((int)(_nrJets));
		}
	}
	else if (_clusterMode == OWN_inclusiveIteration)
	{
		// sanity check: if we have not enough particles, FJ will cause an assert
		if (_reconstructedPars->getNumberOfElements() < (int)_nrJets)
		{
			streamlog_out(WARNING) << "Not enough elements in the input collection to create " << _nrJets << " jets." << endl;
			_statsNrSkippedFixedNrJets++;
		}
		else
		{
			jets = this->doIterativeInclusiveClustering();
		}
	}


	//////////////////////////////
	// save the jets into the lcio stream

	_statsNrEvents++;
	_statsFoundJets += jets.size();
	unsigned nrJets = jets.size();

	// create output collection and save every jet with its particles in it
	IMPL::LCCollectionVec* lccJetsOut = new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	// create output collection and save every particle which contributes to a jet
	IMPL::LCCollectionVec* lccParticlesOut(NULL);
	if (_storeParticlesInJets){
	  lccParticlesOut= new IMPL::LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
	  lccParticlesOut->setSubset(true);
	}

	vector< fastjet::PseudoJet >::iterator it;

	for (it=jets.begin(); it != jets.end(); it++)
	{
		// create a reconstructed particle for this jet, and add all the containing particles to it
		ReconstructedParticle* rec = getRecPar( (*it), cs.constituents(*it) );
		lccJetsOut->addElement( rec );

		if (_storeParticlesInJets) 
		{
			for (unsigned int n = 0; n < cs.constituents(*it).size(); ++n)
			{
				ReconstructedParticle* p = dynamic_cast< ReconstructedParticle* > (_reconstructedPars->getElementAt((cs.constituents(*it))[n].user_index()));
				lccParticlesOut->addElement( p ); 
			}
		}
	}

	evt->addCollection(lccJetsOut, _lcJetOutName);
	if (_storeParticlesInJets) evt->addCollection(lccParticlesOut, _lcParticleOutName);

	// special case for the exclusive jet mode: we can save the transition y_cut value
	if (_clusterMode == FJ_exclusive_nJets && jets.size() == _nrJets)
	{
		// save the dcut value for this algorithm (although it might not be meaningful)
		LCParametersImpl &lccJetParams((LCParametersImpl &)lccJetsOut->parameters());

		lccJetParams.setValue(string("d_{n-1,n}"), (float)cs.exclusive_dmerge(nrJets-1));
		lccJetParams.setValue(string("d_{n,n+1}"), (float)cs.exclusive_dmerge(nrJets));

		lccJetParams.setValue(string("y_{n-1,n}"), (float)cs.exclusive_ymerge(nrJets-1));
		lccJetParams.setValue(string("y_{n,n+1}"), (float)cs.exclusive_ymerge(nrJets));
	}


}

vector < fastjet::PseudoJet > FastJetProcessor::doIterativeInclusiveClustering()
{
	// lets do a iterative procedure until we found the correct number of jets
	// for that we will do inclusive clustering, modifying the R parameter in some kind of minimization
	// this is based on Marco Battaglia's FastJetClustering

	double R = M_PI_4;	// maximum of R is Pi/2, minimum is 0. So we start hat Pi/4
	double RDiff = R / 2;	// the step size we modify the R parameter at each iteration. Its size for the n-th step is R/(2n), i.e. starts with R/2
	vector< fastjet::PseudoJet > jets;
	vector< fastjet::PseudoJet > jetsReturn;
	unsigned nJets;
	int iIter = 0;	// nr of current iteration

	// these variables are only used if the SisCone(Spherical)Plugin is selected
	// This is necessary, as these are plugins and hence use a different constructor than
	// the built in fastjet algorithms

	// here we save pointer to the plugins, so that if they are created, we can delete them again after usage
	fastjet::SISConePlugin* pluginSisCone = NULL;
	fastjet::SISConeSphericalPlugin* pluginSisConeSph = NULL;
	// check if we use the siscones
	bool useSisCone = _jetAlgoName.compare("SISConePlugin") == 0;
	bool useSisConeSph = _jetAlgoName.compare("SISConeSphericalPlugin") == 0;
	// save the 2nd parameter of the SisCones
	double sisConeOverlapThreshold = 0;
	if (useSisCone || useSisConeSph)
		sisConeOverlapThreshold = atof(_jetAlgoNameAndParams[2].c_str());

	// do a maximum of N iterations
	// For each iteration we will modify the Radius R by RDiff. RDiff will be reduced by 2 for each iteration
	// i.e. RDiff = Pi/8, Pi/16, Pi/32, Pi/64, ...
	// e.g.
	// R = Pi/4
	// R = Pi/4 + Pi/8
	// R = Pi/4 + Pi/8 - Pi/16
	// R = Pi/4 + Pi/8 - Pi/16 + Pi/32
	// ...
	for (iIter=0; iIter<ITERATIVE_INCLUSIVE_MAX_ITERATIONS; iIter++)
	{

		// do the clustering for this value of R. For this we need to re-initialize the JetDefinition, as it takes the R parameter
		JetDefinition* def = NULL;

		// unfortunately SisCone(spherical) are being initialized differently, so we have to check for this
		if (useSisCone)
		{
			pluginSisCone = new fastjet::SISConePlugin(R, sisConeOverlapThreshold);
			def = new fastjet::JetDefinition(pluginSisCone);
		}
		else if (useSisConeSph)
		{
			pluginSisConeSph = new fastjet::SISConeSphericalPlugin(R, sisConeOverlapThreshold);
			def = new fastjet::JetDefinition(pluginSisConeSph);
		}
		else
		{
			def = new JetDefinition(_jetAlgoType, R, _jetRecoScheme, _strategy);
		}

		// noe we can finally create the cluster sequence
		ClusterSequence cs(_pjList, *def);

		jets = cs.inclusive_jets(0);	// no pt cut, we will do an energy cut
		jetsReturn.clear();

		// count the number of jets above threshold
		nJets = 0;
		for (unsigned j=0; j<jets.size(); j++)
			if (jets[j].E() > _minE)
				jetsReturn.push_back(jets[j]);
		nJets = jetsReturn.size();

		streamlog_out(DEBUG) << iIter << " " << R << " " << jets.size() << " " << nJets << endl;

		if (nJets == _nrJets)
		{	// if the number of jets is correct: success!
			break;
		}
		else if (nJets < _nrJets)
		{	// if number of jets is too small: we need a smaller Radius per jet (so that we get more jets)
			R -= RDiff;
		}
		else if (nJets > _nrJets)
		{	// if the number of jets is too high: increase the Radius
			R += RDiff;
		}

		RDiff /= 2;

		// clean up
		if (pluginSisCone) delete pluginSisCone;
		if (pluginSisConeSph) delete pluginSisConeSph;
		delete def;

	}

	if (iIter == ITERATIVE_INCLUSIVE_MAX_ITERATIONS)
	{
		streamlog_out(WARNING) << "Maximum number of iterations reached. Canceling" << endl;
		_statsNrSkippedMaxIterations++;
		// Currently we will return the latest results, independent if the number is actually matched
		// jetsReturn.clear();
	}



	return jetsReturn;
}

EVENT::ReconstructedParticle* FastJetProcessor::getRecPar(fastjet::PseudoJet& fj, const vector< fastjet::PseudoJet >& constituents )
{
	// create a ReconstructedParticle that saves the jet
	ReconstructedParticleImpl* jet = new ReconstructedParticleImpl();

	// save the jet's parameters
	jet->setEnergy( fj.E() );
	jet->setMass( fj.m() );

	double mom[3] = {fj.px(), fj.py(), fj.pz()};
	jet->setMomentum( mom );

	// add information about the included particles
	for (unsigned int n = 0; n < constituents.size(); ++n)
	{
		ReconstructedParticle* p = dynamic_cast< ReconstructedParticle* > (_reconstructedPars->getElementAt(constituents[n].user_index())) ;
		jet->addParticle( p );

	}

	return jet;
}

void FastJetProcessor::convertFromRecParticle(LCCollection* recCol)
{
	// foreach RecoParticle in the LCCollection: convert it into a PseudoJet and and save it in our list
	for (int i = 0; i < recCol->getNumberOfElements(); ++i)
	{
		ReconstructedParticle* par = dynamic_cast<ReconstructedParticle*> (recCol->getElementAt(i));
		_pjList.push_back(
				fastjet::PseudoJet( par->getMomentum()[0],
									par->getMomentum()[1],
									par->getMomentum()[2],
									par->getEnergy() ) );
		_pjList.back().set_user_index(i);	// save the id of this recParticle
	}
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
			<< endl;

}




