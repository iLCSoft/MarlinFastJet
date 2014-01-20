/*
 * FastJetProcessor.h
 *
 *  Created on: 25.05.2010
 *      Author: Lars Weuste (MPP Munich) - weuste@mpp.mpg.de
 *      		iterative inclusive algorithm based on design by Marco Battaglia (CERN) - Marco.Battaglia@cern.ch
 */

#ifndef FASTJETPROCESSOR_H_
#define FASTJETPROCESSOR_H_

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

// The enum, name and value of the enum for the Cluster Mode
enum EClusterMode
{
	NONE = 0,
	FJ_exclusive_yCut = 1,		// exclusive clustering mode implemented in FastJet
	FJ_exclusive_nJets = 2,		// exclusive clustering mode implemented in FastJet
	FJ_inclusive = 4,		// inclusive "-"
	OWN_inclusiveIteration = 8	// use FJ inclusive Clustering, but iterate until we have the desired number of jets
};

class FastJetProcessor : Processor
{
public:
	FastJetProcessor();
	virtual ~FastJetProcessor();

	virtual Processor* newProcessor()
	{
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



private:

	// the LC Collection names for input/output
	string	_lcParticleInName;
	string	_lcParticleOutName;
	string	_lcJetOutName;


	// the cluster modes
	StringVec		_clusterModeNameAndParam;
	string			_clusterModeName;
	EClusterMode	_clusterMode;
	unsigned		_nrJets;
	double			_yCut;
	double			_minPt;
	double			_minE;
	void initClusterMode();

	// list of pseudo particles / jets
	typedef vector< fastjet::PseudoJet > pseudojetList;
	pseudojetList	_pjList;

	// helper functions to init the jet algorithms in general
	void initJetAlgo();
	bool isJetAlgo(string algo, int nrParams, int supportedModes);
	StringVec					_jetAlgoNameAndParams;
	string						_jetAlgoName;
	fastjet::JetDefinition*	 	_jetAlgo;
	fastjet::JetAlgorithm		_jetAlgoType;

	void initRecoScheme();
	string							_jetRecoSchemeName;
	fastjet::RecombinationScheme 	_jetRecoScheme;

	void initStrategy();
	string							_strategyName;
	fastjet::Strategy				_strategy;


	vector< fastjet::PseudoJet > doIterativeInclusiveClustering();

	void convertFromRecParticle(LCCollection* recCol);

	EVENT::ReconstructedParticle* getRecPar(fastjet::PseudoJet& fj, const vector< fastjet::PseudoJet >& constituents);

	LCCollection*	_reconstructedPars;

	int	_statsFoundJets;
	int _statsNrEvents;
	int _statsNrSkippedEmptyEvents;
	int _statsNrSkippedFixedNrJets;
	int _statsNrSkippedMaxIterations;
	bool _storeParticlesInJets;

};

ostream& operator<<(ostream& out, EClusterMode& m);

#endif /* FASTJETPROCESSOR_H_ */
