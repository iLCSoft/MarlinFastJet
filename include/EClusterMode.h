#ifndef EClusterMode_h
#define EClusterMode_h 1

// The enum, name and value of the enum for the Cluster Mode
enum EClusterMode {
	NONE = 0,
	FJ_exclusive_yCut = 1,		// exclusive clustering mode implemented in FastJet
	FJ_exclusive_nJets = 2,		// exclusive clustering mode implemented in FastJet
	FJ_inclusive = 4,		// inclusive "-"
	OWN_inclusiveIteration = 8	// use FJ inclusive Clustering, but iterate until we have the desired number of jets
};

#endif // EClusterMode_h
