#include "FastJetClustering.h"
#include <iostream>
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include "IMPL/ReconstructedParticleImpl.h"
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/PIDHandler.h"
#include "EVENT/LCFloatVec.h"
#include "TMath.h"
#include <UTIL/LCTOOLS.h>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace lcio ;
using namespace marlin ;
using namespace std ;

FastJetClustering aFastJetClustering ;

FastJetClustering::FastJetClustering() : Processor("FastJetClustering") {
  
  // Processor description
  _description = "FastJet Clustering ..." ;
 

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			   "InputCollection",
			   "Collection of reconstructed particles",
			   _inputCollection,
			   std::string("Unset") );
    
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			    "OutputCollection",
			    "Name of collection with the found jets",
			    _outputCollection,
			    std::string("Unset") );

  registerProcessorParameter("Debug",
			     "debug printout",
			     _print,
			     int(0)); 

  registerProcessorParameter("Algorithm",
			     "FastJet algorithm",
			     sAlgorithm,
			     std::string("antikt_algorithm")); 

  registerProcessorParameter("R",
			     "R Parameter",
			     _RPar,
			     double(0.7)); 

  registerProcessorParameter("EjetMin",
			     "Ejet",
			     _eJet,
			     double(10.0)); 

  registerProcessorParameter("NJets",
			     "max nb of jets",
			     _nJetMax,
			     int(25)); 

  registerProcessorParameter("FillTree",
			     "tuple",
			     _fillTree,
			     int(0)); 
}

void FastJetClustering::init() { 
  
  printParameters() ;

  if(_fillTree){

    if(_print>0)cout << "FastJetClustering: Making Tuples" << endl;
    _rootfile = new TFile("FastJetClustering.root","RECREATE");
    _rootfile->cd("");
    
    // Declaration of Tree 
    _Etree = new TTree("Events","DST Events");
    _Etree->Branch("NRun",&_nRun,"NRun/I");
    _Etree->Branch("NEvt",&_nEvt,"NEvt/I");
    _Etree->Branch("Ecms",&_eCMS,"Ecms/F");
    _Etree->Branch("NJets",&_nJets,"NJets/I");
    _Etree->Branch("NJetsHE",&_nJetsHE,"NJetsHE/I");
    _Etree->Branch("RPar",&_rp,"RPar/D");
    _Etree->Branch("jetVector",&_jetVector,"jetVector[5][50]/F");

  }

  if(sAlgorithm=="kt_algorithm"){
    SetAlgorithm(fastjet::kt_algorithm);
  }else if(sAlgorithm=="cambridge_algorithm"){ 
    SetAlgorithm(fastjet::cambridge_algorithm);
  }else if(sAlgorithm=="antikt_algorithm"){ 
    SetAlgorithm(fastjet::antikt_algorithm);
  }else if(sAlgorithm=="ee_kt_algorithm"){ 
    SetAlgorithm(fastjet::ee_kt_algorithm);
  }else if(sAlgorithm=="genkt_algorithm"){ 
    SetAlgorithm(fastjet::genkt_algorithm);
  }else if(sAlgorithm=="ee_genkt_algorithm"){ 
    SetAlgorithm(fastjet::ee_genkt_algorithm);
  }
}

void FastJetClustering::processRunHeader( LCRunHeader* run) { 

  _nRun = run->getRunNumber() ;

} 

void FastJetClustering::processEvent( LCEvent * evt ) { 

  _jetsCol= new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

  _nRun = evt->getRunNumber();
  _nEvt = evt->getEventNumber();

  if(_print>0)cout <<"Run " << _nRun << " Evt " << _nEvt << endl;

  for(int i1=0;i1<5;i1++){
    for(int i2=0;i2<50;i2++){
      _jetVector[i1][i2]=0.;
    }
  }

  LCCollection* enflowcol=evt->getCollection(_inputCollection);
  int nenflow =  enflowcol->getNumberOfElements(); 

  double px, py, pz, E;

  fastjet::JetAlgorithm algorithm = GetAlgorithm(); 

  vector<fastjet::PseudoJet> input_particles;

  for ( int ienflow=0; ienflow<nenflow ; ienflow++){
    ReconstructedParticle* enflow = dynamic_cast<ReconstructedParticle*>(enflowcol->getElementAt( ienflow ));

    px = enflow->getMomentum()[0];
    py = enflow->getMomentum()[1];
    pz = enflow->getMomentum()[2];
    E  = enflow->getEnergy();
    
    fastjet::PseudoJet thisPtc(px,py,pz,E);
    thisPtc.set_user_index(ienflow);

    input_particles.push_back(thisPtc);
  }

  _rp = _RPar;

  fastjet::Strategy strategy = fastjet::Best;

  fastjet::RecombinationScheme recomb_scheme = fastjet::E_scheme;

  float momentum[3],energy;

  if(_nJetMax>0 && _eJet>0){

    fastjet::JetDefinition jet_def_0(algorithm, _rp, recomb_scheme, strategy);

    fastjet::ClusterSequence cs_0(input_particles, jet_def_0);
  
    vector<fastjet::PseudoJet> jets_0 = cs_0.inclusive_jets();

    vector<fastjet::PseudoJet> sortedJets_0 = sorted_by_E(cs_0.inclusive_jets());

    if(_print>0)cout << "FastJetClustering: Nb of Jets " << sortedJets_0.size() << endl;
  
    _nJets = sortedJets_0.size();

    _nJetsHE=0;

    for(unsigned ij=0; ij<sortedJets_0.size();ij++){
      if(sortedJets_0[ij].e() > _eJet) _nJetsHE++;
    }
    
    while(_nJetsHE <_nJetMax && _rp>0.35 && _nJetMax>0 && _eJet>0){
      
      _rp-=0.10;
      
      fastjet::JetDefinition jet_def(algorithm, _rp, recomb_scheme, strategy); 
      
      fastjet::ClusterSequence cs(input_particles, jet_def);
      
      vector<fastjet::PseudoJet> jets = cs.inclusive_jets();
      
      vector<fastjet::PseudoJet> sortedJets = sorted_by_E(cs.inclusive_jets());
      
      _nJetsHE=0;
      
      for(unsigned ij=0; ij<sortedJets.size();ij++){
	if(sortedJets[ij].e() > _eJet) _nJetsHE++;
      }
      
      if(_print>1)cout << "RPar " << _rp << " " <<  _nJetsHE << " E_3 " << sortedJets[_nJetMax-2].e() << " E_4 " << sortedJets[_nJetMax-1].e() << " E_5 " << sortedJets[_nJetMax].e() << endl;
      
    }

  }

  fastjet::JetDefinition jet_defF(algorithm, _rp, recomb_scheme, strategy); 

  fastjet::ClusterSequence csF(input_particles, jet_defF);
  
  vector<fastjet::PseudoJet> jetsF = csF.inclusive_jets();
  
  vector<fastjet::PseudoJet> sortedJetsF = sorted_by_E(csF.inclusive_jets());
  
  int nmx = sortedJetsF.size();
  if(nmx > _nJetMax) nmx = _nJetMax;

  _nJets = sortedJetsF.size();
  _nJetsHE=0;

  for(int ij=0; ij<_nJets;ij++){

    momentum[0]= sortedJetsF[ij].px();
    momentum[1]= sortedJetsF[ij].py();
    momentum[2]= sortedJetsF[ij].pz();
    energy = sortedJetsF[ij].e();

    if(_print>1)cout << "Jet " << ij << " En " << energy << endl;

    if(energy>_eJet && energy>5. && ij < 25){
      vector<fastjet::PseudoJet> jetConstituents = csF.constituents(sortedJetsF[ij]);

      _jetVector[0][ij] = sortedJetsF[ij].px();
      _jetVector[1][ij] = sortedJetsF[ij].py();
      _jetVector[2][ij] = sortedJetsF[ij].pz();
      _jetVector[3][ij] = sortedJetsF[ij].e();
      _jetVector[4][ij] = (float)jetConstituents.size();
      
      _nJetsHE++;

      if(ij<_nJetMax){

	ReconstructedParticleImpl* Jets = new ReconstructedParticleImpl;

	for(unsigned ip = 0; ip < jetConstituents.size(); ip++){

	  if(jetConstituents[ip].user_index()>=0 && jetConstituents[ip].user_index() < enflowcol->getNumberOfElements()){
	    
	    ReconstructedParticle* enflow = dynamic_cast<ReconstructedParticle*>(enflowcol->getElementAt(jetConstituents[ip].user_index()));
	    
	    Jets->addParticle(enflow);
	  }
	}

	Jets->setMomentum(momentum);
	Jets->setEnergy(energy);
	
	_jetsCol->addElement(Jets);

      }
    }
  }

  if(_fillTree){
    _Etree->Fill();
  }    

  _jetsCol->parameters().setValue( "RPar", (float)_rp );
  _jetsCol->parameters().setValue( "NJets", (float)_nJets );

  evt->addCollection(_jetsCol ,_outputCollection) ; 
}

void  FastJetClustering::check( LCEvent* ) {

}


void FastJetClustering::end(){ 
  
  std::cout << "FastJetClustering::end()  " 
 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
 	    << std::endl ;

  if(_fillTree){
    if(_print>0)cout << "FastJetClustering: Saving Tuples" << endl;
    _rootfile->cd("");
    _rootfile->Write();
    _rootfile->Close();
  }
}
