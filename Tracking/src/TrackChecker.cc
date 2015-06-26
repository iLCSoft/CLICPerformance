/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "TrackChecker.h"
#include "DDRec/API/IDDecoder.h"

#include <marlinutil/HelixClass.h>
#include "MarlinTrk/MarlinTrkUtils.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/HelixFit.h"
#include "MarlinTrk/HelixFit.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/Factory.h"
#include "MarlinTrk/MarlinTrkDiagnostics.h"
#include "MarlinTrk/IMarlinTrkSystem.h"

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <EVENT/MCParticle.h>
#include <IMPL/TrackImpl.h>

#include <UTIL/CellIDEncoder.h>
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>
#include <UTIL/LCRelationNavigator.h>

#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/SurfaceManager.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "marlin/ProcessorEventSeeder.h"
#include "marlin/Global.h"

#include "CLHEP/Vector/TwoVector.h"

#include <AIDA/IAnalysisFactory.h>
#include <AIDA/IHistogramFactory.h>

#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <climits>
#include <cfloat>
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"

using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace DD4hep ;
using namespace AIDA ;

TrackChecker aTrackChecker;

/*

 Generic code set up to quickly add algorithms
 
*/

TrackChecker::TrackChecker() : Processor("TrackChecker") {
	
	// modify processor description
	_description = "TrackChecker plots pull distributions, track resolution parameters and reconstruction efficiencies" ;

  // Input collections
  registerInputCollection( LCIO::TRACK,
                           "TrackCollectionName",
                           "Track collection name",
                           m_inputTrackCollection,
                           std::string("SiTracks"));

  registerInputCollection( LCIO::LCRELATION,
                          "TrackRelationCollectionName",
                          "Track relation collection name",
                          m_inputTrackRelationCollection,
                          std::string("SiTrackRelations"));
  
  registerInputCollection( LCIO::MCPARTICLE,
                          "MCParticleCollectionName",
                          "Name of the MCParticle input collection",
                          m_inputParticleCollection,
                          std::string("MCParticle"));
  
  registerProcessorParameter( "OutputFileName",
                             "Name for the output histogram file",
                             m_outputFile,
                             std::string("TrackCheckerHistos.root"));
  
}


void TrackChecker::init() {
	
	// Print the initial parameters
	printParameters() ;
	
	// Reset counters
	m_runNumber = 0 ;
	m_eventNumber = 0 ;
  
  // Get the magnetic field
  DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
  const double position[3]={0,0,0}; // position to calculate magnetic field at (the origin in this case)
  double magneticFieldVector[3]={0,0,0}; // initialise object to hold magnetic field
  lcdd.field().magneticField(position,magneticFieldVector); // get the magnetic field vector from DD4hep
  m_magneticField = magneticFieldVector[2]/dd4hep::tesla; // z component at (0,0,0)

	// Register this process
	Global::EVENTSEEDER->registerProcessor(this);
  
  // Create output histograms
  m_omegaPull = new TH1F("OmegaPullPlot","OmegaPullPlot",200,-10,10);
  m_phiPull = new TH1F("PhiPullPlot","PhiPullPlot",200,-10,10);
  m_tanLambdaPull = new TH1F("TanLambdaPullPlot","TanLambdaPullPlot",200,-10,10);
  m_d0Pull = new TH1F("D0PullPlot","D0PullPlot",200,-10,10);
  m_z0Pull = new TH1F("Z0PullPlot","Z0PullPlot",200,-10,10);
  
  m_omegaMCParticle = new TH1F("OmegaMCParticle","OmegaMCParticle",200,-10,10);
  m_phiMCParticle = new TH1F("PhiMCParticle","PhiMCParticle",200,-10,10);
  m_tanLambdaMCParticle = new TH1F("TanLambdaMCParticle","TanLambdaMCParticle",200,-10,10);
  m_d0MCParticle = new TH1F("D0MCParticle","D0MCParticle",200,-10,10);
  m_z0MCParticle = new TH1F("Z0MCParticle","Z0MCParticle",200,-10,10);
  
  m_omegaTrack = new TH1F("OmegaTrack","OmegaTrack",200,-10,10);
  m_phiTrack = new TH1F("PhiTrack","PhiTrack",200,-10,10);
  m_tanLambdaTrack = new TH1F("TanLambdaTrack","TanLambdaTrack",200,-10,10);
  m_d0Track = new TH1F("D0Track","D0Track",200,-10,10);
  m_z0Track = new TH1F("Z0Track","Z0Track",200,-10,10);
  
  m_trackChi2 = new TH1F("TrackChi2","TrackChi2",500,0,10);

  m_omegaResidual = new TH1F("omegaResidual","omegaResidual",200,-0.001,0.001);
  m_phiResidual = new TH1F("phiResidual","phiResidual",100,-0.005,0.005);
  m_tanLambdaResidual = new TH1F("tanLambdaResidual","tanLambdaResidual",100,-0.005,0.005);
  m_d0Residual = new TH1F("d0Residual","d0Residual",200,-0.1,0.1);
  m_z0Residual = new TH1F("z0Residual","z0Residual",200,-0.1,0.1);
}


void TrackChecker::processRunHeader( LCRunHeader* run) {
	++m_runNumber ;
}

void TrackChecker::processEvent( LCEvent* evt ) {
	
  // Get the collection of tracks
  LCCollection* trackCollection = 0 ;
  getCollection(trackCollection, m_inputTrackCollection, evt); if(trackCollection == 0) return;

  // Get the collection of MC particles
  LCCollection* particleCollection = 0 ;
  getCollection(particleCollection, m_inputParticleCollection, evt); if(particleCollection == 0) return;

  // Get the collection of track relations to MC particles
  LCCollection* trackRelationCollection = 0 ;
  getCollection(trackRelationCollection, m_inputTrackRelationCollection, evt); if(trackRelationCollection == 0) return;
  
  // Create the relations navigator
  LCRelationNavigator* relation = new LCRelationNavigator( trackRelationCollection );
  
  /*
   Loop over all tracks, get the MC particle that produced it, and make
   a helical track fit to produce the "true" path of the particle. Then
   calculate the pull distributions and absolute error distributions.
  */
  
  // Loop over all tracks
  int nTracks = trackCollection->getNumberOfElements();
  for(int itTrack=0;itTrack<nTracks;itTrack++){

    // Get the track
    Track* track = dynamic_cast<Track*>( trackCollection->getElementAt(itTrack) ) ;
    
    // Get the related list of MC particles
    const LCObjectVec& mcparticleVector = relation->getRelatedToObjects( track );

    // Take the first MC particle (there should only be one)
    MCParticle* particle = dynamic_cast<MCParticle*>(mcparticleVector.at(0));
    
    /*
     Now apply some cuts - only valid MC particles will be considered
    */
    
//    if(particle->isDecayedInTracker()) continue; // Ignore tracks that decay inside the tracker volume
		
    
    // Construct the true helical trajectory
    HelixClass helix ;
    float pos[3]={particle->getVertex()[0],particle->getVertex()[1],particle->getVertex()[2]}; // Vertex position
    float mom[3]={particle->getMomentum()[0],particle->getMomentum()[1],particle->getMomentum()[2]}; // Particle momentum
    float charge = particle->getCharge() ; // Particle charge
    helix.Initialize_VP(pos,mom,charge,m_magneticField) ; // Initialised helix track

    // Now get the pull distributions
    double d0mcp = helix.getD0() ;
    double phmcp = helix.getPhi0() ;
    double ommcp = helix.getOmega() ;
    double z0mcp = helix.getZ0() ;
    double tLmcp = helix.getTanLambda() ;

    double d0track = track->getD0() ;
    double phtrack = track->getPhi() ;
    double omtrack = track->getOmega() ;
    double z0track = track->getZ0() ;
    double tLtrack = track->getTanLambda() ;

    double d0trackError = track->getCovMatrix()[0] ;
    double phtrackError = track->getCovMatrix()[2] ;
    double omtrackError = track->getCovMatrix()[5] ;
    double z0trackError = track->getCovMatrix()[9] ;
    double tLtrackError = track->getCovMatrix()[14] ;

    // Fill the output histograms
    m_d0Pull->Fill( (d0track-d0mcp)/sqrt(d0trackError) );
    m_phiPull->Fill( (phtrack-phmcp)/sqrt(phtrackError) );
    m_omegaPull->Fill( (omtrack-ommcp)/sqrt(omtrackError) );
    m_z0Pull->Fill( (z0track-z0mcp)/sqrt(z0trackError) );
    m_tanLambdaPull->Fill( (tLtrack-tLmcp)/sqrt(tLtrackError) );
 
    m_d0Residual->Fill( (d0track-d0mcp) );
    m_phiResidual->Fill( (phtrack-phmcp) );
    m_omegaResidual->Fill( (omtrack-ommcp) );
    m_z0Residual->Fill( (z0track-z0mcp) );
    m_tanLambdaResidual->Fill( (tLtrack-tLmcp) );

    m_omegaMCParticle->Fill(ommcp);
    m_phiMCParticle->Fill(phmcp);
    m_tanLambdaMCParticle->Fill(tLmcp);
    m_d0MCParticle->Fill(d0mcp);
    m_z0MCParticle->Fill(z0mcp);
    
    m_omegaTrack->Fill(omtrack);
    m_phiTrack->Fill(phtrack);
    m_tanLambdaTrack->Fill(tLtrack);
    m_d0Track->Fill(d0track);
    m_z0Track->Fill(z0track);
    
    m_trackChi2->Fill(track->getChi2());

//    if(track->getChi2() < 0.05) std::cout<<"EVENT NUMBER "<<m_eventNumber<<std::endl;
  }
  
	// Increment the event number
	m_eventNumber++ ;
	delete relation;
	
}

void TrackChecker::check( LCEvent * evt ) {
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}


void TrackChecker::end(){

	streamlog_out(MESSAGE) << " end()  " << name()
	<< " processed " << m_eventNumber << " events in " << m_runNumber << " runs "
	<< std::endl ;
	
  // Write output root file with histograms
  TFile* output = new TFile(m_outputFile.c_str(),"RECREATE");
  m_omegaPull->Write();
  m_phiPull->Write();
  m_tanLambdaPull->Write();
  m_d0Pull->Write();
  m_z0Pull->Write();
	
	TCanvas* canv = new TCanvas();
	canv->Divide(3,2);
	canv->cd(1);
	m_omegaPull->Draw();
	m_omegaPull->Fit("gaus","");
	gStyle->SetOptFit(1111);
	canv->cd(2);
	m_phiPull->Draw();
	m_phiPull->Fit("gaus","");
	gStyle->SetOptFit(1111);
	canv->cd(3);
	m_tanLambdaPull->Draw();
	m_tanLambdaPull->Fit("gaus","");
	gStyle->SetOptFit(1111);
	canv->cd(4);
	m_d0Pull->Draw();
	m_d0Pull->Fit("gaus","");
	gStyle->SetOptFit(1111);
	canv->cd(5);
	m_z0Pull->Draw();
	m_z0Pull->Fit("gaus","");
	gStyle->SetOptFit(1111);
	canv->Write();

  m_omegaResidual->Write();
  m_phiResidual->Write();
  m_tanLambdaResidual->Write();
  m_d0Residual->Write();
  m_z0Residual->Write();

  m_omegaMCParticle->Write();
  m_phiMCParticle->Write();
  m_tanLambdaMCParticle->Write(); 
  m_d0MCParticle->Write();
  m_z0MCParticle->Write();
  
  m_omegaTrack->Write();
  m_phiTrack->Write();
  m_tanLambdaTrack->Write();
  m_d0Track->Write();
  m_z0Track->Write();
  
  m_trackChi2->Write();

  output->Close();

}

void TrackChecker::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
  try{
    collection = evt->getCollection( collectionName ) ;
  }
  catch(DataNotAvailableException &e){
    std::cout<<"- cannot get collections !!"<<std::endl;
    streamlog_out(DEBUG4) << "Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}


  
  
  
  
  
  
  
