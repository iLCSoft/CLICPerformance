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
#include <UTIL/LCTrackerConf.h>
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
#include "TTree.h"

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
                           std::string("RefittedTracks"));

  registerInputCollection( LCIO::LCRELATION,
                           "TrackRelationCollectionName",
                           "Track relation collection name",
                           m_inputTrackRelationCollection,
                           std::string("RefittedTracksMCTruthLink"));
  
  registerInputCollection( LCIO::MCPARTICLE,
                           "MCParticleCollectionName",
                           "Name of the MCParticle input collection",
                           m_inputParticleCollection,
                           std::string("MCParticle"));
  
  registerProcessorParameter("UseOnlyTree",
                             "Flag to fill only tree variables or also histogram version of pulls and residuals",
                             m_useOnlyTree,
                             bool(true));
  
  registerProcessorParameter( "TreeName",
                             "Name of the root tree",
                             m_treeName,
                             std::string("perftree"));
  
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
  

}


void TrackChecker::processRunHeader( LCRunHeader* ) {
	++m_runNumber ;
}

void TrackChecker::processEvent( LCEvent* evt ) {

  //std::cout<< "---- event << " << m_eventNumber << std::endl;

  if( isFirstEvent() ) { 
    perftree = new TTree(m_treeName.c_str(),m_treeName.c_str());

    int bufsize = 32000; //default buffer size 32KB
          
    perftree->Branch("truePt","std::vector<double >",&truePt,bufsize,0) ;
    perftree->Branch("trueTheta","std::vector<double >",&trueTheta,bufsize,0) ;
    perftree->Branch("truePhi","std::vector<double >",&truePhi,bufsize,0) ;
    perftree->Branch("trueD0","std::vector<double >",&trueD0,bufsize,0) ;
    perftree->Branch("trueZ0","std::vector<double >",&trueZ0,bufsize,0) ;
    perftree->Branch("trueP","std::vector<double >",&trueP,bufsize,0) ;

    perftree->Branch("recoPt","std::vector<double >",&recoPt,bufsize,0) ;
    perftree->Branch("recoTheta","std::vector<double >",&recoTheta,bufsize,0) ;
    perftree->Branch("recoPhi","std::vector<double >",&recoPhi,bufsize,0) ;
    perftree->Branch("recoD0","std::vector<double >",&recoD0,bufsize,0) ;
    perftree->Branch("recoZ0","std::vector<double >",&recoZ0,bufsize,0) ;
    perftree->Branch("recoP","std::vector<double >",&recoP,bufsize,0) ;

    perftree->Branch("recoNhits","std::vector<int >",&recoNhits,bufsize,0) ;
    perftree->Branch("recoChi2OverNDF","std::vector<double >",&recoChi2OverNDF,bufsize,0) ;
    perftree->Branch("recoMinDist","std::vector<double >",&recoMinDist,bufsize,0) ;

    //   if(m_useOnlyTree) {
      perftree->Branch("pullOmega","std::vector<double >",&pullOmega,bufsize,0) ;
      perftree->Branch("pullPhi","std::vector<double >",&pullPhi,bufsize,0) ;
      perftree->Branch("pullTanLambda","std::vector<double >",&pullTanLambda,bufsize,0) ;
      perftree->Branch("pullD0","std::vector<double >",&pullD0,bufsize,0) ;
      perftree->Branch("pullZ0","std::vector<double >",&pullZ0,bufsize,0) ;

      perftree->Branch("resOmega","std::vector<double >",&resOmega,bufsize,0) ;
      perftree->Branch("resPhi","std::vector<double >",&resPhi,bufsize,0) ;
      perftree->Branch("resTanLambda","std::vector<double >",&resTanLambda,bufsize,0) ;
      perftree->Branch("resD0","std::vector<double >",&resD0,bufsize,0) ;
      perftree->Branch("resZ0","std::vector<double >",&resZ0,bufsize,0) ;
      //    } else {
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
      //    }

  }//end is first event


  clearEventVariables();

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
  
  //TEST - just try to run on single muons with only 1 mcparticle
  //int nmcp =  particleCollection->getNumberOfElements();
  //if (nmcp==1) {


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
    float pos[3]={static_cast<float>(particle->getVertex()[0]),static_cast<float>(particle->getVertex()[1]),static_cast<float>(particle->getVertex()[2])}; // Vertex position
    float mom[3]={static_cast<float>(particle->getMomentum()[0]),static_cast<float>(particle->getMomentum()[1]),static_cast<float>(particle->getMomentum()[2])}; // Particle momentum
    float charge = particle->getCharge() ; // Particle charge
    helix.Initialize_VP(pos,mom,charge,m_magneticField) ; // Initialised helix track

    // Now get the pull distributions
    double d0mcp = helix.getD0() ;
    double phmcp = helix.getPhi0() ;
    double ommcp = helix.getOmega() ;
    double z0mcp = helix.getZ0() ;
    double tLmcp = helix.getTanLambda() ;



    streamlog_out(MESSAGE) << " ---- reference point of track fit x,y,z =  " << track->getReferencePoint()[0] << " , " << track->getReferencePoint()[1] << " , " << track->getReferencePoint()[2] << std::endl;



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



    //angle conversions
    angleInFixedRange(phmcp);
	  angleInFixedRange(phtrack);    

    //    if (m_useOnlyTree) {

      resOmega.push_back(omtrack-ommcp);
      resPhi.push_back(phtrack-phmcp);
      resTanLambda.push_back(tLtrack-tLmcp);
      resD0.push_back(d0track-d0mcp);
      resZ0.push_back(z0track-z0mcp);

      pullOmega.push_back(resOmega.back()/sqrt(omtrackError)); 
      pullPhi.push_back(resPhi.back()/sqrt(phtrackError)); 
      pullTanLambda.push_back(resTanLambda.back()/sqrt(tLtrackError)); 
      pullD0.push_back(resD0.back()/sqrt(d0trackError)); 
      pullZ0.push_back(resZ0.back()/sqrt(z0trackError)); 

      //   } else {
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
      //   }

	  truePt.push_back(0.3*m_magneticField/fabs(ommcp*1000.)); //omega in 1/mm -> transformed in 1/m
    double trueThetaRad = M_PI/2 - atan(tLmcp); 
	  trueTheta.push_back(trueThetaRad*180./M_PI);
	  truePhi.push_back(phmcp*180./M_PI);
	  trueD0.push_back(d0mcp);
	  trueZ0.push_back(z0mcp);
	  trueP.push_back(truePt.back()/sin(trueThetaRad)); 

	  recoPt.push_back(0.3*m_magneticField/(fabs(omtrack)*1000.));
    double recoThetaRad =  M_PI/2 - atan(tLtrack);
	  recoTheta.push_back(recoThetaRad*180./M_PI);
	  recoPhi.push_back(phtrack*180./M_PI);
	  recoD0.push_back(d0track);
	  recoZ0.push_back(z0track);
	  recoP.push_back(recoPt.back()/sin(recoThetaRad)); 

    recoNhits.push_back(track->getTrackerHits().size());
    if (track->getNdf()>0) recoChi2OverNDF.push_back(track->getChi2()/track->getNdf());
    else recoChi2OverNDF.push_back(-1.);

    double minDist = -1.;
    for(int jTrack=itTrack+1;jTrack<nTracks;jTrack++){
      Track* trackB = dynamic_cast<Track*>( trackCollection->getElementAt(jTrack) ) ;
      double distAB = getDistanceFromClosestHit(track, trackB);
      if (minDist<0. || minDist>distAB) minDist = distAB;
    }
    recoMinDist.push_back(minDist);
  } //end loop on tracks
  
  perftree->Fill();

  //} // END TEST

	// Increment the event number
	m_eventNumber++ ;
	delete relation;
	
}



void TrackChecker::check( LCEvent * ) {
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}



void TrackChecker::end(){

	streamlog_out(MESSAGE) << " end()  " << name()
                         << " processed " << m_eventNumber << " events in " << m_runNumber << " runs "
                         << std::endl ;
	

	
  gStyle->SetOptStat(111111);


  //  if (!m_useOnlyTree) {



    m_omegaPull->Write();
    m_phiPull->Write();
    m_tanLambdaPull->Write();
    m_d0Pull->Write();
    m_z0Pull->Write();
 
    pulls = new TCanvas("pulls","Pulls of the track parameters",800,800);
    pulls->Divide(3,2);
    pulls->cd(1);
    m_omegaPull->Draw();
    m_omegaPull->Fit("gaus","");
    gStyle->SetOptFit(1111);
    pulls->cd(2);
    m_phiPull->Draw();
    m_phiPull->Fit("gaus","");
    gStyle->SetOptFit(1111);
    pulls->cd(3);
    m_tanLambdaPull->Draw();
    m_tanLambdaPull->Fit("gaus","");
    gStyle->SetOptFit(1111);
    pulls->cd(4);
    m_d0Pull->Draw();
    m_d0Pull->Fit("gaus","");
    gStyle->SetOptFit(1111);
    pulls->cd(5);
    m_z0Pull->Draw();
    m_z0Pull->Fit("gaus","");
    gStyle->SetOptFit(1111);
    pulls->Write();

    res = new TCanvas("res","Residuals of the track parameters",800,800);
    res->Divide(3,2);
    res->cd(1);
    m_omegaResidual->Draw();
    res->cd(2);
    m_phiResidual->Draw();
    res->cd(3);
    m_tanLambdaResidual->Draw();
    res->cd(4);
    m_d0Residual->Draw();
    res->cd(5);
    m_z0Residual->Draw();
    res->Write();

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
    //  }

  perftree->Write();

}




void TrackChecker::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
  try{
    collection = evt->getCollection( collectionName ) ;
  }
  catch(DataNotAvailableException &e){
    std::cout<<"- cannot get collections !!"<<std::endl;
    streamlog_out(ERROR) << "Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}
  


double TrackChecker::getDistanceFromClosestHit(Track*& trackA, Track*& trackB){

  double min_dist = -1.;
  
  TrackerHitVec vecHitsA = trackA->getTrackerHits();
  TrackerHitVec vecHitsB = trackB->getTrackerHits();

  for(size_t i=0; i<vecHitsA.size(); i++){
    for(size_t j=0; j<vecHitsB.size(); j++){
      double dist = getDist(vecHitsA.at(i),vecHitsB.at(j));
      if (min_dist<0 || min_dist>dist) min_dist = dist;
    }
  }

  return min_dist;
}


double TrackChecker::getDist(TrackerHit*& hitA, TrackerHit*& hitB){

  double dx = hitA->getPosition()[0] - hitB->getPosition()[0];
  double dy = hitA->getPosition()[1] - hitB->getPosition()[1];
  double dz = hitA->getPosition()[2] - hitB->getPosition()[2];
  double dist = sqrt(dx*dx + dy*dy + dz*dz);
  
  return dist;
}


void TrackChecker::angleInFixedRange(double& angle){
  
  while (angle <= -M_PI ) angle = angle + 2*M_PI;
  while (angle >   M_PI ) angle = angle - 2*M_PI;
  
  return;
}
 
  

  
void TrackChecker::clearEventVariables(){

  truePt.clear();
  trueTheta.clear();
  truePhi.clear();
  trueD0.clear();
  trueZ0.clear();
  trueP.clear();

  recoPt.clear();
  recoTheta.clear();
  recoPhi.clear();
  recoD0.clear();
  recoZ0.clear();
  recoP.clear();

  recoNhits.clear();
  recoChi2OverNDF.clear();
  recoMinDist.clear();

  pullOmega.clear();
  pullPhi.clear();
  pullTanLambda.clear();
  pullD0.clear();
  pullZ0.clear();

  resOmega.clear();
  resPhi.clear();
  resTanLambda.clear();
  resD0.clear();
  resZ0.clear();

  return;

}  
  
