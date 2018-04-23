/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "ClicEfficiencyCalculator.h"

#include "marlin/AIDAProcessor.h"

#include "MarlinTrk/HelixTrack.h"
#include <marlinutil/HelixClass.h>

#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTrackerConf.h>
#include <UTIL/LCRelationNavigator.h>

#include <marlinutil/GeometryUtil.h>

#include <AIDA/IAnalysisFactory.h>
#include <AIDA/IHistogramFactory.h>

#include <TLorentzVector.h>
#include <TTree.h>


using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace AIDA ;

ClicTrackingValidate aClicTrackingValidate;

/* 

   Processor to validate the tracking. Uses MC information at particle and hit level.

*/

ClicTrackingValidate::ClicTrackingValidate() : Processor("ClicTrackingValidate") {

  // processor description
  _description = "ClicTrackingValidate makes use of MC information at particle and hit level and compares the reconstructed hits and tracks to validate the tracking";

  // Input collections
  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticleCollectionName",
			   _inputParticleCollection,
			   std::string("MCParticle"));

  std::vector<std::string> inputSimTrackerHitCollections;
  inputSimTrackerHitCollections.push_back(std::string("VertexBarrelHits"));

  registerInputCollections( LCIO::SIMTRACKERHIT,
			   "SimTrackerHitCollectionNames",
			   "Name of the SimTrackerHit input collections",
			   _inputSimTrackerHitCollections,
			   inputSimTrackerHitCollections);

  std::vector<std::string> inputTrackerHitCollections;
  inputTrackerHitCollections.push_back(std::string("VXDTrackerHits"));

  registerInputCollections( LCIO::TRACKERHITPLANE,
			    "TrackerHitCollectionNames",
			    "Name of the TrackerHit input collections",
			    _inputTrackerHitCollections,
			    inputTrackerHitCollections);

  std::vector<std::string> inputTrackerHitRelationCollections;
  inputTrackerHitRelationCollections.push_back(std::string("VXDTrackerHitRelations"));

  registerInputCollections( LCIO::LCRELATION,
			    "SimTrackerHitRelCollectionNames",
			    "Name of the TrackerHit - SimTrackerHit relation collections",
			    _inputTrackerHitRelationCollections,
			    inputTrackerHitRelationCollections);

  registerInputCollection( LCIO::TRACK,
			   "SiTrackCollectionName",
			   "Name of the silicon track collection",
			   _inputTrackCollection,
			   std::string("SiTracks"));

  // Output ntuple names
  registerProcessorParameter( "trackingTreeName",
			      "Name of the tracking validation tree",
			      _trackingTreeName,
			      std::string("trackingTree"));


}

void ClicTrackingValidate::init(){
  
  // Print the initial parameters
  printParameters();
  
  // Reset counters
  _runNumber = 0;
  _eventNumber = 0;

  // Get the magnetic field
  _magneticField = MarlinUtil::getBzAtOrigin();

  // Initialise histograms
  AIDAProcessor::histogramFactory(this);

  // Tree for tracking validation information
  _trackingTree = new TTree(_trackingTreeName.c_str(),_trackingTreeName.c_str());
  int bufsize = 32000;
  //_trackingTree->Branch();


}

void ClicTrackingValidate::processRunHeader( LCRunHeader* ){
  ++_runNumber;
}

void ClicTrackingValidate::processEvent( LCEvent* evt ){

  streamlog_out( DEBUG8 ) << "Processing event " << _eventNumber << std::endl;

  // First pick up all of the collections that will be used - MCParticles, tracks, simhits and recohits from relevant subdetectors and their relations

  // Initialise CELLID encoder to get the subdetector ID from a hit
  UTIL::BitField64f _encoder(lcio::LCTrackerCellID::encoding_string());
  std::vector<int> _collections;

  // Get the collection of MCParticles
  LCCollection *particleCollection = 0;
  getCollection(particleCollection, _inputParticleCollection, evt); if(particleCollection == 0) return;

  // Get the collection of tracks
  LCCollection *trackCollection = 0;
  getCollection(trackCollection, _inputTrackCollection, evt); if(trackCollection == 0) return;

  // Make objects to hold all of the tracker hit and relation collections
  std::map<int,LCCollection*> simTrackerHitCollections;
  std::map<int,LCCollection*> trackerHitCollections;
  std::map<int,LCCollection*> trackerHitRelationCollections;
  std::map<int,std::shared_ptr<LCRelationNavigator>> relations;
  std::map<int,bool> activeDetectors;

  // Loop over each sim input collection and get the data
  for(unsigned int collection = 0; collection < _inputSimTrackerHitCollections; collection++){
    
    LCCollection *simTrackerHitCollection = 0;
    getCollection(simTrackerHitCollection, _inputSimTrackerHitCollections[collection], evt); if(simTrackerHitCollection == 0) continue;

    // If there are no hits in the simTrackerHit collection, continue
    if(simTrackerHitCollection->getNumberOfElements() == 0) continue;
    
  }

  // Loop over each reco and relation collection and get the data
  for(unsigned int collection = 0; collection < _inputTrackerHitCollections; collection++){
    
    // Get the collection of tracker hits
    LCCollection *trackerHitCollection = 0;
    getCollection(trackerHitCollection, _inputTrackerHitCollections[collection], evt); if(trackerHitCollection == 0) continue;
    
    // If there are no hits in the trackerHit collection, continue
    if(trackerHitCollection->getNumberOfElements() == 0) continue;

    // Get the collection of tracker hit relations
    LCCollection *trackerHitRelationCollection = 0;
    getCollection(trackerHitRelationCollection, _inputTrackerHitRelationCollections[collection], evt);

    // Create the relations navigator
    std::shared_ptr<LCRelationNavigator> relation = std::make_shared<LCRelationNavigator>(LCRelationNavigator(trackerHitRelationCollection));

    // Get the subdetector for this set of collections
    int subdetector = getSubdetector(trackerHitCollection,_encoder);
    activeDetectors[subdetectors] = true;

    // Save all of the created collections
    trackerHitCollections[subdetector] = trackerHitCollection;
    trackerHitRelationCollections[subdetector] = trackerHitRelationCollection;
    relations[subdetector] = relation;
    _collections.push_back(subdetector);
    
  }

}

void ClicTrackingValidate::check( LCEvent *) {
  //nothing to check at the moment
}

void ClicTrackingValidate::end(){

  streamlog_out( MESSAGE ) << " end() " << name() << " processed " << _eventNumber << " events in " << _runNumber << " runs " << std::endl;

}

void ClicTrackingValidate::getCollection(LCCollection*& collection, std::string collectionName, LCEvent* evt){

  try{
    collection = evt->getCollection(collectionName);
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG5 ) << " - cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

int ClicTrackingValidate::getSubdetector(LCCollection* collection, UTIL::BitField64 &encoder){
  TrackerHitPlane* hit = static_cast<TrackerHitPlane*>(collection->getElementAt(0));
  const int celId = hit->getCellID0();
  encoder.setValue(celId);
  int subdet = encoder[lcio::LCTrackerCellID::subdet()];
  return subdet;
}

