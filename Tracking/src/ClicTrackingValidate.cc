/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "ClicTrackingValidate.h"

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
                           "Name of the MCParticle input collection",
                           _inputParticleCollection,
                           std::string("MCParticle"));

  std::vector<std::string> inputSimTrackerHitCollections;
  inputSimTrackerHitCollections.push_back(std::string("VertexBarrelCollection"));
  inputSimTrackerHitCollections.push_back(std::string("VertexEndcapCollection"));
  inputSimTrackerHitCollections.push_back(std::string("InnerTrackerBarrelCollection"));
  inputSimTrackerHitCollections.push_back(std::string("InnerTrackerEndcapCollection"));
  inputSimTrackerHitCollections.push_back(std::string("OuterTrackerBarrelCollection"));
  inputSimTrackerHitCollections.push_back(std::string("OuterTrackerEndcapCollection"));

  registerInputCollections( LCIO::SIMTRACKERHIT,
                            "SimTrackerHitCollectionNames",
                            "Name of the SimTrackerHit input collections",
                            _inputSimTrackerHitCollections,
                            inputSimTrackerHitCollections);

  std::vector<std::string> inputTrackerHitCollections;
  inputTrackerHitCollections.push_back(std::string("VXDTrackerHits"));
  inputTrackerHitCollections.push_back(std::string("VXDEndcapTrackerHits"));
  inputTrackerHitCollections.push_back(std::string("ITrackerHits"));
  inputTrackerHitCollections.push_back(std::string("OTrackerHits"));
  inputTrackerHitCollections.push_back(std::string("ITrackerEndcapHits"));
  inputTrackerHitCollections.push_back(std::string("OTrackerEndcapHits"));

  registerInputCollections( LCIO::TRACKERHITPLANE,
                            "TrackerHitCollectionNames",
                            "Name of the TrackerHit input collections",
                            _inputTrackerHitCollections,
                            inputTrackerHitCollections);

  std::vector<std::string> inputTrackerHitRelationCollections;
  inputTrackerHitRelationCollections.push_back(std::string("VXDTrackerHitRelations"));
  inputTrackerHitRelationCollections.push_back(std::string("VXDEndcapTrackerHitRelations"));
  inputTrackerHitRelationCollections.push_back(std::string("InnerTrackerBarrelHitsRelations"));
  inputTrackerHitRelationCollections.push_back(std::string("OuterTrackerBarrelHitsRelations"));
  inputTrackerHitRelationCollections.push_back(std::string("InnerTrackerEndcapHitsRelations"));
  inputTrackerHitRelationCollections.push_back(std::string("OuterTrackerEndcapHitsRelations"));

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
  _trackingTree->Branch("mc_pdg","std::vector<int>",&m_vec_pdg,bufsize,0);
  _trackingTree->Branch("mc_pt","std::vector<double>",&m_vec_pt,bufsize,0);
  _trackingTree->Branch("mc_theta","std::vector<double>",&m_vec_theta,bufsize,0);
  _trackingTree->Branch("mc_phi","std::vector<double>",&m_vec_phi,bufsize,0);
  _trackingTree->Branch("mc_vertexR","std::vector<double>",&m_vec_vertexR,bufsize,0);
  _trackingTree->Branch("mc_nhits","std::vector<int>",&m_vec_mchits,bufsize,0);
  _trackingTree->Branch("trk_nhits","std::vector<int>",&m_vec_trkhits,bufsize,0);
  _trackingTree->Branch("trk_purity","std::vector<double>",&m_vec_purity,bufsize,0);
  _trackingTree->Branch("trk_misassocHits","std::vector<int>",&m_vec_misassocHits,bufsize,0);
  _trackingTree->Branch("trk_misassocFirstHit_pos","std::vector<double>",&m_vec_misassocFirstHitPos,bufsize,0);
  _trackingTree->Branch("trk_missedHits","std::vector<int>",&m_vec_missedHits,bufsize,0);
  _trackingTree->Branch("trk_missedFirstHit_pos","std::vector<double>",&m_vec_missedFirstHitPos,bufsize,0);
  

}

void ClicTrackingValidate::processRunHeader( LCRunHeader* ){
  ++_runNumber;
}

void ClicTrackingValidate::processEvent( LCEvent* evt ){

  streamlog_out( DEBUG8 ) << "Processing event " << _eventNumber << std::endl;

  // First pick up all of the collections that will be used - MCParticles, tracks, simhits and recohits from relevant subdetectors and their relations

  // Initialise CELLID encoder to get the subdetector ID from a hit
  UTIL::BitField64 encoder(lcio::LCTrackerCellID::encoding_string());
  std::vector<int> sim_collections;
  std::vector<int> reco_collections;

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
  for(unsigned int collection = 0; collection < _inputSimTrackerHitCollections.size(); collection++){
    
    // Get the collection of sim tracker hits
    LCCollection *simTrackerHitCollection = 0;
    getCollection(simTrackerHitCollection, _inputSimTrackerHitCollections[collection], evt); if(simTrackerHitCollection == 0) continue;

    // If there are no hits in the simTrackerHit collection, continue
    if(simTrackerHitCollection->getNumberOfElements() == 0) continue;

    // Get the subdetector for this set of collections
    int subdetector = getSubdetectorFromSim(simTrackerHitCollection,encoder);

    // Save the created collections
    simTrackerHitCollections[subdetector] = simTrackerHitCollection;
    sim_collections.push_back(subdetector);
    
  }

  // Loop over each reco and relation collection and get the data
  for(unsigned int collection = 0; collection < _inputTrackerHitCollections.size(); collection++){
    
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
    int subdetector = getSubdetector(trackerHitCollection,encoder);
    activeDetectors[subdetector] = true;

    // Save all of the created collections
    trackerHitCollections[subdetector] = trackerHitCollection;
    trackerHitRelationCollections[subdetector] = trackerHitRelationCollection;
    relations[subdetector] = relation;
    reco_collections.push_back(subdetector);
    
  }

  // Containers to store MCParticle information
  std::map<MCParticle*,TLorentzVector> MCParticleProperties;
  std::map<MCParticle*,std::vector<SimTrackerHit*>> MCParticleHits;
  std::map<MCParticle*,std::vector<SimTrackerHit*>>::const_iterator iter;

  //Loop over MCParticle
  int nMC = particleCollection->getNumberOfElements();

  for(int iMC = 0; iMC < nMC; iMC++){
    
    MCParticle *particle = static_cast<MCParticle*> (particleCollection->getElementAt(iMC));
    TLorentzVector mc_helper;
    mc_helper.SetPxPyPzE(particle->getMomentum()[0],particle->getMomentum()[1],particle->getMomentum()[2],particle->getEnergy());
    MCParticleProperties[particle] = mc_helper;

    std::vector<SimTrackerHit*> hitsInParticle;

    // Loop over SimTrackerHits and store the ones belonging to the particle
    for(unsigned int itCollection = 0; itCollection < sim_collections.size(); itCollection++){

      // Get the collection index (subdetector)
      int collection = sim_collections[itCollection];
      int nSimHits = simTrackerHitCollections[collection]->getNumberOfElements();

      for(int itSimHit = 0; itSimHit < nSimHits; itSimHit++){
        // Get the simHit
        SimTrackerHit* simHit = static_cast<SimTrackerHit*>(simTrackerHitCollections[collection]->getElementAt(itSimHit));

        // Get the particle to which the hit belongs
        MCParticle *mcp = simHit->getMCParticle();

        // Fill the 
        if(mcp == particle){
          hitsInParticle.push_back(simHit);
        }
      }
    }

    MCParticleHits[particle] = hitsInParticle;
    
  } // end loop MC


  // Loop over tracks and get associated MCParticle
  int nTracks = trackCollection->getNumberOfElements();

  for(int iTrack = 0; iTrack < nTracks; iTrack++){

    // Containers to store track information
    std::vector<MCParticle*> particlesOnTrack;
    std::map<MCParticle*,int> hitsParticlesOnTrack;

    Track* track = static_cast<Track*> (trackCollection->getElementAt(iTrack));
    const TrackerHitVec& hitVector = track->getTrackerHits();
    int nHits = hitVector.size();
  
    // Loop over tracker hits and store the MCParticle with hits on it
    if(nHits == 0) continue;
    for(int itHit = 0; itHit < nHits; itHit++){
      // Get the tracker hit
      TrackerHitPlane* hit = static_cast<TrackerHitPlane*>(hitVector.at(itHit));
      // Get the subdetector number
      int subdetector = getSubdetector(hit,encoder);
      // Get the simulated hit
      const LCObjectVec& simHitVector = relations[subdetector]->getRelatedToObjects(hit);
      // Take the first hit only
      SimTrackerHit* simHit = static_cast<SimTrackerHit*>(simHitVector.at(0));
      // Get the particle belonging to that hit
      MCParticle *particle = simHit->getMCParticle();
      // Check if this particle already has hits on this track
      if(std::find(particlesOnTrack.begin(),particlesOnTrack.end(),particle) == particlesOnTrack.end()) particlesOnTrack.push_back(particle);
      // Increment the number of hit counts
      hitsParticlesOnTrack[particle]++; 
    }

    // Find the associated particle
    int nParticlesOnTrack = particlesOnTrack.size();
    int maxHits = 0;
    MCParticle* associatedParticle = 0;
    int countMisassocHits = 0;
    for(int itParticle = 0; itParticle < nParticlesOnTrack; itParticle++){
      int nHitsParticleOnTrack = hitsParticlesOnTrack[particlesOnTrack[itParticle]];
      if(nHitsParticleOnTrack > maxHits){
        maxHits = nHitsParticleOnTrack;
        associatedParticle = particlesOnTrack[itParticle];
      }
      if(itParticle > 0){
        countMisassocHits += nHitsParticleOnTrack;
      }
    }


    // Reloop over tracker hits and find the first not belonging to the associated particle. Store its position (2D radius)
    double misassocFirstHitPos = -1;
    for(int itHit = 0; itHit < nHits; itHit++){

      TrackerHitPlane* hit = static_cast<TrackerHitPlane*>(hitVector.at(itHit));
      int subdetector = getSubdetector(hit,encoder);
      const LCObjectVec& simHitVector = relations[subdetector]->getRelatedToObjects(hit);
      SimTrackerHit* simHit = static_cast<SimTrackerHit*>(simHitVector.at(0));
      MCParticle *particle = simHit->getMCParticle();

      if(particle!=associatedParticle){
        misassocFirstHitPos = sqrt(pow(hit->getPosition()[0],2)+pow(hit->getPosition()[1],2));
        break;
      }

    }

    // Loop over simhits of the associated particle on track and compare with the simhits of the associated particles. Find the first mismatch. Store its position (2D radius)
    double missedFirstHitPos = -1;
    std::vector<SimTrackerHit*> simHitsAssociatedParticleOnTrack;
      
    for(int itHit = 0; itHit < nHits; itHit++){

      TrackerHitPlane* hit = static_cast<TrackerHitPlane*>(hitVector.at(itHit));
      int subdetector = getSubdetector(hit,encoder);
      const LCObjectVec& simHitVector = relations[subdetector]->getRelatedToObjects(hit);
      SimTrackerHit* simHit = static_cast<SimTrackerHit*>(simHitVector.at(0));
      simHitsAssociatedParticleOnTrack.push_back(simHit);

    }

    if(isSubsetOrEqual(MCParticleHits[associatedParticle],simHitsAssociatedParticleOnTrack)){
      continue;
    }
    else{
      auto mismatchedHits = std::mismatch(simHitsAssociatedParticleOnTrack.begin(),simHitsAssociatedParticleOnTrack.end(),MCParticleHits[associatedParticle].begin());
      SimTrackerHit *firstMismatchedHit = *mismatchedHits.second;
      missedFirstHitPos = sqrt(pow(firstMismatchedHit->getPosition()[0],2)+pow(firstMismatchedHit->getPosition()[1],2));
    }
    double purity = (double)maxHits/(double)nHits;


    int countMissedHits = MCParticleHits[associatedParticle].size() - nHits;
    
    int pdg = associatedParticle->getPDG();
    TLorentzVector mc_helper;
    mc_helper.SetPxPyPzE(associatedParticle->getMomentum()[0],associatedParticle->getMomentum()[1],associatedParticle->getMomentum()[2],associatedParticle->getEnergy());
    double pt = mc_helper.Pt();
    double theta = mc_helper.Theta();
    double phi = mc_helper.Phi();
    double vertexR = sqrt(pow(associatedParticle->getVertex()[0],2)+pow(associatedParticle->getVertex()[1],2));

    m_vec_pdg.push_back(pdg);
    m_vec_pt.push_back(pt);
    m_vec_theta.push_back(theta);
    m_vec_phi.push_back(phi);
    m_vec_vertexR.push_back(vertexR);
    m_vec_mchits.push_back(MCParticleHits[associatedParticle].size());
    m_vec_trkhits.push_back(nHits);
    m_vec_purity.push_back(purity);
    m_vec_misassocHits.push_back(countMisassocHits);
    m_vec_misassocFirstHitPos.push_back(misassocFirstHitPos);
    m_vec_missedHits.push_back(countMissedHits);
    m_vec_missedFirstHitPos.push_back(missedFirstHitPos);

  } // end loop tracks
  
  _trackingTree->Fill();

  _eventNumber++;
  
  //Clear the trees
  clearTreeVar();
  

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

int ClicTrackingValidate::getSubdetector(TrackerHit* hit, UTIL::BitField64 &encoder){
  const int celId = hit->getCellID0() ;
  encoder.setValue(celId) ;
  int subdet = encoder[lcio::LCTrackerCellID::subdet()];
  return subdet;
}

int ClicTrackingValidate::getSubdetectorFromSim(LCCollection* collection, UTIL::BitField64 &encoder){
  SimTrackerHit* hit = static_cast<SimTrackerHit*>(collection->getElementAt(0));
  const int celId = hit->getCellID0();
  encoder.setValue(celId);
  int subdet = encoder[lcio::LCTrackerCellID::subdet()];
  return subdet;
}

void ClicTrackingValidate::clearTreeVar(){

  m_vec_pdg.clear();
  m_vec_pt.clear();
  m_vec_theta.clear();
  m_vec_phi.clear();
  m_vec_vertexR.clear();
  m_vec_mchits.clear();
  m_vec_trkhits.clear();
  m_vec_purity.clear();
  m_vec_misassocHits.clear();
  m_vec_misassocFirstHitPos.clear();
  m_vec_missedHits.clear();
  m_vec_missedFirstHitPos.clear();
  
}
