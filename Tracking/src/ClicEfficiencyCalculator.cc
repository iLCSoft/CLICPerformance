/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "ClicEfficiencyCalculator.h"
#include "DDRec/API/IDDecoder.h"

#include "MarlinTrk/MarlinTrkUtils.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/HelixFit.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/Factory.h"
#include "MarlinTrk/MarlinTrkDiagnostics.h"
#include "MarlinTrk/IMarlinTrkSystem.h"

#include <UTIL/CellIDDecoder.h>
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

using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace DD4hep ;
using namespace AIDA ;

ClicEfficiencyCalculator aClicEfficiencyCalculator;

/*

 Efficiency calculator for tracking. WARNING: Only works with 1 input collection at present
 
*/

ClicEfficiencyCalculator::ClicEfficiencyCalculator() : Processor("ClicEfficiencyCalculator") {
	
	// modify processor description
	_description = "ClicEfficiencyCalculator calculates the tracking efficiency and makes some performance plots" ;

	// Input collections
	registerInputCollection( LCIO::TRACK,
													"TrackCollectionName",
													"Track collection name",
													m_inputTrackCollection,
													std::string("CATracks"));
	
//	registerInputCollection( LCIO::LCRELATION,
//													"TrackRelationCollectionName",
//													"Track relation collection name",
//													m_inputTrackRelationCollection,
//													std::string("SiTrackRelations"));
	
	registerInputCollection( LCIO::MCPARTICLE,
													"MCParticleCollectionName",
													"Name of the MCParticle input collection",
													m_inputParticleCollection,
													std::string("MCParticle"));
	
	// All tracker hit collections
	std::vector<std::string> inputTrackerHitCollections;
	inputTrackerHitCollections.push_back(std::string("VXDTrackerHits"));
	
	registerInputCollections( LCIO::TRACKERHITPLANE,
													 "TrackerHitCollectionNames" ,
													 "Name of the TrackerHit input collections"  ,
													 m_inputTrackerHitCollections ,
													 inputTrackerHitCollections ) ;
	
	// All tracker hit relation collections
	std::vector<std::string> inputTrackerHitRelationCollections;
	inputTrackerHitRelationCollections.push_back(std::string("VXDTrackerHitRelations"));
	
	registerInputCollections( LCIO::LCRELATION,
													 "TrackerHitRelCollectionNames",
													 "Name of TrackerHit relation collections",
													 m_inputTrackerHitRelationCollections,
													 inputTrackerHitRelationCollections );
	
}


void ClicEfficiencyCalculator::init() {
	
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
	
	// Define the purity cut - should be made configurable
	m_purity = 0.5;
	
	// Register this process
	Global::EVENTSEEDER->registerProcessor(this);
	
}


void ClicEfficiencyCalculator::processRunHeader( LCRunHeader* run) {
	++m_runNumber ;
}

void ClicEfficiencyCalculator::processEvent( LCEvent* evt ) {
	
	std::cout<<"Processing event "<<m_eventNumber<<std::endl;
	// First pick up all of the collections that will be used - tracks, MCparticles, hits from relevent subdetectors - and their relations
	
	// Initialise CELLID encoder to get the subdetector ID from a hit
	UTIL::BitField64 m_encoder( lcio::ILDCellID0::encoder_string ) ;
	std::vector<int> m_collections;

	// Get the collection of tracks
	LCCollection* trackCollection = 0 ;
	getCollection(trackCollection, m_inputTrackCollection, evt); if(trackCollection == 0) return;
	
	// Get the collection of MC particles
	LCCollection* particleCollection = 0 ;
	getCollection(particleCollection, m_inputParticleCollection, evt); if(particleCollection == 0) return;
	
	// Get the collection of track relations to MC particles
//	LCCollection* trackRelationCollection = 0 ;
//	getCollection(trackRelationCollection, m_inputTrackRelationCollection, evt); if(trackRelationCollection == 0) return;
	
	// Create the relations navigator
//	LCRelationNavigator* trackRelation = new LCRelationNavigator( trackRelationCollection );
	
	// Make objects to hold all of the tracker hit and relation collections
	std::map<int,LCCollection*> trackerHitCollections;
	std::map<int,LCCollection*> trackerHitRelationCollections;
	std::map<int,LCRelationNavigator*> relations;
	std::map<int,bool> activeDetectors;
	
	// Loop over each input collection and get the data
	for(unsigned int collection=0; collection<m_inputTrackerHitCollections.size();collection++){
		
		// Get the collection of tracker hits
		LCCollection* trackerHitCollection = 0 ;
		getCollection(trackerHitCollection, m_inputTrackerHitCollections[collection], evt); if(trackerHitCollection == 0) continue;
		
		// If there are no hits in the hit collection, continue
		if(trackerHitCollection->getNumberOfElements() == 0) continue;

		// Get the collection of tracker hit relations
		LCCollection* trackerHitRelationCollection = 0 ;
		getCollection(trackerHitRelationCollection, m_inputTrackerHitRelationCollections[collection], evt);
  
		// Create the relations navigator
		LCRelationNavigator* relation = new LCRelationNavigator( trackerHitRelationCollection );
		
		// Get the subdetector for this set of collections
		int subdetector = getSubdetector(trackerHitCollection,m_encoder);
		activeDetectors[subdetector] = true;

		// Save all of the created collections
		trackerHitCollections[subdetector] = trackerHitCollection;
		trackerHitRelationCollections[subdetector] = trackerHitRelationCollection;
		relations[subdetector] = relation;
		m_collections.push_back(subdetector);
  
	}

	std::map<MCParticle*,int> particleTracks;
	
	/*
		Look at all tracks that were reconstructed and calculate their purity. This can be used to point to all
	 	reconstructed MC particles. Then loop over MC particles and for each that was not reconstructed check the
	  particle properties. Use these to define efficiency etc. For this we also need the list of hits belonging 
	  to each particle.
  */
	
	// Loop over all tracks
	int nTracks = trackCollection->getNumberOfElements();
	for(int itTrack=0;itTrack<nTracks;itTrack++){
		
		// Get the track
		Track* track = dynamic_cast<Track*>( trackCollection->getElementAt(itTrack) ) ;
		
		// Get the hits
		const TrackerHitVec& hitVector = track->getTrackerHits();
		
		// Some storage to keep track of which particles have hits on this track
		std::vector<MCParticle*> trackParticles;
		std::map<MCParticle*,int> trackParticleHits;
		
		// Loop over hits and check the MC particle that they correspond to, and which subdetector they are on
		int nHits = hitVector.size();
		int nExcluded = 0;
		for(int itHit=0;itHit<nHits;itHit++){
			// Get the tracker hit
			TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(hitVector.at(itHit));
			// Get the subdetector number
			int subdetector = getSubdetector(hit,m_encoder);
			// Check if this subdetector is to be included in the efficiency calculation
			if( activeDetectors.count(subdetector) == 0 ){nExcluded++; continue;}
			// Get the simulated hit
			const LCObjectVec& simHitVector = relations[subdetector]->getRelatedToObjects( hit );
			// Take the first hit only (this should be changed? Yes - loop over all related simHits and add an entry for each mcparticle so that this hit is in each fit)
			SimTrackerHit* simHit = dynamic_cast<SimTrackerHit*>(simHitVector.at(0));
			// Get the particle belonging to that hit
			MCParticle* particle = simHit->getMCParticle();
			// Check if this particle already has hits on this track
			if( std::find(trackParticles.begin(), trackParticles.end(), particle) == trackParticles.end() ) trackParticles.push_back(particle);
			// Increment the number of counts
			trackParticleHits[particle]++;
		}
		
		// Check the purity of the track, treat ghosts
		int nTrackParticles = trackParticles.size(); int maxHits=0;
		MCParticle* associatedParticle = 0;
		for(int itParticle=0; itParticle<nTrackParticles; itParticle++){
			int nParticleHits = trackParticleHits[trackParticles[itParticle]];
			if(nParticleHits>maxHits){ maxHits=nParticleHits; associatedParticle = trackParticles[itParticle]; }
		}
		double purity = (double)maxHits/(double)(nHits-nExcluded);
		if(purity < m_purity) continue; // do something with ghosts here

		// Now have a track which is associated to a particle
		particleTracks[associatedParticle]++;

	}
	
	
	/*
	 Now for each MC particle we want the list of hits belonging to it. The most
	 efficient way is to loop over all hits once, and store the pointers in a
	 map, with the key a pointer to the MC particle. We can then loop over each
	 MC particle at the end and get all of the hits, before making a track.
  */
	
	// Make the container
	std::map<MCParticle*, std::vector<TrackerHit*> > particleHits;
	
	// Loop over all input collections
	for(unsigned int itCollection=0; itCollection<m_collections.size();itCollection++){
	
		// Get the collection ID (subdetector number)
		int collection = m_collections[itCollection];

		// Loop over tracker hits
		int nHits = trackerHitCollections[collection]->getNumberOfElements();
		for(int itHit=0;itHit<nHits;itHit++){
			
			// Get the hit
			TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>( trackerHitCollections[collection]->getElementAt(itHit) ) ;
			
			// Get the related simulated hit(s)
			const LCObjectVec& simHitVector = relations[collection]->getRelatedToObjects( hit );
			
			// Take the first hit only (this should be changed? Yes - loop over all related simHits and add an entry for each mcparticle so that this hit is in each fit)
			SimTrackerHit* simHit = dynamic_cast<SimTrackerHit*>(simHitVector.at(0));
			
			// Get the particle belonging to that hit
			MCParticle* particle = simHit->getMCParticle();
			
			// Push back the element into the container
			particleHits[particle].push_back(hit);
			
		}
	}
	
	/*
	   Now there is a map containing the number of tracks associated to each MC particle, and a 
	   map of all hits on each particle. Loop over all MC particles, decide if they are 
	   reconstructable or not, and check if they were reconstructed more than once (clones).
	 */
	
	int nReconstructed=0;
	int nReconstructable=0;
	// Loop over particles
	int nParticles = particleCollection->getNumberOfElements();
	for(int itParticle=0;itParticle<nParticles;itParticle++){
		
		// Get the particle
		MCParticle* particle = dynamic_cast<MCParticle*>( particleCollection->getElementAt(itParticle) ) ;
		
		// Check if it was reconstructed
		if(particleTracks.count(particle)){
			// Assumption: if particle was reconstructed then it is reconstructable!
			m_particles["all"]++; m_particles["reconstructable"]++; nReconstructable++;
			m_reconstructedParticles["all"]++; nReconstructed++;
			// Check if clones were produced (1 particle, more than 1 track)
			if(particleTracks[particle] > 1) m_reconstructedParticles["clones"]+=(particleTracks[particle]-1);
			continue;
		}
		
		// Now decide on criteria for different particle types/classifications
		m_particles["all"]++; // all particles
		
		// No hits in the input collections
		if(particleHits.count(particle) == 0) continue; // No hits in the input collections

		// Exclude particles which are not "real"
		if(particle->vertexIsNotEndpointOfParent()) continue;
		
		// Only make tracks with 3 or more hits
		std::vector<TrackerHit*> trackHits = particleHits[particle];
		if(trackHits.size() < 3) continue;

		m_particles["reconstructable"]++; // reconstructable particles
		nReconstructable++;
		
	}
	
	// Increment the event number
	m_eventNumber++ ;
	std::cout<<"For this event reconstructed "<<100.*(double)nReconstructed/(double)nReconstructable<<" % ("<<nReconstructed<<"/"<<nReconstructable<<")"<<std::endl;
	
}

void ClicEfficiencyCalculator::check( LCEvent * evt ) {
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ClicEfficiencyCalculator::end(){

	streamlog_out(MESSAGE) << " end()  " << name()
	<< " processed " << m_eventNumber << " events in " << m_runNumber << " runs "
	<< std::endl ;
	
	// Calculate efficiency results
	streamlog_out(MESSAGE)<<std::fixed<<std::setprecision(2)<<"Reconstructable particle efficiency: "<<100.*m_reconstructedParticles["all"]/m_particles["reconstructable"]<<" % ("<<std::setprecision(0)<<m_reconstructedParticles["all"]<<"/"<<m_particles["reconstructable"]<<")"<<std::endl;
	
}

void ClicEfficiencyCalculator::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
  try{
    collection = evt->getCollection( collectionName ) ;
  }
  catch(DataNotAvailableException &e){
		streamlog_out( DEBUG5 )<<"- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

int ClicEfficiencyCalculator::getSubdetector(LCCollection* collection, UTIL::BitField64 &encoder){
	TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>( collection->getElementAt(0) ) ;
	const int celId = hit->getCellID0() ;
	encoder.setValue(celId) ;
	int subdet = encoder[lcio::ILDCellID0::subdet];
	return subdet;
}

int ClicEfficiencyCalculator::getSubdetector(TrackerHitPlane* hit, UTIL::BitField64 &encoder){
	const int celId = hit->getCellID0() ;
	encoder.setValue(celId) ;
	int subdet = encoder[lcio::ILDCellID0::subdet];
	return subdet;
}


  
  
  
  
  
  
