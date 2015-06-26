#ifndef ClicEfficiencyCalculator_h
#define ClicEfficiencyCalculator_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <string>
#include <vector>
#include <map>

#include <gsl/gsl_rng.h>
#include "DDRec/Surface.h"
#include <EVENT/LCCollection.h>
#include "MarlinTrk/IMarlinTrkSystem.h"
#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <AIDA/AIDA.h>

using namespace lcio ;
using namespace marlin ;
using namespace AIDA ;

class ClicEfficiencyCalculator : public Processor {
		
public:
	
	virtual Processor*  newProcessor() { return new ClicEfficiencyCalculator ; }
	
	ClicEfficiencyCalculator() ;
	
	// Initialisation - run at the beginning to start histograms, etc.
	virtual void init() ;
	
	// Called at the beginning of every run
	virtual void processRunHeader( LCRunHeader* run ) ;
	
	// Run over each event - the main algorithm
	virtual void processEvent( LCEvent * evt ) ;
	
  // Run at the end of each event
	virtual void check( LCEvent * evt ) ;
	
	// Called at the very end for cleanup, histogram saving, etc.
	virtual void end() ;
	
  // Call to get collections
  void getCollection(LCCollection*&, std::string, LCEvent*);

	// Get the subdetector ID from a collection or hit
	int getSubdetector(LCCollection*, UTIL::BitField64&);
	int getSubdetector(TrackerHitPlane*, UTIL::BitField64&);

	
protected:
	
	// Collection names for (in/out)put
	std::vector<std::string> m_inputTrackerHitCollections ;
	std::vector<std::string> m_inputTrackerHitRelationCollections ;
	std::string m_inputParticleCollection ;
	std::string m_inputTrackCollection ;
	std::string m_inputTrackRelationCollection;
	
	// Run and event counters
	int m_eventNumber ;
	int m_runNumber ;
	
	// Track fit factory
	MarlinTrk::IMarlinTrkSystem* trackFactory;
	
	// Parameters
	double m_purity;
	double m_magneticField;
	std::map<std::string,double> m_particles;
	std::map<std::string,double> m_reconstructedParticles;
	
} ;

#endif



