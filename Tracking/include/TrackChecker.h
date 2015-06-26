#ifndef TrackChecker_h
#define TrackChecker_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <string>
#include <vector>
#include <map>

#include <EVENT/LCCollection.h>

#include <AIDA/AIDA.h>
#include "TH1F.h"

using namespace lcio ;
using namespace marlin ;
using namespace AIDA ;

class TrackChecker : public Processor {
		
public:
	
	virtual Processor*  newProcessor() { return new TrackChecker ; }
	
	TrackChecker() ;
	
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

	
protected:
	
  // Collection names for (in/out)put
  std::string m_inputTrackCollection ;
  std::string m_inputTrackRelationCollection ;
  std::string m_inputParticleCollection ;
 	
  // Run and event counters
	int m_eventNumber ;
	int m_runNumber ;
  
  // Magnetic field
  float m_magneticField;
	
  // Output histograms
  TH1F* m_omegaPull;
  TH1F* m_phiPull;
  TH1F* m_tanLambdaPull;
  TH1F* m_d0Pull;
  TH1F* m_z0Pull;

  TH1F* m_omegaResidual;
  TH1F* m_phiResidual;
  TH1F* m_tanLambdaResidual;
  TH1F* m_d0Residual;
  TH1F* m_z0Residual;
  
  TH1F* m_omegaMCParticle;
  TH1F* m_phiMCParticle;
  TH1F* m_tanLambdaMCParticle;
  TH1F* m_d0MCParticle;
  TH1F* m_z0MCParticle;

  TH1F* m_omegaTrack;
  TH1F* m_phiTrack;
  TH1F* m_tanLambdaTrack;
  TH1F* m_d0Track;
  TH1F* m_z0Track;
  
  TH1F* m_trackChi2;
  
  // Output file name
  std::string m_outputFile;
} ;

#endif



