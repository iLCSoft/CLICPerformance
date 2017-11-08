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
#include <EVENT/MCParticle.h>
#include "MarlinTrk/IMarlinTrkSystem.h"
#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <AIDA/AIDA.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>


class TTree;
/* class TFile; */

using namespace lcio ;
using namespace marlin ;
using namespace AIDA ;

class ClicEfficiencyCalculator : public Processor {
		
public:
	
	virtual Processor*  newProcessor() { return new ClicEfficiencyCalculator ; }
	
	ClicEfficiencyCalculator() ;
	ClicEfficiencyCalculator(const ClicEfficiencyCalculator&) = delete;
	ClicEfficiencyCalculator& operator=(const ClicEfficiencyCalculator&) = delete;

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
	int getSubdetector(TrackerHit*, UTIL::BitField64&);
	int getLayer(TrackerHit*, UTIL::BitField64&);
  int getUniqueHits(std::vector<TrackerHit*>, UTIL::BitField64&);

	bool isReconstructable(MCParticle*& particle, std::string cut, UTIL::BitField64&);

	void clearTreeVar();

	
protected:
	
  // Collection names for (in/out)put
  std::vector<std::string> m_inputTrackerHitCollections = {};
  std::vector<std::string> m_inputTrackerHitRelationCollections = {};
  std::string m_inputParticleCollection = "";
  std::string m_inputPhysicsParticleCollection = "";
  std::string m_inputTrackCollection = "";
  std::string m_inputTrackRelationCollection = "";
  std::string m_outputEfficientMCParticleCollection = "";
  std::string m_outputInefficientMCParticleCollection = "";
	
  // Run and event counters
  int m_eventNumber = 0;
  int m_runNumber = 0;
	
  // Track fit factory
  MarlinTrk::IMarlinTrkSystem* trackFactory = NULL;
	
  // Parameters
  bool m_fullOutput = false;
  bool m_simpleOutput = false;
  double m_magneticField = 0.0;
  std::string m_cuts = "";
  std::string m_mcTreeName = "";
  std::string m_purityTreeName = "";
  std::string m_efficiencyTreeName = "";
  std::map<std::string,double> m_particles = {};
  std::map<std::string,double> m_reconstructedParticles = {};
  std::map<MCParticle*, std::vector<TrackerHit*> > particleHits = {};

  // Histograms
  TH2F* m_thetaPtMCParticle = NULL;
  TH2F* m_thetaPtMCParticleReconstructed = NULL;
  TH2F* m_thetaPtMCParticleEfficiency = NULL;

  // Trees
  TTree *m_mctree = NULL;
  std::vector<int > m_mcCat = {};
  std::vector<double > m_mcTheta = {};
  std::vector<double > m_mcPt = {};
  std::vector<int > m_mcIsDecayedInTracker = {};
  std::vector<int > m_mcNHitsTot = {};
  std::vector<int > m_mcNHitsVXD = {};
  std::vector<int > m_mcNTracks = {};
  std::vector<int > m_mcNTrkHits = {};
  std::vector<int > m_mcThetaTrk = {};
  std::vector<int > m_mcPtTrk = {};
  std::vector<int > m_mcPhiTrk = {};
  std::vector<int > m_mcNTracksCone = {};

  TTree *m_trackTree = NULL;
  std::vector<double > m_vec_vx_reconstructable = {};
  std::vector<double > m_vec_vy_reconstructable = {};
  std::vector<double > m_vec_vz_reconstructable = {};
  std::vector<double > m_vec_vr_reconstructable = {};
  std::vector<double > m_vec_pt_reconstructable = {};
  std::vector<double > m_vec_theta_reconstructable = {};
  std::vector<bool > m_vec_is_reconstructed = {};

  TTree *m_purityTree = NULL;
  std::vector<int > m_vec_nhits_vtx = {};
  std::vector<int > m_vec_nhits_trk = {};
  std::vector<int > m_vec_nhits = {};
  std::vector<double > m_vec_purity = {};
  std::vector<int > m_vec_pdg = {};
  std::vector<double > m_vec_theta = {};
  std::vector<double > m_vec_phi = {};
  std::vector<double > m_vec_p = {};
  
  TTree *m_simplifiedTree = NULL;
  double m_type = 0.0, m_pt = 0.0, m_theta = 0.0, m_phi = 0.0, m_vertexR = 0.0, m_distClosestMCPart = 0.0,  m_closeTracks = 0.0, m_purity = 0.0;
  int m_nHits = 0, m_nHitsMC = 0;
  bool m_reconstructed = false, m_signal = false;


} ;

#endif



