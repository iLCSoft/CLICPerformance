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
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>


class TTree;

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
	int getSubdetector(TrackerHit*, UTIL::BitField64&);

	bool isReconstructable(MCParticle*& particle, std::string cut);

	void clearTreeVar();

	
protected:
	
	// Collection names for (in/out)put
	std::vector<std::string> m_inputTrackerHitCollections ;
	std::vector<std::string> m_inputTrackerHitRelationCollections ;
	std::string m_inputParticleCollection ;
	std::string m_inputTrackCollection ;
	std::string m_inputTrackRelationCollection;
	std::string m_notRecoMCColName;
	
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

	int m_vertexBarrelID;
	std::string m_cuts;
	std::map<MCParticle*, std::vector<TrackerHit*> > particleHits;
	bool m_morePlots;

	// Plots 

	TCanvas *eff_vs_theta;
	TGraphAsymmErrors *g_eff_vs_theta;
	TH1F *h_theta_reconstructed ;
	TH1F *h_theta_reconstructable ;

	TCanvas *eff_vs_pt;
	TGraphAsymmErrors *g_eff_vs_pt;
	TH1F *h_pt_reconstructed ;
	TH1F *h_pt_reconstructable ;

	// Tree

	TTree *mctree ;
	std::vector<int > m_mcCat;
	std::vector<double > m_mcTheta;
	std::vector<double > m_mcPt;
	std::vector<int > m_mcIsDecayedInTracker;
	//std::vector<std::vector<int > > m_mcNHits;
	//std::vector<int > m_mcNHits_helper;
	std::vector<int > m_mcNHitsTot;
	std::vector<int > m_mcNHitsVXD;

	std::vector<int > m_mcNTracks; 
	std::vector<int > m_mcNTrkHits; 
	std::vector<int > m_mcThetaTrk; 
	std::vector<int > m_mcPtTrk; 
	std::vector<int > m_mcPhiTrk; 
	std::vector<int > m_mcNTracksCone; 



} ;

#endif



