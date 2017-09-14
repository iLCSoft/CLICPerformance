#ifndef TrackChecker_h
#define TrackChecker_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <string>
#include <vector>
#include <map>

#include <EVENT/LCCollection.h>
#include <EVENT/Track.h>

#include <AIDA/AIDA.h>
#include "TH1F.h"

class TTree;
class TCanvas;

using namespace lcio ;
using namespace marlin ;
using namespace AIDA ;

class TrackChecker : public Processor {
		
 public:
	
  virtual Processor*  newProcessor() { return new TrackChecker ; }
	
  TrackChecker() ;

  TrackChecker(const TrackChecker&) = delete;
  TrackChecker& operator=(const TrackChecker&) = delete;

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
	
  void angleInFixedRange(double& angle);

  void clearEventVariables();

  double getDistanceFromClosestHit(Track*& trackA, Track*& trackB);

  double getDist(TrackerHit*& hitA, TrackerHit*& hitB);



 protected:
	
  // Collection names for (in/out)put
  std::string m_inputTrackCollection = "";
  std::string m_inputTrackRelationCollection = "";
  std::string m_inputParticleCollection = "" ;
 	
  // Run and event counters
  int m_eventNumber = 0;
  int m_runNumber = 0;
  
  // Magnetic field
  float m_magneticField = 0;
	
  //Flags
  bool m_useOnlyTree = false;

  // Output histograms

  TCanvas* pulls = NULL;
  TCanvas* res = NULL;

  TH1F* m_omegaPull = NULL;
  TH1F* m_phiPull = NULL;
  TH1F* m_tanLambdaPull = NULL;
  TH1F* m_d0Pull = NULL;
  TH1F* m_z0Pull = NULL;

  TH1F* m_omegaResidual = NULL;
  TH1F* m_phiResidual = NULL;
  TH1F* m_tanLambdaResidual = NULL;
  TH1F* m_d0Residual = NULL;
  TH1F* m_z0Residual = NULL;
  
  TH1F* m_omegaMCParticle = NULL;
  TH1F* m_phiMCParticle = NULL;
  TH1F* m_tanLambdaMCParticle = NULL;
  TH1F* m_d0MCParticle = NULL;
  TH1F* m_z0MCParticle = NULL;

  TH1F* m_omegaTrack = NULL;
  TH1F* m_phiTrack = NULL;
  TH1F* m_tanLambdaTrack = NULL;
  TH1F* m_d0Track = NULL;
  TH1F* m_z0Track = NULL;
  
  TH1F* m_trackChi2 = NULL;
  

  std::string m_treeName = "";

  TTree *perftree = NULL;

  std::vector<double > truePt = {};
  std::vector<double > trueTheta = {};
  std::vector<double > truePhi = {};
  std::vector<double > trueD0 = {};
  std::vector<double > trueZ0 = {};
  std::vector<double > trueP = {};
  std::vector<int >    trueID = {};

  std::vector<double > recoPt = {};
  std::vector<double > recoTheta = {};
  std::vector<double > recoPhi = {};
  std::vector<double > recoD0 = {};
  std::vector<double > recoZ0 = {};
  std::vector<double > recoP = {};

  std::vector<int > recoNhits = {};
  std::vector<double > recoChi2OverNDF = {};
  std::vector<double > recoMinDist = {};

  std::vector<double > pullOmega = {};
  std::vector<double > pullPhi = {};
  std::vector<double > pullTanLambda = {};
  std::vector<double > pullD0 = {};
  std::vector<double > pullZ0 = {};

  std::vector<double > resOmega = {};
  std::vector<double > resPhi = {};
  std::vector<double > resTanLambda = {};
  std::vector<double > resD0 = {};
  std::vector<double > resZ0 = {};



} ;

#endif



