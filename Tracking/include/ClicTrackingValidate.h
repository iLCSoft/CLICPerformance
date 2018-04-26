#ifndef ClicTrackingValidate_h
#define ClicTrackingValidate_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <string>
#include <vector>
#include <map>

#include "DDRec/Surface.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <EVENT/SimTrackerHit.h>
#include <UTIL/CellIDDecoder.h>
#include <gsl/gsl_rng.h>

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

class ClicTrackingValidate : public Processor {

 public:
  virtual Processor* newProcessor() {return new ClicTrackingValidate;}
  
  ClicTrackingValidate();
  ClicTrackingValidate(const ClicTrackingValidate&) = delete;
  ClicTrackingValidate& operator=(const ClicTrackingValidate&) = delete;

  // Initialisation - run at the beginning to start histograms, etc.
  virtual void init();

  // Called at the beginning of every run
  virtual void processRunHeader( LCRunHeader* run );

  // Run over each event - the main algorithm
  virtual void processEvent( LCEvent* evt );
  
  // Run at the end of each event
  virtual void check( LCEvent* evt );

  // Called at the very end for cleanup, histogram saving, etc.
  virtual void end();

  // Clear trees
  virtual void clearTreeVar();

  // Call to get collections
  void getCollection(LCCollection*&, std::string, LCEvent*);
  
  // Get the subdetector ID from a collection 
  int getSubdetector(LCCollection*, UTIL::BitField64&);
  int getSubdetector(TrackerHit*, UTIL::BitField64&);
  int getSubdetectorFromSim(LCCollection*, UTIL::BitField64&);

 protected:
  
  // Collection names for input
  std::string _inputParticleCollection = "";
  std::string _inputTrackCollection = "";
  std::vector<std::string> _inputSimTrackerHitCollections = {};
  std::vector<std::string> _inputTrackerHitCollections = {};
  std::vector<std::string> _inputTrackerHitRelationCollections = {};

  // Run and event counters
  int _eventNumber = 0;
  int _runNumber = 0;

  // Parameters
  double _magneticField = 0.0;
  std::string _trackingTreeName = "";

  // Trees
  TTree *_trackingTree = NULL;
  std::vector<int> m_vec_pdg = {};
  std::vector<double> m_vec_pt = {};
  std::vector<double> m_vec_theta = {};
  std::vector<double> m_vec_phi = {};
  std::vector<double> m_vec_vertexR = {};
  std::vector<int> m_vec_mchits = {};
  std::vector<int> m_vec_trkhits = {};
  std::vector<double> m_vec_purity = {};
  std::vector<int> m_vec_misassocHits = {};
  std::vector<double> m_vec_misassocFirstHitPos = {};
  std::vector<int> m_vec_missedHits = {};
  std::vector<double> m_vec_missedFirstHitPos = {};

};

// Compare content of two std::vector - checks if subset or equal
bool isSubsetOrEqual(std::vector<SimTrackerHit*> const& a,std::vector<SimTrackerHit*> const& b){

  for(auto const& x:a){
    if(std::find(b.begin(),b.end(),x)==b.end()) return false;

  }
  return true;
 
}
#endif
