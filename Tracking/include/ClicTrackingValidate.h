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

  // Call to get collections
  void getCollection(LCCollection*&, std::string, LCEvent*);
  
  // Get the subdetector ID from a collection 
  int getSubdetector(LCCollection*, UTIL::BitField64&);

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

};

#endif
