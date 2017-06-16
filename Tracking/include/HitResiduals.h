#ifndef HitResiduals_h
#define HitResiduals_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include "MarlinTrk/Factory.h"
#include "marlin/Global.h"
#include <EVENT/Track.h>
#include "DDRec/Surface.h"
#include "DDRec/DetectorSurfaces.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/SurfaceHelper.h"

#include <string>


using namespace lcio ;
using namespace marlin ;
namespace MarlinTrk{
  class IMarlinTrkSystem ;
}
using namespace MarlinTrk ;
using namespace DD4hep::DDRec;



class TFile;
class TTree;

class HitResiduals : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new HitResiduals ; }
  
  
  HitResiduals() ;

  HitResiduals(const HitResiduals&) = delete;
  HitResiduals& operator=(const HitResiduals&) = delete;
  
  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;
  
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;
  
  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 
  
  
  virtual void check( LCEvent * evt ) ; 
  
  
  /** Called after data processing for clean up.
   */
  virtual void end() ;
  

  LCCollection* GetCollection( LCEvent*& evt, std::string colName );
  int FitInit2( Track*& track, MarlinTrk::IMarlinTrack*& _marlinTrk );
  int FitInitFromLCIOTrackState( Track*& track, MarlinTrk::IMarlinTrack*& _marlinTrk );

  
 protected:

  /** Input collection name.
   */
  std::string _inputTrkColName = "";

  std::string _outFileName = "";
  std::string _treeName = "";

  TFile* _out = NULL;
  TTree* _tree = NULL;

  int _nRun = 0;
  int _nEvt = 0;

  std::vector<double > _resU = {};
  std::vector<double > _resV = {};
  std::vector<int > _subdet = {};
  std::vector<int > _layer = {};
  std::vector<int > _side = {};

  MarlinTrk::IMarlinTrkSystem* _trksystem = NULL;
  SurfaceMap _surfMap = {};
  bool _MSOn = true;
  bool _ElossOn = true;
  bool _SmoothOn = false;
  double _Max_Chi2_Incr = 0.0;
  double _bField = 0.0;

} ;

#endif



