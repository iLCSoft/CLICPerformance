#ifndef SkimSecondaryVertices_h
#define SkimSecondaryVertices_h 1

#include "marlin/Processor.h"
#include "lcio.h"

#include <string>
#include <vector>
#include <map>

#include "DDRec/Surface.h"
#include "MarlinTrk/IMarlinTrkSystem.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Vertex.h>
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

class SkimSecondaryVertices : public Processor {
 public:
       
        virtual Processor*  newProcessor() { return new SkimSecondaryVertices ; }
       
        SkimSecondaryVertices() ;
        SkimSecondaryVertices(const SkimSecondaryVertices&) = delete;
        SkimSecondaryVertices& operator=(const SkimSecondaryVertices&) = delete;

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

        // Clear the tree
        void clearTreeVar();

 protected:
        std::string _inputBuildUpVertices = "";
        std::string _outputBuildUpVerticesSkimmed = "";
        std::string _vertexTreeName = "";
        std::string _vertexSkimmedTreeName = "";
        int _eventNumber = 0;
        int _runNumber = 0;

        // Geometry constraints
        std::vector<double > _barrel_radii = {};
        double _barrel_halfLength = 0.0;
        std::vector<double > _endcap_z = {};

        // Trees
        TTree *_vertexTree = NULL;
        int _vtxN = 0;
        std::vector<double > _vtxX = {};
        std::vector<double > _vtxY = {};
        std::vector<double > _vtxR = {};
        std::vector<double > _vtxZ = {};
        TTree *_vertexSkimmedTree = NULL;
        int _vtxskimN = 0;
        std::vector<double > _vtxskimR = {};
        std::vector<double > _vtxskimZ = {};

};

#endif
