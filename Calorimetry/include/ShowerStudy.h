#ifndef ShowerStudy_h
#define ShowerStudy_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>


#include "TH2F.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TFile.h"
#include <map>

using namespace lcio ;
using namespace marlin ;

class TTree;


class ShowerStudy : public Processor {
                
public:
        
    virtual Processor*  newProcessor() { return new ShowerStudy ; }
    
    ShowerStudy() ;
    
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
        
protected:
        
    // Collection names for (in/out)put
    std::string m_inputCalorimeterHitCollection ;
    std::string m_inputMCParticleCollection ;

    std::string m_rootFileName;
    
    // Run and event counters
    int m_eventNumber ;
    int m_runNumber ;
    
    // Call to get collections
    void getCollection(LCCollection*&, std::string, LCEvent*);
    
    // Plots 
    TProfile2D * m_showerProfileLayers;
    TProfile2D * m_showerProfile;
    TH2F * m_totalEnergyHist;
    TTree * m_outputTree;
    TFile * m_rootFile;
    //Tree branch variables
    
    float m_trueEnergy;
    float m_totalEnergy;
    unsigned int m_nhits;
    std::vector<float> * m_hitEnergies;
    std::vector<float> * m_hitDistances;
    std::vector<int> * m_hitLayers;
    std::vector<float> * m_hitLayerThicknesses;
    std::vector<float> * m_hitLayerRadiationLengths;


    
    

} ;

#endif



