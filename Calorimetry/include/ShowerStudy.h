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
    std::string m_inputLeakageCalorimeterHitCollection;
    std::string m_rootFileName;
    
    // Run and event counters
    int m_eventNumber ;
    int m_runNumber ;
    
    // Call to get collections
    void getCollection(LCCollection*&, std::string, LCEvent*);
    
    // Plots 
    TH2F * m_showerHistLayers;
    TH2F * m_showerHist;
    TH2F * m_showerHistX0;
    TProfile * m_leakageProfile;
    

    TH2F * m_raw_showerHistLayers;
    TH2F * m_raw_showerHist;
    TH2F * m_raw_showerHistX0;

    TH2F * m_totalEnergyHist;
    TTree * m_outputTree;
    TFile * m_rootFile;
    //Tree branch variables
    
    float m_trueEnergy;
    float m_totalEnergy;
    float m_totalLeakEnergy;

    unsigned int m_nhits;
    unsigned int m_leak_nhits;

    
    std::vector<float> * m_hitEnergies;
    std::vector<float> * m_raw_hitEnergies;

    std::vector<float> * m_hit_x;
    std::vector<float> * m_hit_y;
    std::vector<float> * m_hit_z;

    std::vector<int> * m_hitLayers;
    std::vector<float> * m_hitLayerThicknesses;
    std::vector<float> * m_hitLayerRadiationLengths;
    std::vector<float> * m_hitLayerIntRadiationLengths;
    std::vector<float> * m_hitLayerDistances;


    std::vector<float> * m_leak_hitEnergies;
    std::vector<float> * m_leak_hit_x;
    std::vector<float> * m_leak_hit_y;
    std::vector<float> * m_leak_hit_z;

    std::vector<int> * m_leak_hitLayers;
    std::vector<float> * m_leak_hitLayerThicknesses;
    std::vector<float> * m_leak_hitLayerRadiationLengths;
    std::vector<float> * m_leak_hitLayerIntRadiationLengths;
    std::vector<float> * m_leak_hitLayerDistances;
    std::vector<float> * m_leak_raw_hitEnergies;


} ;

#endif



