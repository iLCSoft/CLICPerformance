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

    ShowerStudy(const ShowerStudy&) = delete;
    ShowerStudy& operator=(const ShowerStudy&) = delete;

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
        
    
  // Call to get collections
  void getCollection(LCCollection*&, std::string, LCEvent*);

  // Collection names for (in/out)put
  std::string m_inputCalorimeterHitCollection = "";
  std::string m_inputMCParticleCollection = "";
  std::string m_inputLeakageCalorimeterHitCollection = "";
  std::string m_rootFileName = "";
    
  // Run and event counters
  int m_eventNumber = 0;
  int m_runNumber = 0;
    
  // Plots
  TH2F * m_showerHistLayers = NULL;
  TH2F * m_showerHist = NULL;
  TH2F * m_showerHistX0 = NULL;
  TProfile * m_leakageProfile = NULL;
    

  TH2F * m_raw_showerHistLayers = NULL;
  TH2F * m_raw_showerHist = NULL;
  TH2F * m_raw_showerHistX0 = NULL;

  TH2F * m_totalEnergyHist = NULL;
  TTree * m_outputTree = NULL;
  TFile * m_rootFile = NULL;
  //Tree branch variables
    
  float m_trueEnergy = 0.0;
  float m_totalEnergy = 0.0;
  float m_totalLeakEnergy = 0.0;
    
  std::vector<float> * m_totalEnergyInLayerGroup = NULL;

  unsigned int m_nhits = 0;
  unsigned int m_leak_nhits = 0;

  int m_layerThreshold = 0;

  std::vector<float> * m_hitEnergies = NULL;
  std::vector<float> * m_raw_hitEnergies = NULL;

  std::vector<float> * m_hit_x = NULL;
  std::vector<float> * m_hit_y = NULL;
  std::vector<float> * m_hit_z = NULL;

  std::vector<int> * m_hitLayers = NULL;
  std::vector<float> * m_hitLayerThicknesses = NULL;
  std::vector<float> * m_hitLayerRadiationLengths = NULL;
  std::vector<float> * m_hitLayerIntRadiationLengths = NULL;
  std::vector<float> * m_hitLayerDistances = NULL;
  std::vector<float> * m_hitLayerSensitiveThicknesses = NULL;


  std::vector<float> * m_leak_hitEnergies = NULL;
  std::vector<float> * m_leak_hit_x = NULL;
  std::vector<float> * m_leak_hit_y = NULL;
  std::vector<float> * m_leak_hit_z = NULL;

  std::vector<int> * m_leak_hitLayers = NULL;
  std::vector<float> * m_leak_hitLayerThicknesses = NULL;
  std::vector<float> * m_leak_hitLayerRadiationLengths = NULL;
  std::vector<float> * m_leak_hitLayerIntRadiationLengths = NULL;
  std::vector<float> * m_leak_hitLayerDistances = NULL;
  std::vector<float> * m_leak_raw_hitEnergies = NULL;
  std::vector<float> * m_leak_hitLayerSensitiveThicknesses = NULL;


} ;

#endif



