#ifndef TrueMCintoRecoForJets_h
#define TrueMCintoRecoForJets_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <set>

using namespace lcio ;
using namespace marlin ;


class TrueMCintoRecoForJets : public Processor {
                
public:
        
    virtual Processor*  newProcessor() { return new TrueMCintoRecoForJets ; }
    
    TrueMCintoRecoForJets() ;

    TrueMCintoRecoForJets(const TrueMCintoRecoForJets&) = delete;
    TrueMCintoRecoForJets& operator=(const TrueMCintoRecoForJets&) = delete;
    
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
        
    std::string m_inputMCParticleCollection="";
    std::string m_outputRECOParticleCollection="";

    std::string m_inputRecoParticleCollection="";
    std::string m_outputRecoNoLepParticleCollection="";
 
    // Call to get collections
    void getCollection(LCCollection*&, std::string, LCEvent*);

    bool m_vetoBosonLeptons=false;
    bool m_vetoBosonLeptonsOnReco=false;
    void fillStableDaughterSet(MCParticle*, std::set<MCParticle*> &);

} ;


#endif



