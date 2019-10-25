#ifndef JetAnalyzer_h
#define JetAnalyzer_h 1

#include "marlin/Processor.h"

#include "lcio.h"

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>


#include "TH2F.h"
#include "TH1D.h"
#include "TProfile2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <map>

using namespace lcio ;
using namespace marlin ;

class TTree;


class JetAnalyzer : public Processor {
  
public:
  
  virtual Processor*  newProcessor() { return new JetAnalyzer ; }
  
  JetAnalyzer() ;

  JetAnalyzer(const JetAnalyzer&) = delete;
  JetAnalyzer& operator=(const JetAnalyzer&) = delete;
  
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
  std::string m_inputRECOParticleCollection="";
  std::string m_genjetColName="";
  std::string m_recojetColName="";


  std::string m_rootFileName="";
  
  // Run and event counters
  int m_eventNumber=0;
  int m_runNumber=0;
  //parton level information, save first quark in history
  int m_d1_mcPDGID=0;
  float m_d1_mcE=0;
  float m_d1_mcPx=0;
  float m_d1_mcPy=0;
  float m_d1_mcPz=0;
  //parton level information, save second quark in history
  int m_d2_mcPDGID=0;
  float m_d2_mcE=0;
  float m_d2_mcPx=0;
  float m_d2_mcPy=0;
  float m_d2_mcPz=0;
  //particle level information sums --> visible and invisible energy

  //only neutrinos
  float m_E_trueInv=0;
  float m_px_trueInv=0;
  float m_py_trueInv=0;
  float m_pz_trueInv=0;

  //all visible particles
  float m_E_trueVis=0;
  float m_px_trueVis=0;
  float m_py_trueVis=0;
  float m_pz_trueVis=0;

  //all reconstructed particles
  float m_E_totPFO=0;
  float m_px_totPFO=0;
  float m_py_totPFO=0;
  float m_pz_totPFO=0;

  //if extended parton quantities are saved: dibosons and their decay products
  std::vector<float>* m_trueME_E=NULL;
  std::vector<float>* m_trueME_Px=NULL;
  std::vector<float>* m_trueME_Py=NULL;
  std::vector<float>* m_trueME_Pz=NULL;
  std::vector<int>* m_trueME_PDGID=NULL;

  std::vector<float>* m_genJet_E=NULL;
  std::vector<float>* m_genJet_Px=NULL;
  std::vector<float>* m_genJet_Py=NULL;
  std::vector<float>* m_genJet_Pz=NULL;

  std::vector<float>* m_recoJet_E=NULL;
  std::vector<float>* m_recoJet_Px=NULL;
  std::vector<float>* m_recoJet_Py=NULL;
  std::vector<float>* m_recoJet_Pz=NULL;

  
  // Call to get collections
  void getCollection(LCCollection*&, std::string, LCEvent*);
  
  
  TFile * m_rootFile=NULL;
  TTree* m_outputTree=NULL;
  bool m_fillMEInfo=false;
  bool m_performDiBosonChecks=false;

} ;

#endif



