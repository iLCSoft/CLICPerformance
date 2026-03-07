#include "JetAnalyzer.h"
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>

#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DetectorSelector.h"

#include "UTIL/LCRelationNavigator.h"

#include "TFile.h"

using namespace lcio ;
using namespace marlin ;


JetAnalyzer aJetAnalyzer;

JetAnalyzer::JetAnalyzer() : Processor("JetAnalyzer") {
  
  // modify processor description
  _description = "JetAnalyzer calculates properties of calorimeter showers" ;
   

  registerInputCollection( LCIO::MCPARTICLE,
              "MCParticleCollectionName",
              "Name of the MCParticle input collection",
              m_inputMCParticleCollection,
              std::string("MCPhysicsParticles"));
  
  registerProcessorParameter( "OutputRootFileName",
              "ROOT File name to collect plots",
              m_rootFileName,
              std::string("JetAnalyzer.root"));

  registerProcessorParameter( "ProcessName",
              "Which process to run",
              m_processName,
              std::string("ProcessName"));

/*  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
              "RECOParticleCollectionName",
              "Name of the RECOParticle input collection",
              m_inputRECOParticleCollection,
              std::string("PandoraPFOs"));
*/
/*   registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
              "GenJetCollection" ,  
              "Name of the GenJet collection",
              m_genjetColName,
              std::string("GenJet_VLC"));

 registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
              "RecoJetCollection" ,  
              "Name of the RecoJet collection",
              m_recojetColName,
              std::string("RecoJet_VLC"));
*/
  registerProcessorParameter("doDiBosonChecks",
              "MC truth diboson checks for on-shell bosons",
              m_performDiBosonChecks,
              bool(false));

  registerProcessorParameter(
              "fillMEInfo" , 
              "save matrix element information",
              m_fillMEInfo,
              bool(false));

}  


void JetAnalyzer::init() {

  m_rootFile = new TFile(m_rootFileName.c_str(),"RECREATE");

  m_outputTree = new TTree("showerData","showerData");

  
  // Print the initial parameters
  printParameters() ;
  
  // Reset counters
  m_runNumber = 0 ;
  m_eventNumber = 0 ;

  if(m_fillMEInfo){
    m_trueME_E=new std::vector<float>();
    m_trueME_Px=new std::vector<float>();
    m_trueME_Py=new std::vector<float>();
    m_trueME_Pz=new std::vector<float>();
    m_trueME_PDGID=new std::vector<int>();
  }

/*  m_genJet_E   = new std::vector<float>(); 
  m_genJet_Px  = new std::vector<float>(); 
  m_genJet_Py  = new std::vector<float>(); 
  m_genJet_Pz  = new std::vector<float>(); 

  m_recoJet_E   = new std::vector<float>(); 
  m_recoJet_Px  = new std::vector<float>(); 
  m_recoJet_Py  = new std::vector<float>(); 
  m_recoJet_Pz  = new std::vector<float>(); 
*/
  m_E_trueInv  = 0;
  m_px_trueInv = 0;
  m_py_trueInv = 0;
  m_pz_trueInv = 0;

  m_E_trueVis  = 0;
  m_px_trueVis = 0;
  m_py_trueVis = 0;
  m_pz_trueVis = 0;

 /* m_E_totPFO  = 0;
  m_px_totPFO = 0;
  m_py_totPFO = 0;
  m_pz_totPFO = 0;
*/
  m_d1_mcPDGID = 0;
  m_d1_mcE   = 0;
  m_d1_mcPx  = 0;
  m_d1_mcPy  = 0;
  m_d1_mcPz  = 0;
  
  m_d2_mcPDGID = 0;
  m_d2_mcE   = 0;
  m_d2_mcPx  = 0;
  m_d2_mcPy  = 0;
  m_d2_mcPz  = 0;

  m_H1_mcPDGID = 0;
  m_H1_mcE   = 0;
  m_H1_mcPx  = 0;
  m_H1_mcPy  = 0;
  m_H1_mcPz  = 0;

  m_H2_mcPDGID = 0;
  m_H2_mcE   = 0;
  m_H2_mcPx  = 0;
  m_H2_mcPy  = 0;
  m_H2_mcPz  = 0;

  m_b1_mcPDGID = 0;
  m_b1_mcE   = 0;
  m_b1_mcPx  = 0;
  m_b1_mcPy  = 0;
  m_b1_mcPz  = 0;

  m_b2_mcPDGID = 0;
  m_b2_mcE   = 0;
  m_b2_mcPx  = 0;
  m_b2_mcPy  = 0;
  m_b2_mcPz  = 0;

  m_b3_mcPDGID = 0;
  m_b3_mcE   = 0;
  m_b3_mcPx  = 0;
  m_b3_mcPy  = 0;
  m_b3_mcPz  = 0;

  m_b4_mcPDGID = 0;
  m_b4_mcE   = 0;
  m_b4_mcPx  = 0;
  m_b4_mcPy  = 0;
  m_b4_mcPz  = 0;

  m_b5_mcPDGID = 0;
  m_b5_mcE   = 0;
  m_b5_mcPx  = 0;
  m_b5_mcPy  = 0;
  m_b5_mcPz  = 0;

  m_b6_mcPDGID = 0;
  m_b6_mcE   = 0;
  m_b6_mcPx  = 0;
  m_b6_mcPy  = 0;
  m_b6_mcPz  = 0;

  m_b7_mcPDGID = 0;
  m_b7_mcE   = 0;
  m_b7_mcPx  = 0;
  m_b7_mcPy  = 0;
  m_b7_mcPz  = 0;

  m_b8_mcPDGID = 0;
  m_b8_mcE   = 0;
  m_b8_mcPx  = 0;
  m_b8_mcPy  = 0;
  m_b8_mcPz  = 0;

//#######################################
/*  m_genJet_E ->clear(); 
  m_genJet_Px ->clear(); 
  m_genJet_Py ->clear(); 
  m_genJet_Pz ->clear(); 

  m_recoJet_E ->clear(); 
  m_recoJet_Px ->clear(); 
  m_recoJet_Py ->clear(); 
  m_recoJet_Pz ->clear(); 
*/
  if(m_fillMEInfo){
    m_trueME_E->clear();
    m_trueME_Px->clear();
    m_trueME_Py->clear();
    m_trueME_Pz->clear();
    m_trueME_PDGID->clear();
  }

  if(m_fillMEInfo){
    m_outputTree->Branch("trueME_Px", "std::vector< float >", &m_trueME_Px); 
    m_outputTree->Branch("trueME_Py", "std::vector< float >", &m_trueME_Py); 
    m_outputTree->Branch("trueME_Pz", "std::vector< float >", &m_trueME_Pz); 
    m_outputTree->Branch("trueME_E", "std::vector< float >", &m_trueME_E); 
    m_outputTree->Branch("trueME_PDGID", "std::vector< int >", &m_trueME_PDGID); 
  }

  m_outputTree->Branch("d1_mcPDGID",&m_d1_mcPDGID,"d1_mcPDGID/I");
  m_outputTree->Branch("d1_mcE",&m_d1_mcE,"d1_mcE/F");
  m_outputTree->Branch("d1_mcPx",&m_d1_mcPx,"d1_mcPx/F");
  m_outputTree->Branch("d1_mcPy",&m_d1_mcPy,"d1_mcPy/F");
  m_outputTree->Branch("d1_mcPz",&m_d1_mcPz,"d1_mcPz/F");
  
  m_outputTree->Branch("d2_mcPDGID",&m_d2_mcPDGID,"d2_mcPDGID/I");
  m_outputTree->Branch("d2_mcE",&m_d2_mcE,"d2_mcE/F");
  m_outputTree->Branch("d2_mcPx",&m_d2_mcPx,"d2_mcPx/F");
  m_outputTree->Branch("d2_mcPy",&m_d2_mcPy,"d2_mcPy/F");
  m_outputTree->Branch("d2_mcPz",&m_d2_mcPz,"d2_mcPz/F");

  //########################NEW FOR DOUBLE HIGGS ##########à
  m_outputTree->Branch("H1_mcPDGID",&m_H1_mcPDGID,"H1_mcPDGID/I");
  m_outputTree->Branch("H1_mcE",&m_H1_mcE,"H1_mcE/F");
  m_outputTree->Branch("H1_mcPx",&m_H1_mcPx,"H1_mcPx/F");
  m_outputTree->Branch("H1_mcPy",&m_H1_mcPy,"H1_mcPy/F");
  m_outputTree->Branch("H1_mcPz",&m_H1_mcPz,"H1_mcPz/F");

  m_outputTree->Branch("H2_mcPDGID",&m_H2_mcPDGID,"H2_mcPDGID/I");
  m_outputTree->Branch("H2_mcE",&m_H2_mcE,"H2_mcE/F");
  m_outputTree->Branch("H2_mcPx",&m_H2_mcPx,"H2_mcPx/F");
  m_outputTree->Branch("H2_mcPy",&m_H2_mcPy,"H2_mcPy/F");
  m_outputTree->Branch("H2_mcPz",&m_H2_mcPz,"H2_mcPz/F");

  m_outputTree->Branch("b1_mcPDGID",&m_b1_mcPDGID,"b1_mcPDGID/I");
  m_outputTree->Branch("b1_mcE",&m_b1_mcE,"b1_mcE/F");
  m_outputTree->Branch("b1_mcPx",&m_b1_mcPx,"b1_mcPx/F");
  m_outputTree->Branch("b1_mcPy",&m_b1_mcPy,"b1_mcPy/F");
  m_outputTree->Branch("b1_mcPz",&m_b1_mcPz,"b1_mcPz/F");

  m_outputTree->Branch("b2_mcPDGID",&m_b2_mcPDGID,"b2_mcPDGID/I");
  m_outputTree->Branch("b2_mcE",&m_b2_mcE,"b2_mcE/F");
  m_outputTree->Branch("b2_mcPx",&m_b2_mcPx,"b2_mcPx/F");
  m_outputTree->Branch("b2_mcPy",&m_b2_mcPy,"b2_mcPy/F");
  m_outputTree->Branch("b2_mcPz",&m_b2_mcPz,"b2_mcPz/F");

  m_outputTree->Branch("b3_mcPDGID",&m_b3_mcPDGID,"b3_mcPDGID/I");
  m_outputTree->Branch("b3_mcE",&m_b3_mcE,"b3_mcE/F");
  m_outputTree->Branch("b3_mcPx",&m_b3_mcPx,"b3_mcPx/F");
  m_outputTree->Branch("b3_mcPy",&m_b3_mcPy,"b3_mcPy/F");
  m_outputTree->Branch("b3_mcPz",&m_b3_mcPz,"b3_mcPz/F");

  m_outputTree->Branch("b4_mcPDGID",&m_b4_mcPDGID,"b4_mcPDGID/I");
  m_outputTree->Branch("b4_mcE",&m_b4_mcE,"b4_mcE/F");
  m_outputTree->Branch("b4_mcPx",&m_b4_mcPx,"b4_mcPx/F");
  m_outputTree->Branch("b4_mcPy",&m_b4_mcPy,"b4_mcPy/F");
  m_outputTree->Branch("b4_mcPz",&m_b4_mcPz,"b4_mcPz/F");

  m_outputTree->Branch("b5_mcPDGID",&m_b5_mcPDGID,"b5_mcPDGID/I");
  m_outputTree->Branch("b5_mcE",&m_b5_mcE,"b5_mcE/F");
  m_outputTree->Branch("b5_mcPx",&m_b5_mcPx,"b5_mcPx/F");
  m_outputTree->Branch("b5_mcPy",&m_b5_mcPy,"b5_mcPy/F");
  m_outputTree->Branch("b5_mcPz",&m_b5_mcPz,"b5_mcPz/F");

  m_outputTree->Branch("b6_mcPDGID",&m_b6_mcPDGID,"b6_mcPDGID/I");
  m_outputTree->Branch("b6_mcE",&m_b6_mcE,"b6_mcE/F");
  m_outputTree->Branch("b6_mcPx",&m_b6_mcPx,"b6_mcPx/F");
  m_outputTree->Branch("b6_mcPy",&m_b6_mcPy,"b6_mcPy/F");
  m_outputTree->Branch("b6_mcPz",&m_b6_mcPz,"b6_mcPz/F");
  
  m_outputTree->Branch("b7_mcPDGID",&m_b7_mcPDGID,"b7_mcPDGID/I");
  m_outputTree->Branch("b7_mcE",&m_b7_mcE,"b7_mcE/F");
  m_outputTree->Branch("b7_mcPx",&m_b7_mcPx,"b7_mcPx/F");
  m_outputTree->Branch("b7_mcPy",&m_b7_mcPy,"b7_mcPy/F");
  m_outputTree->Branch("b7_mcPz",&m_b7_mcPz,"b7_mcPz/F");

  m_outputTree->Branch("b8_mcPDGID",&m_b8_mcPDGID,"b8_mcPDGID/I");
  m_outputTree->Branch("b8_mcE",&m_b8_mcE,"b8_mcE/F");
  m_outputTree->Branch("b8_mcPx",&m_b8_mcPx,"b8_mcPx/F");
  m_outputTree->Branch("b8_mcPy",&m_b8_mcPy,"b8_mcPy/F");
  m_outputTree->Branch("b8_mcPz",&m_b8_mcPz,"b8_mcPz/F");

//###############################ENDNEW

  //true particle level, exclude neutrinos, true visible content
  m_outputTree->Branch("E_trueVis" ,&m_E_trueVis, "E_trueVis/F");
  m_outputTree->Branch("Px_trueVis",&m_px_trueVis,"Px_trueVis/F");
  m_outputTree->Branch("Py_trueVis",&m_py_trueVis,"Py_trueVis/F");
  m_outputTree->Branch("Pz_trueVis",&m_pz_trueVis,"Pz_trueVis/F");

  //true particle level, only neutrinos, true invisible content
  m_outputTree->Branch("E_trueInv" ,&m_E_trueInv, "E_trueInv/F");
  m_outputTree->Branch("Px_trueInv",&m_px_trueInv,"Px_trueInv/F");
  m_outputTree->Branch("Py_trueInv",&m_py_trueInv,"Py_trueInv/F");
  m_outputTree->Branch("Pz_trueInv",&m_pz_trueInv,"Pz_trueInv/F");
  
  //reconstructed leve
/*  m_outputTree->Branch("E_totPFO" ,&m_E_totPFO, "E_totPFO/F");
  m_outputTree->Branch("Px_totPFO",&m_px_totPFO,"Px_totPFO/F");
  m_outputTree->Branch("Py_totPFO",&m_py_totPFO,"Py_totPFO/F");
  m_outputTree->Branch("Pz_totPFO",&m_pz_totPFO,"Pz_totPFO/F");

  m_outputTree->Branch("genJetE",  "std::vector< float >", &m_genJet_E); 
  m_outputTree->Branch("genJetPx", "std::vector< float >", &m_genJet_Px); 
  m_outputTree->Branch("genJetPy", "std::vector< float >", &m_genJet_Py); 
  m_outputTree->Branch("genJetPz", "std::vector< float >", &m_genJet_Pz); 


  m_outputTree->Branch("recoJetE",  "std::vector< float >", &m_recoJet_E); 
  m_outputTree->Branch("recoJetPx", "std::vector< float >", &m_recoJet_Px); 
  m_outputTree->Branch("recoJetPy", "std::vector< float >", &m_recoJet_Py); 
  m_outputTree->Branch("recoJetPz", "std::vector< float >", &m_recoJet_Pz); 
*/

}


void JetAnalyzer::processRunHeader( LCRunHeader*) {
    ++m_runNumber ;
}

void JetAnalyzer::processEvent( LCEvent* evt ) {

  m_eventNumber=evt->getEventNumber();
  m_runNumber=evt->getRunNumber();

  m_E_trueInv  = 0;
  m_px_trueInv = 0;
  m_py_trueInv = 0;
  m_pz_trueInv = 0;

  m_E_trueVis  = 0;
  m_px_trueVis = 0;
  m_py_trueVis = 0;
  m_pz_trueVis = 0;

/*  m_E_totPFO  = 0;
  m_px_totPFO = 0;
  m_py_totPFO = 0;
  m_pz_totPFO = 0;
*/
  m_d1_mcPDGID=-10;
  m_d1_mcE=-10;
  m_d1_mcPx=-10;
  m_d1_mcPy=-10;
  m_d1_mcPz=-10;
  
  m_d2_mcPDGID=-10;
  m_d2_mcE=-10;
  m_d2_mcPx=-10;
  m_d2_mcPy=-10;
  m_d2_mcPz=-10;

//#######NEW#####################
 m_H1_mcPDGID = -10;
  m_H1_mcE   = -10;
  m_H1_mcPx  = -10;
  m_H1_mcPy  = -10;
  m_H1_mcPz  = -10;

  m_H2_mcPDGID = -10;
  m_H2_mcE   = -10;
  m_H2_mcPx  = -10;
  m_H2_mcPy  = -10;
  m_H2_mcPz  = -10;


  m_b1_mcPDGID = -10;
  m_b1_mcE   = -10;
  m_b1_mcPx  = -10;
  m_b1_mcPy  = -10;
  m_b1_mcPz  = -10;

  m_b2_mcPDGID = -10;
  m_b2_mcE   = -10;
  m_b2_mcPx  = -10;
  m_b2_mcPy  = -10;
  m_b2_mcPz  = -10;


  m_b3_mcPDGID = -10;
  m_b3_mcE   = -10;
  m_b3_mcPx  = -10;
  m_b3_mcPy  = -10;
  m_b3_mcPz  = -10;

  m_b4_mcPDGID = -10;
  m_b4_mcE   = -10;
  m_b4_mcPx  = -10;
  m_b4_mcPy  = -10;
  m_b4_mcPz  = -10;

  m_b5_mcPDGID = -10;
  m_b5_mcE   = -10;
  m_b5_mcPx  = -10;
  m_b5_mcPy  = -10;
  m_b5_mcPz  = -10;

  m_b6_mcPDGID = -10;
  m_b6_mcE   = -10;
  m_b6_mcPx  = -10;
  m_b6_mcPy  = -10;
  m_b6_mcPz  = -10;

  m_b7_mcPDGID = -10;
  m_b7_mcE   = -10;
  m_b7_mcPx  = -10;
  m_b7_mcPy  = -10;
  m_b7_mcPz  = -10;

  m_b8_mcPDGID = -10;
  m_b8_mcE   = -10;
  m_b8_mcPx  = -10;
  m_b8_mcPy  = -10;
  m_b8_mcPz  = -10;
//#######ENDNEW#####################
/*  m_genJet_E ->clear(); 
  m_genJet_Px ->clear(); 
  m_genJet_Py ->clear(); 
  m_genJet_Pz ->clear(); 

  m_recoJet_E ->clear(); 
  m_recoJet_Px ->clear(); 
  m_recoJet_Py ->clear(); 
  m_recoJet_Pz ->clear(); 
*/

  if(m_fillMEInfo){
    m_trueME_E->clear();
    m_trueME_Px->clear();
    m_trueME_Py->clear();
    m_trueME_Pz->clear();
    m_trueME_PDGID->clear();
  }


  LCCollection * mcColl =0;
  getCollection(mcColl,m_inputMCParticleCollection,evt);

  bool pass_W_boson_mass=true;
  bool pass_Z_boson_mass=true;

    if(mcColl!=NULL){
    int boson_counter=0;
    int ind_second_boson=-1;
    MCParticle* index_H1=NULL; //NEW
    MCParticle* index_H2=NULL;
    MCParticle* index_b1=NULL;
    MCParticle* index_b2=NULL;
    MCParticle* index_b3=NULL;
    MCParticle* index_b4=NULL;
    MCParticle* index_b5=NULL;
    MCParticle* index_b7=NULL;
    bool first_filled=true;
    for(int m =0; m< mcColl->getNumberOfElements(); m++){
      MCParticle* mcp= dynamic_cast<MCParticle*>( mcColl->getElementAt(m) ) ;
      //boson checks have to be done for WW,WZ and ZZ configurations
      if(m_performDiBosonChecks){
        if(abs(mcp->getPDG())==24 || mcp->getPDG()==23){
          boson_counter+=1;
          if(abs(mcp->getPDG())==24 && boson_counter<3){
            if(fabs(mcp->getMass()-80.4)>20.0){
              pass_W_boson_mass=false;
            }
          }
          if(mcp->getPDG()==23 && boson_counter<3){
            if(fabs(mcp->getMass()-91.2)>20.0){
              pass_Z_boson_mass=false;
            }
          }
        }
      }
      //first boson is the index immediately in front of the first boson
      if(boson_counter==2 && pass_W_boson_mass && ind_second_boson==-1 && m_fillMEInfo){//for ZZ passed anyway, for WW or WZ passsed if real W
        ind_second_boson=m;
        MCParticle* mcp_1st= dynamic_cast<MCParticle*>( mcColl->getElementAt(ind_second_boson-1) ) ;
        m_trueME_E->push_back(mcp_1st->getEnergy());
        m_trueME_Px->push_back(mcp_1st->getMomentum()[0]);
        m_trueME_Py->push_back(mcp_1st->getMomentum()[1]);
        m_trueME_Pz->push_back(mcp_1st->getMomentum()[2]);
        m_trueME_PDGID->push_back(mcp->getPDG());
      }
      //save all boson daughters
      if(ind_second_boson>-1 && m<(ind_second_boson+5)){
        m_trueME_E->push_back(mcp->getEnergy());
        m_trueME_Px->push_back(mcp->getMomentum()[0]);
        m_trueME_Py->push_back(mcp->getMomentum()[1]);
        m_trueME_Pz->push_back(mcp->getMomentum()[2]);
        m_trueME_PDGID->push_back(mcp->getPDG());
      }
      //fill first quark anti quark pair which appears in the history
      //this IS the hadronic W and Z (in case of a mixed di-boson event)
      if (m_processName == "dijet") {
        if(abs(mcp->getPDG())<7 && m_d1_mcE<0) {
          m_d1_mcPDGID=mcp->getPDG();
          m_d1_mcE=mcp->getEnergy();
          m_d1_mcPx=mcp->getMomentum()[0];
          m_d1_mcPy=mcp->getMomentum()[1];
          m_d1_mcPz=mcp->getMomentum()[2];
        }
        if(m_d2_mcE<0 && (abs(mcp->getPDG())<7 && mcp->getPDG()==(-m_d1_mcPDGID))){
          m_d2_mcPDGID=mcp->getPDG();
          m_d2_mcE=mcp->getEnergy();
          m_d2_mcPx=mcp->getMomentum()[0];
          m_d2_mcPy=mcp->getMomentum()[1];
          m_d2_mcPz=mcp->getMomentum()[2];
        }
      }

      if (m_processName == "HHbbbb") {
        if(mcp->getPDG()==25 && m_H1_mcE<0) {
          index_H1=mcp;
          m_H1_mcPDGID=mcp->getPDG();
          m_H1_mcE=mcp->getEnergy(); 
          m_H1_mcPx=mcp->getMomentum()[0];
          m_H1_mcPy=mcp->getMomentum()[1];
          m_H1_mcPz=mcp->getMomentum()[2];
        }
        if(m_H2_mcE<0 && (mcp->getPDG()==25 && index_H1!=mcp )){ 
          index_H2=mcp;
          m_H2_mcPDGID=mcp->getPDG();
          m_H2_mcE=mcp->getEnergy();
          m_H2_mcPx=mcp->getMomentum()[0];
          m_H2_mcPy=mcp->getMomentum()[1];
          m_H2_mcPz=mcp->getMomentum()[2];
        }
        if(abs(mcp->getPDG())==5 && m_b1_mcE<0 && mcp->getParents()[0]==index_H1) {
          m_b1_mcPDGID=mcp->getPDG();
          m_b1_mcE=mcp->getEnergy();
          m_b1_mcPx=mcp->getMomentum()[0];
          m_b1_mcPy=mcp->getMomentum()[1];
          m_b1_mcPz=mcp->getMomentum()[2];
        }
        if(m_b2_mcE<0 && (abs(mcp->getPDG())==5 && mcp->getPDG()==(-m_b1_mcPDGID)) && mcp->getParents()[0]==index_H1 ){
          m_b2_mcPDGID=mcp->getPDG();
          m_b2_mcE=mcp->getEnergy();
          m_b2_mcPx=mcp->getMomentum()[0];
          m_b2_mcPy=mcp->getMomentum()[1];
          m_b2_mcPz=mcp->getMomentum()[2];
        }
        if(abs(mcp->getPDG())==5 && m_b3_mcE<0 && mcp->getParents()[0]==index_H2) {
          m_b3_mcPDGID=mcp->getPDG();
          m_b3_mcE=mcp->getEnergy();
          m_b3_mcPx=mcp->getMomentum()[0];
          m_b3_mcPy=mcp->getMomentum()[1];
          m_b3_mcPz=mcp->getMomentum()[2];
        }
        if(m_b4_mcE<0 && (abs(mcp->getPDG())==5 && mcp->getPDG()==(-m_b3_mcPDGID))&& mcp->getParents()[0]==index_H2){
          m_b4_mcPDGID=mcp->getPDG();
          m_b4_mcE=mcp->getEnergy();
          m_b4_mcPx=mcp->getMomentum()[0];
          m_b4_mcPy=mcp->getMomentum()[1];
          m_b4_mcPz=mcp->getMomentum()[2];
        }
      }

      // if (m_processName == "bb") {
      //   if(mcp->getPDG()==5 ) {
      //     m_b1_mcPDGID=mcp->getPDG();
      //     m_b1_mcE=mcp->getEnergy();
      //     m_b1_mcPx=mcp->getMomentum()[0];
      //     m_b1_mcPy=mcp->getMomentum()[1];
      //     m_b1_mcPz=mcp->getMomentum()[2];
      //   }
      //   if(mcp->getPDG()==-5  ){
      //     m_b2_mcPDGID=mcp->getPDG();
      //     m_b2_mcE=mcp->getEnergy();
      //     m_b2_mcPx=mcp->getMomentum()[0];
      //     m_b2_mcPy=mcp->getMomentum()[1];
      //     m_b2_mcPz=mcp->getMomentum()[2];
      //   }
      // }

      // if (m_processName == "cc") {
      //   if(mcp->getPDG()==4 ) {
      //     m_b1_mcPDGID=mcp->getPDG();
      //     m_b1_mcE=mcp->getEnergy();
      //     m_b1_mcPx=mcp->getMomentum()[0];
      //     m_b1_mcPy=mcp->getMomentum()[1];
      //     m_b1_mcPz=mcp->getMomentum()[2];
      //   }
      //   if(mcp->getPDG()==-4  ){
      //     m_b2_mcPDGID=mcp->getPDG();
      //     m_b2_mcE=mcp->getEnergy();
      //     m_b2_mcPx=mcp->getMomentum()[0];
      //     m_b2_mcPy=mcp->getMomentum()[1];
      //     m_b2_mcPz=mcp->getMomentum()[2];
      //   }
      // }
      
      // if (m_processName == "qq") {
      //   if(abs(mcp->getPDG())==3 || abs(mcp->getPDG())==2 || abs(mcp->getPDG())==1) {
      //     m_b1_mcPDGID=mcp->getPDG();
      //     m_b1_mcE=mcp->getEnergy();
      //     m_b1_mcPx=mcp->getMomentum()[0];
      //     m_b1_mcPy=mcp->getMomentum()[1];
      //     m_b1_mcPz=mcp->getMomentum()[2];
      //   }
      //   if(abs(mcp->getPDG())==3 || abs(mcp->getPDG())==2 || abs(mcp->getPDG())==1){
      //     m_b2_mcPDGID=mcp->getPDG();
      //     m_b2_mcE=mcp->getEnergy();
      //     m_b2_mcPx=mcp->getMomentum()[0];
      //     m_b2_mcPy=mcp->getMomentum()[1];
      //     m_b2_mcPz=mcp->getMomentum()[2];
      //   }
      // }


      if (m_processName == "bbbb") {
        if(mcp->getPDG()==13 && m_H1_mcE<0 && mcp->getMomentum()[0]==0 && mcp->getMomentum()[1]==0 && abs(mcp->getMomentum()[2])==5000 ) {
          index_H1=mcp; 
          m_H1_mcE=mcp->getEnergy(); 
          m_H1_mcPx=mcp->getMomentum()[0];
          m_H1_mcPy=mcp->getMomentum()[1];
          m_H1_mcPz=mcp->getMomentum()[2];
        }
        if(mcp->getPDG()==-13 && m_H2_mcE<0 && mcp->getMomentum()[0]==0 && mcp->getMomentum()[1]==0 && abs(mcp->getMomentum()[2])==5000  ) { 
          index_H2=mcp; 
          m_H2_mcPDGID=mcp->getPDG();
          m_H2_mcE=mcp->getEnergy();
          m_H2_mcPx=mcp->getMomentum()[0];
          m_H2_mcPy=mcp->getMomentum()[1];
          m_H2_mcPz=mcp->getMomentum()[2];
        }
        if(mcp->getPDG()==5 && m_b1_mcE<0 && mcp->getParents()[0]->getMomentum()[0]==0 && mcp->getParents()[0]->getMomentum()[1]==0 && abs(mcp->getParents()[0]->getMomentum()[2])==5000 && (mcp->getParents()[0]==index_H1 || mcp->getParents()[0]==index_H2) ) { 
          m_b1_mcPDGID=mcp->getPDG();
          m_b1_mcE=mcp->getEnergy();
          m_b1_mcPx=mcp->getMomentum()[0];
          m_b1_mcPy=mcp->getMomentum()[1];
          m_b1_mcPz=mcp->getMomentum()[2];
          continue;
        }
        if(mcp->getPDG()==-5 && m_b2_mcE<0 && (mcp->getParents()[0]==index_H1 || mcp->getParents()[0]==index_H2) ) { 
          m_b2_mcPDGID=mcp->getPDG();
          m_b2_mcE=mcp->getEnergy();
          m_b2_mcPx=mcp->getMomentum()[0];
          m_b2_mcPy=mcp->getMomentum()[1];
          m_b2_mcPz=mcp->getMomentum()[2];
          continue;
        }
        if(mcp->getPDG()==5 && m_b3_mcE<0&& (mcp->getParents()[0]==index_H1 || mcp->getParents()[0]==index_H2)) { 
            m_b3_mcPDGID=mcp->getPDG();
            m_b3_mcE=mcp->getEnergy();
            m_b3_mcPx=mcp->getMomentum()[0];
            m_b3_mcPy=mcp->getMomentum()[1];
            m_b3_mcPz=mcp->getMomentum()[2];
            continue;
        }
        if(mcp->getPDG()==-5 && m_b4_mcE<0&& (mcp->getParents()[0]==index_H1 || mcp->getParents()[0]==index_H2) ) { 
          m_b4_mcPDGID=mcp->getPDG();
          m_b4_mcE=mcp->getEnergy();
          m_b4_mcPx=mcp->getMomentum()[0];
          m_b4_mcPy=mcp->getMomentum()[1];
          m_b4_mcPz=mcp->getMomentum()[2];
          continue;
        }
    }

    if (m_processName == "bbH") {
      if(mcp->getPDG()==13 && m_H1_mcE<0 && mcp->getMomentum()[0]==0 && mcp->getMomentum()[1]==0 && abs(mcp->getMomentum()[2])==5000 ) { 
        index_H1=mcp; 
        m_H1_mcE=mcp->getEnergy(); 
        m_H1_mcPx=mcp->getMomentum()[0];
        m_H1_mcPy=mcp->getMomentum()[1];
        m_H1_mcPz=mcp->getMomentum()[2];
      }
      if(mcp->getPDG()==25 && m_H2_mcE<0 && mcp->getParents()[0]->getMomentum()[0]==0 && mcp->getParents()[0]->getMomentum()[1]==0 && abs(mcp->getParents()[0]->getMomentum()[2])==5000  ) {    
        index_H2=mcp; 
        m_H2_mcPDGID=mcp->getPDG();
        m_H2_mcE=mcp->getEnergy();
        m_H2_mcPx=mcp->getMomentum()[0];
        m_H2_mcPy=mcp->getMomentum()[1];
        m_H2_mcPz=mcp->getMomentum()[2];
      }
      if(mcp->getPDG()==5 && m_b1_mcE<0 && mcp->getParents()[0]->getMomentum()[0]==0 && mcp->getParents()[0]->getMomentum()[1]==0 && abs(mcp->getParents()[0]->getMomentum()[2])==5000 && (mcp->getParents()[0]==index_H1 || mcp->getParents()[1]==index_H1) ) { 
        m_b1_mcPDGID=mcp->getPDG();
        m_b1_mcE=mcp->getEnergy();
        m_b1_mcPx=mcp->getMomentum()[0];
        m_b1_mcPy=mcp->getMomentum()[1];
        m_b1_mcPz=mcp->getMomentum()[2];
        continue;
      }
      if(mcp->getPDG()==-5 && m_b2_mcE<0 && mcp->getParents()[0]->getMomentum()[0]==0 && mcp->getParents()[0]->getMomentum()[1]==0 && abs(mcp->getParents()[0]->getMomentum()[2])==5000 && (mcp->getParents()[0]==index_H1 || mcp->getParents()[1]==index_H1) ) { 
        m_b2_mcPDGID=mcp->getPDG();
        m_b2_mcE=mcp->getEnergy();
        m_b2_mcPx=mcp->getMomentum()[0];
        m_b2_mcPy=mcp->getMomentum()[1];
        m_b2_mcPz=mcp->getMomentum()[2];
        continue;
      }
      if(mcp->getPDG()==5 && m_b3_mcE<0&& (mcp->getParents()[0]==index_H2)) { //&& mcp->getGeneratorStatus()==23 
        m_b3_mcPDGID=mcp->getPDG();
        m_b3_mcE=mcp->getEnergy();
        m_b3_mcPx=mcp->getMomentum()[0];
        m_b3_mcPy=mcp->getMomentum()[1];
        m_b3_mcPz=mcp->getMomentum()[2];
        continue;
      }
      if(mcp->getPDG()==-5 && m_b4_mcE<0&& (mcp->getParents()[0]==index_H2) ){ //&& mcp->getGeneratorStatus()==23
        m_b4_mcPDGID=mcp->getPDG();
        m_b4_mcE=mcp->getEnergy();
        m_b4_mcPx=mcp->getMomentum()[0];
        m_b4_mcPy=mcp->getMomentum()[1];
        m_b4_mcPz=mcp->getMomentum()[2];
        continue;
      }
    }
    if (m_processName == "bc") {
      if((abs(mcp->getPDG())==4 || abs(mcp->getPDG())==5) && mcp->getGeneratorStatus()==23 && first_filled==true) {
        m_b1_mcPDGID=mcp->getPDG();
        m_b1_mcE=mcp->getEnergy();
        m_b1_mcPx=mcp->getMomentum()[0];
        m_b1_mcPy=mcp->getMomentum()[1];
        m_b1_mcPz=mcp->getMomentum()[2];
        first_filled=false;
        continue;
      }
      if((abs(mcp->getPDG())==4 || abs(mcp->getPDG())==5) && mcp->getGeneratorStatus()==23 && first_filled==false ){
        m_b2_mcPDGID=mcp->getPDG();
        m_b2_mcE=mcp->getEnergy();
        m_b2_mcPx=mcp->getMomentum()[0];
        m_b2_mcPy=mcp->getMomentum()[1];
        m_b2_mcPz=mcp->getMomentum()[2];
      }
    }

    if (m_processName == "qqH") {
      if(mcp->getPDG()==13 && m_H1_mcE<0 && mcp->getMomentum()[0]==0 && mcp->getMomentum()[1]==0 && abs(mcp->getMomentum()[2])==1500 ){ 
        index_H1=mcp; 
        m_H1_mcE=mcp->getEnergy(); 
        m_H1_mcPx=mcp->getMomentum()[0];
        m_H1_mcPy=mcp->getMomentum()[1];
        m_H1_mcPz=mcp->getMomentum()[2];
      }
      if(mcp->getPDG()==25 && m_H2_mcE<0 && mcp->getParents()[0]->getMomentum()[0]==0 && mcp->getParents()[0]->getMomentum()[1]==0 && abs(mcp->getParents()[0]->getMomentum()[2])==1500  ) {    
        index_H2=mcp; 
        m_H2_mcPDGID=mcp->getPDG();
        m_H2_mcE=mcp->getEnergy();
        m_H2_mcPx=mcp->getMomentum()[0];
        m_H2_mcPy=mcp->getMomentum()[1];
        m_H2_mcPz=mcp->getMomentum()[2];
      }
      if((abs(mcp->getPDG())==5||abs(mcp->getPDG())==4) && m_b1_mcE<0 && mcp->getParents()[0]->getMomentum()[0]==0 && mcp->getParents()[0]->getMomentum()[1]==0 && abs(mcp->getParents()[0]->getMomentum()[2])==1500 && (mcp->getParents()[0]==index_H1 || mcp->getParents()[1]==index_H1) ) { 
        index_b1=mcp;
        m_b1_mcPDGID=mcp->getPDG();
        m_b1_mcE=mcp->getEnergy();
        m_b1_mcPx=mcp->getMomentum()[0];
        m_b1_mcPy=mcp->getMomentum()[1];
        m_b1_mcPz=mcp->getMomentum()[2];
        continue;
      }
      if((abs(mcp->getPDG())==5||abs(mcp->getPDG())==4) && m_b2_mcE<0 && mcp->getParents()[0]->getMomentum()[0]==0 && mcp->getParents()[0]->getMomentum()[1]==0 && abs(mcp->getParents()[0]->getMomentum()[2])==1500 && (mcp->getParents()[0]==index_H1 || mcp->getParents()[1]==index_H1) && mcp!=index_b1 ) { 
        m_b2_mcPDGID=mcp->getPDG();
        m_b2_mcE=mcp->getEnergy();
        m_b2_mcPx=mcp->getMomentum()[0];
        m_b2_mcPy=mcp->getMomentum()[1];
        m_b2_mcPz=mcp->getMomentum()[2];
        continue;
      }
      if((abs(mcp->getPDG())==5||abs(mcp->getPDG())==4) && m_b5_mcE<0 && mcp->getParents()[0]->getMomentum()[0]==0 && mcp->getParents()[0]->getMomentum()[1]==0 && abs(mcp->getParents()[0]->getMomentum()[2])==1500 && (mcp->getParents()[0]==index_H1 || mcp->getParents()[1]==index_H1) && mcp!=index_b1 && mcp!=index_b2 ) { 
        index_b5=mcp;
        m_b5_mcPDGID=mcp->getPDG();
        m_b5_mcE=mcp->getEnergy();
        m_b5_mcPx=mcp->getMomentum()[0];
        m_b5_mcPy=mcp->getMomentum()[1];
        m_b5_mcPz=mcp->getMomentum()[2];
        continue;
      }
      if((abs(mcp->getPDG())==5||abs(mcp->getPDG())==4) && m_b6_mcE<0 && mcp->getParents()[0]->getMomentum()[0]==0 && mcp->getParents()[0]->getMomentum()[1]==0 && abs(mcp->getParents()[0]->getMomentum()[2])==1500 && (mcp->getParents()[0]==index_H1 || mcp->getParents()[1]==index_H1) && mcp!=index_b1 && mcp!=index_b2 && mcp!=index_b5) { 
        m_b6_mcPDGID=mcp->getPDG();
        m_b6_mcE=mcp->getEnergy();
        m_b6_mcPx=mcp->getMomentum()[0];
        m_b6_mcPy=mcp->getMomentum()[1];
        m_b6_mcPz=mcp->getMomentum()[2];
        continue;
      }
      if(mcp->getPDG()==5 && m_b3_mcE<0&& (mcp->getParents()[0]==index_H2)) {
        index_b3=mcp;
        m_b3_mcPDGID=mcp->getPDG();
        m_b3_mcE=mcp->getEnergy();
        m_b3_mcPx=mcp->getMomentum()[0];
        m_b3_mcPy=mcp->getMomentum()[1];
        m_b3_mcPz=mcp->getMomentum()[2];
        continue;
      }
      if(mcp->getPDG()==-5 && m_b4_mcE<0&& (mcp->getParents()[0]==index_H2) ){
        index_b4=mcp;
        m_b4_mcPDGID=mcp->getPDG();
        m_b4_mcE=mcp->getEnergy();
        m_b4_mcPx=mcp->getMomentum()[0];
        m_b4_mcPy=mcp->getMomentum()[1];
        m_b4_mcPz=mcp->getMomentum()[2];
        continue;
      }
      if((abs(mcp->getPDG())==5||abs(mcp->getPDG())==4) && m_b7_mcE<0&& (mcp->getParents()[0]==index_H2) && mcp!= index_b3 && mcp!= index_b4 ) {
        index_b7=mcp;
        m_b7_mcPDGID=mcp->getPDG();
        m_b7_mcE=mcp->getEnergy();
        m_b7_mcPx=mcp->getMomentum()[0];
        m_b7_mcPy=mcp->getMomentum()[1];
        m_b7_mcPz=mcp->getMomentum()[2];
        continue;
      }
      if((abs(mcp->getPDG())==5||abs(mcp->getPDG())==4) && m_b8_mcE<0&& (mcp->getParents()[0]==index_H2)&& mcp!= index_b3 && mcp!= index_b4 && mcp!=index_b7 ){
        m_b8_mcPDGID=mcp->getPDG();
        m_b8_mcE=mcp->getEnergy();
        m_b8_mcPx=mcp->getMomentum()[0];
        m_b8_mcPy=mcp->getMomentum()[1];
        m_b8_mcPz=mcp->getMomentum()[2];
        continue;
      }
    }

    if (m_processName == "qqqqlnu") {
      if((abs(mcp->getPDG())==5||abs(mcp->getPDG())==4||abs(mcp->getPDG())==3||abs(mcp->getPDG())==2||abs(mcp->getPDG())==1) && m_b1_mcPDGID == -10) {
        index_b1=mcp;
        m_b1_mcPDGID=mcp->getPDG();
        m_b1_mcE=mcp->getEnergy();
        m_b1_mcPx=mcp->getMomentum()[0];
        m_b1_mcPy=mcp->getMomentum()[1];
        m_b1_mcPz=mcp->getMomentum()[2];
      }
      if((abs(mcp->getPDG())==5||abs(mcp->getPDG())==4||abs(mcp->getPDG())==3||abs(mcp->getPDG())==2||abs(mcp->getPDG())==1) &&
        mcp!=index_b1 && m_b2_mcPDGID == -10) {
          index_b2=mcp;
        m_b2_mcPDGID=mcp->getPDG();
        m_b2_mcE=mcp->getEnergy();
        m_b2_mcPx=mcp->getMomentum()[0];
        m_b2_mcPy=mcp->getMomentum()[1];
        m_b2_mcPz=mcp->getMomentum()[2];
      }
      if((abs(mcp->getPDG())==5||abs(mcp->getPDG())==4||abs(mcp->getPDG())==3||abs(mcp->getPDG())==2||abs(mcp->getPDG())==1) && 
        mcp!=index_b1 && mcp!=index_b2 && m_b3_mcPDGID == -10) {
        index_b3=mcp; 
        m_b3_mcPDGID=mcp->getPDG();
        m_b3_mcE=mcp->getEnergy();
        m_b3_mcPx=mcp->getMomentum()[0];
        m_b3_mcPy=mcp->getMomentum()[1];
        m_b3_mcPz=mcp->getMomentum()[2];
      }
      if (m_b3_mcE>0 && (abs(mcp->getPDG())==5||abs(mcp->getPDG())==4||abs(mcp->getPDG())==3||abs(mcp->getPDG())==2||abs(mcp->getPDG())==1) && 
        mcp!=index_b1 && mcp!=index_b2 && mcp!=index_b3 && m_b4_mcPDGID == -10) {
        index_b4=mcp;     
	      m_b4_mcPDGID=mcp->getPDG();
        m_b4_mcE=mcp->getEnergy();
        m_b4_mcPx=mcp->getMomentum()[0];
        m_b4_mcPy=mcp->getMomentum()[1];
        m_b4_mcPz=mcp->getMomentum()[2];
      }
    }

      if(mcp->getGeneratorStatus()==1){//visible sum of stable particles --> take neutrinos out
        if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 && abs(mcp->getPDG())!=16){
          m_E_trueVis+=mcp->getEnergy();
          m_px_trueVis+=mcp->getMomentum()[0];
          m_py_trueVis+=mcp->getMomentum()[1];
          m_pz_trueVis+=mcp->getMomentum()[2];
        }else{
          m_E_trueInv+=mcp->getEnergy();
          m_px_trueInv+=mcp->getMomentum()[0];
          m_py_trueInv+=mcp->getMomentum()[1];
          m_pz_trueInv+=mcp->getMomentum()[2];
        }
      }
    }
  }

/*  LCCollection* recoparticlecol = NULL;
  recoparticlecol = evt->getCollection(m_inputRECOParticleCollection) ;
  if(recoparticlecol!=NULL){
    //PandoraCandidate loop
    for(int i=0;i<recoparticlecol->getNumberOfElements();i++){
      ReconstructedParticle* pandorapart = dynamic_cast<ReconstructedParticle*>(recoparticlecol->getElementAt(i));
      m_E_totPFO+=pandorapart->getEnergy();
      m_px_totPFO+=pandorapart->getMomentum()[0];
      m_py_totPFO+=pandorapart->getMomentum()[1];
      m_pz_totPFO+=pandorapart->getMomentum()[2];
    }
  }
*/
/*  LCCollection* recojets = NULL;
  recojets = evt->getCollection(m_recojetColName) ;
  if(recojets!=NULL){
    //PandoraCandidate loop
    for(int i=0;i<recojets->getNumberOfElements();i++){
      ReconstructedParticle* recojet = dynamic_cast<ReconstructedParticle*>(recojets->getElementAt(i));
      m_recoJet_E->push_back(recojet->getEnergy());
      m_recoJet_Px->push_back(recojet->getMomentum()[0]);
      m_recoJet_Py->push_back(recojet->getMomentum()[1]);
      m_recoJet_Pz->push_back(recojet->getMomentum()[2]);
    }
  }
*/
/*
  LCCollection* genjets = NULL;
  genjets = evt->getCollection(m_genjetColName) ;
  if(genjets!=NULL){
    //PandoraCandidate loop
    for(int i=0;i<genjets->getNumberOfElements();i++){
      ReconstructedParticle* genjet = dynamic_cast<ReconstructedParticle*>(genjets->getElementAt(i));
      m_genJet_E->push_back(genjet->getEnergy());
      m_genJet_Px->push_back(genjet->getMomentum()[0]);
      m_genJet_Py->push_back(genjet->getMomentum()[1]);
      m_genJet_Pz->push_back(genjet->getMomentum()[2]);
    }
  }
*/

  if(pass_W_boson_mass && pass_Z_boson_mass){
    m_outputTree->Fill();
  }
}

void JetAnalyzer::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
  try{
    collection = evt->getCollection( collectionName ) ;
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG5 )<<"- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

void JetAnalyzer::check(LCEvent*){
}

void JetAnalyzer::end(){
  
  m_rootFile->cd();  
 
  m_rootFile->Write();
  m_rootFile->Close();

}

