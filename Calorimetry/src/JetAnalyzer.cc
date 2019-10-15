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

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
              "RECOParticleCollectionName",
              "Name of the RECOParticle input collection",
              m_inputRECOParticleCollection,
              std::string("PandoraPFOs"));

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
              "GenJetCollection" ,  
              "Name of the GenJet collection",
              m_genjetColName,
              std::string("GenJet_VLC"));

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
              "RecoJetCollection" ,  
              "Name of the RecoJet collection",
              m_recojetColName,
              std::string("RecoJet_VLC"));

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

  m_genJet_E   = new std::vector<float>(); 
  m_genJet_Px  = new std::vector<float>(); 
  m_genJet_Py  = new std::vector<float>(); 
  m_genJet_Pz  = new std::vector<float>(); 

  m_recoJet_E   = new std::vector<float>(); 
  m_recoJet_Px  = new std::vector<float>(); 
  m_recoJet_Py  = new std::vector<float>(); 
  m_recoJet_Pz  = new std::vector<float>(); 

  m_E_trueInv  = 0;
  m_px_trueInv = 0;
  m_py_trueInv = 0;
  m_pz_trueInv = 0;

  m_E_trueAll  = 0;
  m_px_trueAll = 0;
  m_py_trueAll = 0;
  m_pz_trueAll = 0;

  m_E_totPFO  = 0;
  m_px_totPFO = 0;
  m_py_totPFO = 0;
  m_pz_totPFO = 0;

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

  m_genJet_E ->clear(); 
  m_genJet_Px ->clear(); 
  m_genJet_Py ->clear(); 
  m_genJet_Pz ->clear(); 

  m_recoJet_E ->clear(); 
  m_recoJet_Px ->clear(); 
  m_recoJet_Py ->clear(); 
  m_recoJet_Pz ->clear(); 

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

  //true particle level, exclude neutrinos
  m_outputTree->Branch("E_trueAll" ,&m_E_trueAll, "E_trueAll/F");
  m_outputTree->Branch("Px_trueAll",&m_px_trueAll,"Px_trueAll/F");
  m_outputTree->Branch("Py_trueAll",&m_py_trueAll,"Py_trueAll/F");
  m_outputTree->Branch("Pz_trueAll",&m_pz_trueAll,"Pz_trueAll/F");

  //true particle level, only neutrinos
  m_outputTree->Branch("E_trueInv" ,&m_E_trueInv, "E_trueInv/F");
  m_outputTree->Branch("Px_trueInv",&m_px_trueInv,"Px_trueInv/F");
  m_outputTree->Branch("Py_trueInv",&m_py_trueInv,"Py_trueInv/F");
  m_outputTree->Branch("Pz_trueInv",&m_pz_trueInv,"Pz_trueInv/F");
  
  //reconstructed leve
  m_outputTree->Branch("E_totPFO" ,&m_E_totPFO, "E_totPFO/F");
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

  m_E_trueAll  = 0;
  m_px_trueAll = 0;
  m_py_trueAll = 0;
  m_pz_trueAll = 0;

  m_E_totPFO  = 0;
  m_px_totPFO = 0;
  m_py_totPFO = 0;
  m_pz_totPFO = 0;

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

  m_genJet_E ->clear(); 
  m_genJet_Px ->clear(); 
  m_genJet_Py ->clear(); 
  m_genJet_Pz ->clear(); 

  m_recoJet_E ->clear(); 
  m_recoJet_Px ->clear(); 
  m_recoJet_Py ->clear(); 
  m_recoJet_Pz ->clear(); 


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
      if(mcp->getGeneratorStatus()==1){//visible sum of stable particles --> take neutrinos out
        if(abs(mcp->getPDG())!=12 && abs(mcp->getPDG())!=14 && abs(mcp->getPDG())!=16){
          m_E_trueAll+=mcp->getEnergy();
          m_px_trueAll+=mcp->getMomentum()[0];
          m_py_trueAll+=mcp->getMomentum()[1];
          m_pz_trueAll+=mcp->getMomentum()[2];
        }else{
          m_E_trueInv+=mcp->getEnergy();
          m_px_trueInv+=mcp->getMomentum()[0];
          m_py_trueInv+=mcp->getMomentum()[1];
          m_pz_trueInv+=mcp->getMomentum()[2];
        }
      }
    }
  }

  LCCollection* recoparticlecol = NULL;
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

  LCCollection* recojets = NULL;
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

