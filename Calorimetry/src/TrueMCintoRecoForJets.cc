#include "TrueMCintoRecoForJets.h"
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/ReconstructedParticleImpl.h>

#include "DDRec/DetectorData.h"

#include "TLorentzVector.h"

using namespace lcio ;
using namespace marlin ;


TrueMCintoRecoForJets aTrueMCintoRecoForJets;


TrueMCintoRecoForJets::TrueMCintoRecoForJets() : Processor("TrueMCintoRecoForJets") {
  
    // modify processor description
    _description = "TrueMCintoRecoForJets calculates properties of calorimeter showers" ;
   

    registerInputCollection( LCIO::MCPARTICLE,
                            "MCParticleInputCollectionName",
                            "Name of the MCParticle input collection",
                            m_inputMCParticleCollection,
                            std::string("MCPhysicsParticles"));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "RECOParticleCollectionName",
                            "Name of the RECOParticle output collection",
                            m_outputRECOParticleCollection,
                            std::string("MCParticlePandoraPFOs"));

    registerProcessorParameter("vetoBosonLeptons",
			       "leptons from boson veto flag for MC truth",
			       m_vetoBosonLeptons,
			       bool(false));

    registerProcessorParameter("vetoBosonLeptonsOnReco",
			       "leptons from boson veto flag for Reco",
			       m_vetoBosonLeptonsOnReco,
			       bool(false));

    registerProcessorParameter("ignoreNeutrinosInMCJets",
			       "remove neutrinos prior to MC Jet Clustering",
			       m_ignoreNeutrinosInMCJets,
			       bool(true));

    registerProcessorParameter("cosAngle_pfo_lepton",
			       "matching angle between PFOs and MC leptons from vector bosons",
			       m_cosAngle_pfo_lepton,
			       float(0.995));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "RecoParticleInputCollectionName",
                            "Name of the RecoParticle input collection",
                            m_inputRecoParticleCollection,
                            std::string("TightSelectedPandoraPFOs"));

    registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                            "RecoParticleNoLeptonCollectionName",
                            "Name of the RecoParticleNoLepton output collection",
			     m_outputRecoNoLepParticleCollection,
                            std::string("PandoraPFOsNoLeptons"));
  
}  


void TrueMCintoRecoForJets::init() {
}


void TrueMCintoRecoForJets::processRunHeader( LCRunHeader*) {

}

void TrueMCintoRecoForJets::processEvent( LCEvent* evt ) {
  
    LCCollection * mcColl =0;
    getCollection(mcColl,m_inputMCParticleCollection,evt);
    LCCollectionVec *reccolForJets = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    LCCollectionVec *reccolNoLep = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    
    if( mcColl != 0 ){
      int nMCP = mcColl->getNumberOfElements();
      
      std::set<MCParticle*> boson_daughtersFunc;
 
      int ind_MCLep=-1;//semi leptonically decaying W-lepton (e,mu), or first lepton (e,mu) for Z
      int ind_MCLep2=-1;// in case of Z second lepton (e,mu)

      for(int m=0;m<nMCP;m++){
	MCParticle *mcp = static_cast<MCParticle*>(mcColl->getElementAt(m));
	//in order to catch events where the lepton branches off an FSR photon
	if(m_vetoBosonLeptons && (abs(mcp->getPDG())==24 || mcp->getPDG()==23) && mcp->getDaughters().size()>1 && (abs(mcp->getDaughters()[0]->getPDG())>6 && abs(mcp->getDaughters()[0]->getPDG())<17)){
	  fillStableDaughterSet(mcp, boson_daughtersFunc);
	}
	if(mcp->getGeneratorStatus()!=1){
	  continue;
	}
	if (boson_daughtersFunc.count(mcp) != 0)
	{
	  if(ind_MCLep==-1 && (abs(mcp->getPDG())==11 || abs(mcp->getPDG())==13)){
	    ind_MCLep=m;	  
	  }
	  //only check for second lepton in case it is not already the first lepton
	  if((ind_MCLep!=m && ind_MCLep2==-1) && (abs(mcp->getPDG())==11 || abs(mcp->getPDG())==13)){
	    ind_MCLep2=m;
	  }
	  continue;
	}
	if(m_ignoreNeutrinosInMCJets && (fabs(mcp->getPDG())==12 || fabs(mcp->getPDG())==14 || fabs(mcp->getPDG())==16)){
	  continue;
	}
	ReconstructedParticleImpl* truePartIntoReco = new ReconstructedParticleImpl;
	truePartIntoReco->setMomentum(mcp->getMomentum());
	truePartIntoReco->setType(mcp->getPDG());
	truePartIntoReco->setEnergy(mcp->getEnergy());
	truePartIntoReco->setMass(mcp->getMass());
	truePartIntoReco->setCharge(mcp->getCharge());
	reccolForJets->addElement(truePartIntoReco); 
      }

      
      LCCollection * recoColl =0;
      getCollection(recoColl,m_inputRecoParticleCollection,evt);
      int nReco = recoColl->getNumberOfElements();
      for(int r=0;r<nReco;r++){
	ReconstructedParticle *reco = static_cast<ReconstructedParticle*>(recoColl->getElementAt(r));
	if( ind_MCLep!=-1){
	  MCParticle *mcpLep = static_cast<MCParticle*>(mcColl->getElementAt(ind_MCLep));
	  TLorentzVector trueLep(0,0,0,0);
	  trueLep.SetPxPyPzE(mcpLep->getMomentum()[0],mcpLep->getMomentum()[1],mcpLep->getMomentum()[2],mcpLep->getEnergy());	   
	  //we don't want to check lepton reconstruction/lepton finders
	  //veto particles close to the true lepton from the bosons, ignore tau's for now
	  //take true MC lepton direction, veto all reconstructed particles around that direction within acos(alpha)>0.90
	  TLorentzVector recoVec(0,0,0,0);
	  recoVec.SetPxPyPzE(reco->getMomentum()[0],reco->getMomentum()[1],reco->getMomentum()[2],reco->getEnergy());
	  if(cos(recoVec.Angle(trueLep.Vect()))>m_cosAngle_pfo_lepton && m_vetoBosonLeptonsOnReco){
	    continue;
	  }
	}
	if( ind_MCLep2!=-1){
	  MCParticle *mcpLep2 = static_cast<MCParticle*>(mcColl->getElementAt(ind_MCLep2));
	  TLorentzVector trueLep2(0,0,0,0);
	  trueLep2.SetPxPyPzE(mcpLep2->getMomentum()[0],mcpLep2->getMomentum()[1],mcpLep2->getMomentum()[2],mcpLep2->getEnergy());	   
	  //we don't want to check lepton reconstruction/lepton finders
	  //veto particles close to the true lepton from the bosons, ignore tau's for now
	  //take true MC lepton direction, veto all reconstructed particles around that direction within acos(alpha)>0.90
	  TLorentzVector recoVec(0,0,0,0);
	  recoVec.SetPxPyPzE(reco->getMomentum()[0],reco->getMomentum()[1],reco->getMomentum()[2],reco->getEnergy());
	  if(cos(recoVec.Angle(trueLep2.Vect()))>m_cosAngle_pfo_lepton && m_vetoBosonLeptonsOnReco){
	    continue;
	  }
	}
	ReconstructedParticleImpl* recoNoLep = new ReconstructedParticleImpl;
	recoNoLep->setMomentum(reco->getMomentum());
	recoNoLep->setEnergy(reco->getEnergy());
	recoNoLep->setType(reco->getType());
	recoNoLep->setCovMatrix(reco->getCovMatrix());
	recoNoLep->setMass(reco->getMass());
	recoNoLep->setCharge(reco->getCharge());
	recoNoLep->setParticleIDUsed(reco->getParticleIDUsed());
	recoNoLep->setGoodnessOfPID(reco->getGoodnessOfPID());
	recoNoLep->setStartVertex(reco->getStartVertex());
	reccolNoLep->addElement(recoNoLep); 
      }
    }
    
    evt->addCollection(reccolNoLep,m_outputRecoNoLepParticleCollection);
    evt->addCollection(reccolForJets,m_outputRECOParticleCollection);
    
}

void TrueMCintoRecoForJets::fillStableDaughterSet(MCParticle* mcp, std::set<MCParticle*> &stableDaughterSet){
  if(mcp->getGeneratorStatus()==1){
    stableDaughterSet.insert(mcp);
  }else if (mcp->getGeneratorStatus()==0){
    return;
  }
  for(unsigned int d=0;d<mcp->getDaughters().size();d++){
    fillStableDaughterSet(mcp->getDaughters()[d], stableDaughterSet);
  }
}

void TrueMCintoRecoForJets::check(LCEvent*){
}

void TrueMCintoRecoForJets::end(){

}


void TrueMCintoRecoForJets::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
  try{
    collection = evt->getCollection( collectionName ) ;
  }
  catch(DataNotAvailableException &e){
    
    streamlog_out( DEBUG5 )<<"- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}
