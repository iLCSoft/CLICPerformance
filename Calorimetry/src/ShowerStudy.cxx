#include "ShowerStudy.h"
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/CalorimeterHit.h>

#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DetectorSelector.h"

#include "TFile.h"

using namespace lcio ;
using namespace marlin ;

DD4hep::DDRec::LayeredCalorimeterData * getExtension(unsigned int includeFlag, unsigned int excludeFlag=0) {
  
  
  DD4hep::DDRec::LayeredCalorimeterData * theExtension = 0;
  
  DD4hep::Geometry::LCDD & lcdd = DD4hep::Geometry::LCDD::getInstance();
  const std::vector< DD4hep::Geometry::DetElement>& theDetectors = DD4hep::Geometry::DetectorSelector(lcdd).detectors(  includeFlag, excludeFlag ) ;
  
  
  streamlog_out( DEBUG2 ) << " getExtension :  includeFlag: " << DD4hep::DetType( includeFlag ) << " excludeFlag: " << DD4hep::DetType( excludeFlag ) 
			  << "  found : " << theDetectors.size() << "  - first det: " << theDetectors.at(0).name() << std::endl ;
  
  if( theDetectors.size()  != 1 ){
    
    std::stringstream es ;
    es << " getExtension: selection is not unique (or empty)  includeFlag: " << DD4hep::DetType( includeFlag ) << " excludeFlag: " << DD4hep::DetType( excludeFlag ) 
       << " --- found detectors : " ;
    for( unsigned i=0, N= theDetectors.size(); i<N ; ++i ){
      es << theDetectors.at(i).name() << ", " ; 
    }
    throw std::runtime_error( es.str() ) ;
  }
  
  theExtension = theDetectors.at(0).extension<DD4hep::DDRec::LayeredCalorimeterData>();
  
  return theExtension;
}

ShowerStudy aShowerStudy;

ShowerStudy::ShowerStudy() : Processor("ShowerStudy") {
    
    // modify processor description
    _description = "ShowerStudy calculates properties of calorimeter showers" ;

    // Input collections
    registerInputCollection( LCIO::CALORIMETERHIT,
                            "CalorimeterHitCollection",
                            "CalorimeterHit collection name",
                            m_inputCalorimeterHitCollection,
                            std::string("ECALBarrel"));
    
    registerInputCollection( LCIO::MCPARTICLE,
                            "MCParticleCollectionName",
                            "Name of the MCParticle input collection",
                            m_inputMCParticleCollection,
                            std::string("MCParticle"));
    
    registerProcessorParameter( "OutputRootFileName",
                                "ROOT File name to collect plots",
                                m_rootFileName,
                                std::string("showerStudy.root"));
}


void ShowerStudy::init() {

    // Print the initial parameters
    printParameters() ;

    // Reset counters
    m_runNumber = 0 ;
    m_eventNumber = 0 ;
    
    double energyBins[11] = {0.,2.,11.,22.,51.,101.,151.,201,501,1001.,1501};
    int nbins = 10;
    m_showerProfile = new TProfile2D("showerProfile","Shower profile; Distance [mm]; E_{true} [GeV]; E/0.5 mm [GeV]",1000,1500,2000,nbins,energyBins);
    m_totalEnergyHist = new TH2F("totalEnergyHist","Total Energy;E [GeV];E_{true} [GeV];Entries/0.5 GeV",6000,0,3000,nbins,energyBins);
    
    m_showerProfileLayers = new TProfile2D("showerProfileLayers","Shower profile; Layer #; E_{true} [GeV]; E/0.5 mm [GeV]",50,0,50,nbins,energyBins);
    
    m_rootFile = new TFile(m_rootFileName.c_str(),"RECREATE");
    
    m_outputTree = new TTree("showerData","showerData");
    m_hitEnergies = new std::vector<float>();
    m_hitDistances = new std::vector<float>();
    m_hitLayers = new std::vector<int>();
    m_hitLayerThicknesses = new std::vector<float>();
    m_hitLayerRadiationLengths = new std::vector<float>();

    m_outputTree->Branch("trueEnergy",&m_trueEnergy,"trueEnergy/F");
    m_outputTree->Branch("totalEnergy",&m_totalEnergy,"totalEnergy/F");
    m_outputTree->Branch("hit_n",&m_nhits,"hit_n/i");

    m_outputTree->Branch("hit_E","std::vector< float >", m_hitEnergies);
    m_outputTree->Branch("hit_r","std::vector< float >", m_hitDistances);
    m_outputTree->Branch("hit_layer","std::vector< int >", m_hitLayers);
    m_outputTree->Branch("hit_layer_thickness","std::vector< float >", m_hitLayerThicknesses);
    m_outputTree->Branch("hit_layer_X0","std::vector< float >", m_hitLayerRadiationLengths);

    
}


void ShowerStudy::processRunHeader( LCRunHeader* run) {
	++m_runNumber ;
}

void ShowerStudy::processEvent( LCEvent* evt ) {
    
    m_nhits=0;
    m_trueEnergy=0.;
    m_totalEnergy=0.;
    m_hitEnergies->clear(); 
    m_hitDistances->clear();
    m_hitLayers->clear();
    m_hitLayerThicknesses->clear();
    m_hitLayerRadiationLengths->clear();

    //Get ECal Barrel extension by type, ignore plugs and rings 
    const DD4hep::DDRec::LayeredCalorimeterData * eCalBarrelExtension= getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::ELECTROMAGNETIC | DD4hep::DetType::BARREL), 
										     ( DD4hep::DetType::AUXILIARY  |  DD4hep::DetType::FORWARD ) );

    
    const std::vector<DD4hep::DDRec::LayeredCalorimeterStruct::Layer> & layers = eCalBarrelExtension->layers;
    
    LCCollection * mcColl =0;
    
    getCollection(mcColl,m_inputMCParticleCollection,evt);
    
    MCParticle * photonMcParticle = 0;
    //Look for a photon
    for(int m =0; m< mcColl->getNumberOfElements(); m++){
        MCParticle* mcp= dynamic_cast<MCParticle*>( mcColl->getElementAt(m) ) ;
        if(mcp->getPDG() == 22 && mcp->getGeneratorStatus()==1){
            photonMcParticle = mcp;
            break;
        }
    }

    if (photonMcParticle==0){
        streamlog_out( WARNING)<<"No generator status 1 photon found in " << m_inputMCParticleCollection.c_str() << " collection! Skipping event!" << std::endl;
        return;
    }
    
    //Check if in barrel
    double cosTheta = fabs(photonMcParticle->getMomentum()[2])/sqrt(photonMcParticle->getMomentum()[0]*photonMcParticle->getMomentum()[0]+photonMcParticle->getMomentum()[1]*photonMcParticle->getMomentum()[1]+photonMcParticle->getMomentum()[2]*photonMcParticle->getMomentum()[2]);
    
    if (cosTheta>0.7)
        return;
        
    double trueEnergy = photonMcParticle->getEnergy();
    
    LCCollection * caloColl =0;
    
    getCollection(caloColl,m_inputCalorimeterHitCollection,evt);
    UTIL::BitField64 _encoder(caloColl->getParameters().getStringVal( EVENT::LCIO::CellIDEncoding ));

    double totalEnergy=0.;
    for(int i=0; i< caloColl->getNumberOfElements(); i++){ 
	CalorimeterHit* hit = dynamic_cast<CalorimeterHit*>( caloColl->getElementAt(i) ) ;
        
        double hit_x = hit->getPosition()[0];
        double hit_y = hit->getPosition()[1];
//         double hit_z = hit->getPosition()[2];
        double radius= sqrt(hit_x*hit_x+hit_y*hit_y);
        double hitEnergy = hit->getEnergy();
        lcio::long64 cellId = long( hit->getCellID0() & 0xffffffff ) | ( long( hit->getCellID1() ) << 32 );
        _encoder.setValue(cellId);
        int layer=_encoder["layer"].value();
        
        totalEnergy = totalEnergy + hitEnergy;
        m_showerProfile->Fill(radius, trueEnergy,hitEnergy);
        m_showerProfileLayers->Fill(layer, trueEnergy,hitEnergy);
    
        m_hitDistances->push_back(radius);
        m_hitEnergies->push_back(hitEnergy);
        m_hitLayers->push_back(layer);
        
        double thickness = (layers[layer].inner_thickness+layers[layer].outer_thickness)/dd4hep::mm;
        double nRadLengths = layers[layer].inner_nRadiationLengths+layers[layer].inner_nRadiationLengths;
        
        m_hitLayerThicknesses->push_back(thickness);
        m_hitLayerRadiationLengths->push_back(nRadLengths);
        m_nhits++;
        
    }
    
    
    m_totalEnergyHist->Fill(totalEnergy,trueEnergy);

    m_totalEnergy= totalEnergy;
    m_trueEnergy = trueEnergy;
    m_outputTree->Fill();
    

}

void ShowerStudy::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
  try{
    collection = evt->getCollection( collectionName ) ;
  }
  catch(DataNotAvailableException &e){
    
    streamlog_out( DEBUG5 )<<"- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

void ShowerStudy::check(LCEvent* evt){
}

void ShowerStudy::end(){
    
    
    
    m_totalEnergyHist->Write();
    m_showerProfile->Write();
    m_showerProfileLayers->Write();
    m_outputTree->Write();
    m_rootFile->Close();
    
    
    
}
