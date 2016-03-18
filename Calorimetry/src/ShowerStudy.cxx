#include "ShowerStudy.h"
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/CalorimeterHit.h>
#include <EVENT/SimCalorimeterHit.h>

#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DetectorSelector.h"

#include "UTIL/LCRelationNavigator.h"

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
    
    registerInputCollection( LCIO::CALORIMETERHIT,
                            "LeakageCalorimeterHitCollection",
                            "CalorimeterHit collection name for leakage energy",
                            m_inputLeakageCalorimeterHitCollection,
                            std::string("HCALBarrel"));
    
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
    
//     double energyBins[11] = {0.,2.,11.,22.,51.,101.,151.,201,501,1001.,1501};
    double energyBins[22] = {0.,2.,9.,11.,19.,21.,49.,51.,99.,101.,149,151.,199.,201.,499.,501.,999,1001.,1499.,1501,1999.,2001.};

//     int nbins = 10;
    int nbins = 21;

    m_totalEnergyHist = new TH2F("totalEnergyHist","Total Energy;E [GeV];E_{true} [GeV];Entries/0.5 GeV",6000,0,3000,nbins,energyBins);
    
    m_showerHist = new TH2F("showerHist","Shower profile; Distance [mm]; E_{true} [GeV]; E/1 mm [GeV]",2500,1500,4000,nbins,energyBins);
    m_showerHistLayers = new TH2F("showerHistLayers","Shower profile; Layer No; E_{true} [GeV]; E/0.5 mm [GeV]",80,0,80,nbins,energyBins);
    m_showerHistX0 = new TH2F("showerHistX0","Shower profile; Number of X0; E_{true} [GeV]; E/0.1 [GeV]",100,0,50,nbins,energyBins);
    
    m_raw_showerHist = new TH2F("rawShowerHist","Shower profile; Distance [mm]; E_{true} [GeV]; E/1 mm [GeV]",2500,1500,4000,nbins,energyBins);
    m_raw_showerHistLayers = new TH2F("rawShowerHistLayers","Shower profile; Layer No; E_{true} [GeV]; E/0.5 mm [GeV]",80,0,80,nbins,energyBins);
    m_raw_showerHistX0 = new TH2F("rawShowerHistX0","Shower profile; Number of X0; E_{true} [GeV]; E/0.1 [GeV]",100,0,50,nbins,energyBins);

    m_leakageProfile = new TProfile("leakageProfile","Leakage vs E_{true};  E_{true} [GeV]; E_{HCal}/(E_{ECal}+E_{HCal})",nbins,energyBins);

    
    m_rootFile = new TFile(m_rootFileName.c_str(),"RECREATE");
    
    m_outputTree = new TTree("showerData","showerData");
    m_hitEnergies = new std::vector<float>();
    m_raw_hitEnergies = new std::vector<float>();


    m_hit_x = new std::vector<float>();
    m_hit_y = new std::vector<float>();
    m_hit_z = new std::vector<float>();

    m_hitLayers = new std::vector<int>();
    m_hitLayerThicknesses = new std::vector<float>();
    m_hitLayerRadiationLengths = new std::vector<float>();
    m_hitLayerIntRadiationLengths = new std::vector<float>();
    m_hitLayerDistances = new std::vector<float>();
    
    m_leak_hit_x = new std::vector<float>();
    m_leak_hit_y = new std::vector<float>();
    m_leak_hit_z = new std::vector<float>();
    
    m_leak_hitEnergies = new std::vector<float>();
    m_leak_raw_hitEnergies = new std::vector<float>();
    
    m_leak_hitLayers = new std::vector<int>();
    m_leak_hitLayerThicknesses = new std::vector<float>();
    m_leak_hitLayerRadiationLengths = new std::vector<float>();
    m_leak_hitLayerIntRadiationLengths = new std::vector<float>();
    m_leak_hitLayerDistances = new std::vector<float>();

    m_outputTree->Branch("trueEnergy",&m_trueEnergy,"trueEnergy/F");
    m_outputTree->Branch("totalEnergy",&m_totalEnergy,"totalEnergy/F");
    m_outputTree->Branch("hit_n",&m_nhits,"hit_n/i");
    m_outputTree->Branch("totalLeakEnergy",&m_totalLeakEnergy,"totalLeakEnergy/F");
    m_outputTree->Branch("hit_leak_n",&m_leak_nhits,"hit_leak_n/i");

    m_outputTree->Branch("hit_E","std::vector< float >", m_hitEnergies);
    m_outputTree->Branch("hit_rawE","std::vector< float >", m_raw_hitEnergies);

    m_outputTree->Branch("hit_x","std::vector< float >", m_hit_x);
    m_outputTree->Branch("hit_y","std::vector< float >", m_hit_y);
    m_outputTree->Branch("hit_z","std::vector< float >", m_hit_z);

    m_outputTree->Branch("hit_layer","std::vector< int >", m_hitLayers);
    m_outputTree->Branch("hit_layer_thickness","std::vector< float >", m_hitLayerThicknesses);
    m_outputTree->Branch("hit_layer_X0","std::vector< float >", m_hitLayerRadiationLengths);
    m_outputTree->Branch("hit_layer_intX0","std::vector< float >", m_hitLayerIntRadiationLengths);
    m_outputTree->Branch("hit_layer_distances","std::vector< float >", m_hitLayerDistances);

    m_outputTree->Branch("hit_leak_E","std::vector< float >", m_leak_hitEnergies);
    m_outputTree->Branch("hit_leak_rawE","std::vector< float >", m_leak_raw_hitEnergies);

    m_outputTree->Branch("hit_leak_x","std::vector< float >", m_leak_hit_x);
    m_outputTree->Branch("hit_leak_y","std::vector< float >", m_leak_hit_y);
    m_outputTree->Branch("hit_leak_z","std::vector< float >", m_leak_hit_z);

    m_outputTree->Branch("hit_leak_layer","std::vector< int >", m_leak_hitLayers);
    m_outputTree->Branch("hit_leak_layer_thickness","std::vector< float >", m_leak_hitLayerThicknesses);
    m_outputTree->Branch("hit_leak_layer_X0","std::vector< float >", m_leak_hitLayerRadiationLengths);
    m_outputTree->Branch("hit_leak_layer_intX0","std::vector< float >", m_leak_hitLayerIntRadiationLengths);
    m_outputTree->Branch("hit_leak_layer_distances","std::vector< float >", m_leak_hitLayerDistances);
    
}


void ShowerStudy::processRunHeader( LCRunHeader* run) {
	++m_runNumber ;
}

void ShowerStudy::processEvent( LCEvent* evt ) {
    
    m_nhits=0;
    m_trueEnergy=0.;
    m_totalEnergy=0.;
    m_totalLeakEnergy=0.;
    m_leak_nhits=0;

    m_hitEnergies->clear(); 
    m_raw_hitEnergies->clear(); 
    m_hit_x->clear();
    m_hit_y->clear();
    m_hit_z->clear();
    m_hitLayers->clear();
    m_hitLayerThicknesses->clear();
    m_hitLayerRadiationLengths->clear();
    m_hitLayerIntRadiationLengths->clear();
    m_hitLayerDistances->clear();
    
    m_leak_hitEnergies->clear(); 
    m_leak_raw_hitEnergies->clear(); 

    m_leak_hit_x->clear();
    m_leak_hit_y->clear();
    m_leak_hit_z->clear();
    m_leak_hitLayers->clear();
    m_leak_hitLayerThicknesses->clear();
    m_leak_hitLayerRadiationLengths->clear();
    m_leak_hitLayerIntRadiationLengths->clear();
    m_leak_hitLayerDistances->clear();

    //Get ECal Barrel extension by type, ignore plugs and rings 
    const DD4hep::DDRec::LayeredCalorimeterData * eCalBarrelExtension= getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::ELECTROMAGNETIC | DD4hep::DetType::BARREL), 
										     ( DD4hep::DetType::AUXILIARY  |  DD4hep::DetType::FORWARD ) );

    
    const std::vector<DD4hep::DDRec::LayeredCalorimeterStruct::Layer> & layers = eCalBarrelExtension->layers;
    
    LCCollection * mcColl =0;
    getCollection(mcColl,m_inputMCParticleCollection,evt);
    

    LCCollection * relColl =0;
    getCollection(relColl,"RelationCaloHit",evt);

    LCRelationNavigator * nav = new LCRelationNavigator(relColl);
    
    

    
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
        double hit_z = hit->getPosition()[2];
//         double radius= sqrt(hit_x*hit_x+hit_y*hit_y);
        double hitEnergy = hit->getEnergy();
        lcio::long64 cellId = long( hit->getCellID0() & 0xffffffff ) | ( long( hit->getCellID1() ) << 32 );
        _encoder.setValue(cellId);
        int layer=_encoder["layer"].value();
        
  
    
        m_hit_x->push_back(hit_x);
        m_hit_y->push_back(hit_y);
        m_hit_z->push_back(hit_z);
        
        m_hitEnergies->push_back(hitEnergy);
        m_hitLayers->push_back(layer);
        
        double thickness = (layers[layer].inner_thickness+layers[layer].outer_thickness)/dd4hep::mm;
        double nRadLengths = layers[layer].inner_nRadiationLengths+layers[layer].outer_nRadiationLengths;
        double distance = layers[layer].distance/dd4hep::mm;
        
        double intX0=layers[layer].inner_nRadiationLengths;
        for(int l=layer-1; l>=0; l--)
            intX0=intX0+layers[l].inner_nRadiationLengths+layers[l].outer_nRadiationLengths;
        
        
        
        m_hitLayerDistances->push_back(distance);
        m_hitLayerIntRadiationLengths->push_back(intX0);

        m_hitLayerThicknesses->push_back(thickness);
        m_hitLayerRadiationLengths->push_back(nRadLengths);
        m_nhits++;
        
        
        totalEnergy = totalEnergy + hitEnergy;
        m_showerHist->Fill(distance, trueEnergy,hitEnergy/nRadLengths);
        m_showerHistLayers->Fill(layer, trueEnergy,hitEnergy/nRadLengths);
        m_showerHistX0->Fill(intX0, trueEnergy,hitEnergy/nRadLengths);
        
        const EVENT::LCObjectVec & rawHitVec = nav->getRelatedToObjects(hit);
        
        //There should be only one hit
        float rawHitEnergy = dynamic_cast<SimCalorimeterHit*>(rawHitVec.at(0))->getEnergy();
        
        m_raw_showerHist->Fill(distance, trueEnergy,rawHitEnergy);
        m_raw_showerHistLayers->Fill(layer, trueEnergy,rawHitEnergy);
        m_raw_showerHistX0->Fill(intX0, trueEnergy,rawHitEnergy);
        
        m_raw_hitEnergies->push_back(hitEnergy);

    }
    
    
    m_totalEnergyHist->Fill(totalEnergy,trueEnergy);

    m_totalEnergy= totalEnergy;
    m_trueEnergy = trueEnergy;
    
    caloColl =0;
    
    //All X0 before HCal (due to the ECal)
    double startX0=0.;
    
    for(int l=layers.size()-1; l>=0; l--)
        startX0=startX0+layers[l].inner_nRadiationLengths+layers[l].outer_nRadiationLengths;
    
    int startLayer= layers.size(); //Use largest layer number +1 from ECal
    
    getCollection(caloColl,m_inputLeakageCalorimeterHitCollection,evt);
    UTIL::BitField64 _leak_encoder(caloColl->getParameters().getStringVal( EVENT::LCIO::CellIDEncoding ));
    
    
    //Get HCal Barrel extension by type, ignore plugs and rings 
    const DD4hep::DDRec::LayeredCalorimeterData * hCalBarrelExtension= getExtension( ( DD4hep::DetType::CALORIMETER | DD4hep::DetType::HADRONIC| DD4hep::DetType::BARREL), 
										     ( DD4hep::DetType::AUXILIARY  |  DD4hep::DetType::FORWARD ) );

    
    const std::vector<DD4hep::DDRec::LayeredCalorimeterStruct::Layer> & leak_layers = hCalBarrelExtension->layers;
    double totalLeakEnergy=0.;
    for(int i=0; i< caloColl->getNumberOfElements(); i++){ 
	CalorimeterHit* hit = dynamic_cast<CalorimeterHit*>( caloColl->getElementAt(i) ) ;
        
        double hit_x = hit->getPosition()[0];
        double hit_y = hit->getPosition()[1];
        double hit_z = hit->getPosition()[2];
//         double radius= sqrt(hit_x*hit_x+hit_y*hit_y);
        double hitEnergy = hit->getEnergy();
        lcio::long64 cellId = long( hit->getCellID0() & 0xffffffff ) | ( long( hit->getCellID1() ) << 32 );
        _leak_encoder.setValue(cellId);
        int layer=_leak_encoder["layer"].value();
        
  
    
        m_leak_hit_x->push_back(hit_x);
        m_leak_hit_y->push_back(hit_y);
        m_leak_hit_z->push_back(hit_z);

        m_leak_hitEnergies->push_back(hitEnergy);
        m_leak_hitLayers->push_back(layer);
        
        double thickness = (leak_layers[layer].inner_thickness+leak_layers[layer].outer_thickness)/dd4hep::mm;
        double nRadLengths = leak_layers[layer].inner_nRadiationLengths+leak_layers[layer].outer_nRadiationLengths;
        double distance = leak_layers[layer].distance/dd4hep::mm;
        
        double intX0=leak_layers[layer].inner_nRadiationLengths + startX0;
        for(int l=layer-1; l>=0; l--)
            intX0=intX0+leak_layers[l].inner_nRadiationLengths+leak_layers[l].outer_nRadiationLengths;
        
        
        
        m_leak_hitLayerDistances->push_back(distance);
        m_leak_hitLayerIntRadiationLengths->push_back(intX0);

        m_leak_hitLayerThicknesses->push_back(thickness);
        m_leak_hitLayerRadiationLengths->push_back(nRadLengths);

        m_leak_nhits++;
        
        
        totalLeakEnergy = totalLeakEnergy + hitEnergy;
        m_showerHist->Fill(distance, trueEnergy,hitEnergy/nRadLengths);
        m_showerHistLayers->Fill(startLayer+layer, trueEnergy,hitEnergy/nRadLengths);
        m_showerHistX0->Fill(intX0, trueEnergy,hitEnergy/nRadLengths);
        
        const EVENT::LCObjectVec & rawHitVec = nav->getRelatedToObjects(hit);
        
        //There should be only one hit
        float rawHitEnergy = dynamic_cast<SimCalorimeterHit*>(rawHitVec.at(0))->getEnergy();
        
        m_raw_showerHist->Fill(distance, trueEnergy,rawHitEnergy);
        m_raw_showerHistLayers->Fill(startLayer+layer, trueEnergy,rawHitEnergy);
        m_raw_showerHistX0->Fill(intX0, trueEnergy,rawHitEnergy);
        
        m_leak_raw_hitEnergies->push_back(hitEnergy);

    }
    
    
    m_leakageProfile->Fill(trueEnergy,(totalEnergy+totalLeakEnergy>0?totalLeakEnergy/(totalEnergy+totalLeakEnergy):0));
    m_totalLeakEnergy= totalLeakEnergy;

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
    m_showerHist->Write();
    m_showerHistLayers->Write();
    m_showerHistX0->Write();
    m_leakageProfile->Write();
    m_raw_showerHist->Write();
    m_raw_showerHistLayers->Write();
    m_raw_showerHistX0->Write();
    
    m_outputTree->Write();
    m_rootFile->Close();
}
