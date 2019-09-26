/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "SkimSecondaryVertices.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DD4hep/DetType.h"
#include "DDRec/DetectorData.h"
#include "DD4hep/DetectorSelector.h"

#include "marlin/AIDAProcessor.h"

#include "MarlinTrk/HelixTrack.h"
#include <marlinutil/HelixClass.h>

#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Vertex.h>

#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTrackerConf.h>
#include <UTIL/LCRelationNavigator.h>

#include <marlinutil/GeometryUtil.h>

#include <AIDA/IAnalysisFactory.h>
#include <AIDA/IHistogramFactory.h>

#include <TLorentzVector.h>
#include <TTree.h>


using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace AIDA ;

using dd4hep::DetType;
using dd4hep::rec::ZPlanarData;
using dd4hep::rec::ZDiskPetalsData;

SkimSecondaryVertices aSkimSecondaryVertices;
/*

   Processor to remove secondary vertices from material interactions and validate the skimmed vertices used for Flavour Tagging. Uses BuildUpVertices info.

*/
SkimSecondaryVertices::SkimSecondaryVertices() : Processor("SkimSecondaryVertices") {

  // processor description
  _description = "SkimSecondaryVertices removes secondary vertices from material interactions and creates a skimmed secondary vertices collection to use for the ntuple making for flavour tagging";

  // Input collections
  registerInputCollection( LCIO::VERTEX,
                           "BuildUpVertices",
                           "Name of the Vertex input collection",
                           _inputBuildUpVertices,
                           std::string("BuildUpVertices"));

  // Output collection
  registerOutputCollection( LCIO::VERTEX,
                            "BuildUpVerticesSkimmed",
                            "Name of the Vertex skimmed output collection",
                            _outputBuildUpVerticesSkimmed,
                            std::string("BuildUpVerticesSkimmed"));
 
  // Output ntuple names
  registerProcessorParameter( "vertexTreeName",
                              "Name of the output ntuple",
                              _vertexTreeName,
                              std::string("vertexTree"));
 
  registerProcessorParameter( "vertexSkimmedTreeName",
                              "Name of the output skimmed ntuple",
                              _vertexSkimmedTreeName,
                              std::string("vertexSkimmedTree"));

}

void SkimSecondaryVertices::init(){

  // Print the initial parameters
  printParameters();

  // Reset counters
  _runNumber = 0;
  _eventNumber = 0;

  // Initialise histograms
  AIDAProcessor::histogramFactory(this);

  // Output trees
  _vertexTree = new TTree(_vertexTreeName.c_str(), _vertexTreeName.c_str());
  int bufsize = 32000;
  _vertexTree->Branch("_vtxN",&_vtxN,"_vtxN/I");
  _vertexTree->Branch("_vtxX","std::vector<double >",&_vtxX,bufsize,0);
  _vertexTree->Branch("_vtxY","std::vector<double >",&_vtxY,bufsize,0);
  _vertexTree->Branch("_vtxR","std::vector<double >",&_vtxR,bufsize,0);
  _vertexTree->Branch("_vtxZ","std::vector<double >",&_vtxZ,bufsize,0);

  _vertexSkimmedTree = new TTree(_vertexSkimmedTreeName.c_str(), _vertexSkimmedTreeName.c_str());
  _vertexSkimmedTree->Branch("_vtxskimN",&_vtxskimN,"_vtxskimN/I");
  _vertexSkimmedTree->Branch("_vtxskimR","std::vector<double >",&_vtxskimR,bufsize,0);
  _vertexSkimmedTree->Branch("_vtxskimZ","std::vector<double >",&_vtxskimZ,bufsize,0);
 
  // Get detector data

  dd4hep::Detector &mainDetector = dd4hep::Detector::getInstance();

  // - barrels

  const std::vector<dd4hep::DetElement>& barrelDets = dd4hep::DetectorSelector(mainDetector).detectors((dd4hep::DetType::VERTEX && dd4hep::DetType::BARREL));
 
  streamlog_out( DEBUG0 ) << "Number of vertex barrel dets: " << barrelDets.size() << endl;

  for(unsigned int i=0; i<barrelDets.size(); i++){

    try{
      dd4hep::rec::ZPlanarData* barrels = barrelDets.at(i).extension<dd4hep::rec::ZPlanarData>();
      std::vector<ZPlanarData::LayerLayout> layers = barrels->layers;

      for(unsigned int l=0; l<layers.size(); l++){
        ZPlanarData::LayerLayout layer = layers[l];
        streamlog_out( DEBUG0 ) << "***BARREL LAYER " << l << std::endl;
       
        streamlog_out( DEBUG0 ) << "distanceSupport: " << layer.distanceSupport / dd4hep::mm
                                << "\t thicknessSupport: " << layer.thicknessSupport / dd4hep::mm
                                << "\n offsetSupport: " << layer.offsetSupport / dd4hep::mm
                                << "\t widthSupport: " << layer.widthSupport / dd4hep::mm
                                << "\t zHalfSupport: " << layer.zHalfSupport / dd4hep::mm
                                << "\n distanceSensitive:  " << layer.distanceSensitive / dd4hep::mm
                                << "\t thicknessSensitive: " << layer.thicknessSensitive / dd4hep::mm
                                << "\n offsetSensitive: " << layer.offsetSensitive / dd4hep::mm
                                << "\t widthSensitive: " << layer.widthSensitive / dd4hep::mm
                                << "\t zHalfSensitive: " << layer.zHalfSensitive / dd4hep::mm
                                << std::endl;
        _barrel_radii.push_back(layer.distanceSupport / dd4hep::mm);
        _barrel_halfLength = layer.zHalfSupport / dd4hep::mm ;

      }
    } catch (std::exception &e) {
      if( _barrel_radii.size() != barrelDets.size() )_barrel_radii.clear();
      streamlog_out( DEBUG0 ) << "Caught exception " << e.what() << std::endl;
    }

  }

  // - disks

  const std::vector<dd4hep::DetElement>& endcapDets = dd4hep::DetectorSelector(mainDetector).detectors((dd4hep::DetType::VERTEX && dd4hep::DetType::ENDCAP));

  streamlog_out( DEBUG0 ) << "Number of vertex endcap dets: " << endcapDets.size() << endl;

  for(unsigned int i=0; i<endcapDets.size(); i++){

    try{
      dd4hep::rec::ZDiskPetalsData* endcaps = endcapDets.at(i).extension<dd4hep::rec::ZDiskPetalsData>();
      std::vector<ZDiskPetalsData::LayerLayout> disks = endcaps->layers;

      for(unsigned int d=0; d<disks.size(); d++){
        ZDiskPetalsData::LayerLayout disk = disks[d];
        streamlog_out( DEBUG0 ) << "***ENDCAP LAYER " << d << std::endl;
        streamlog_out( DEBUG0 ) << "zPosition: " << disk.zPosition / dd4hep::mm << std::endl;
        streamlog_out( DEBUG0 ) << "distanceSupport: " << disk.distanceSupport  / dd4hep::mm
                                << "\t thicknessSupport: " << disk.thicknessSupport / dd4hep::mm
                                << "\n zOffsetSupport: " << disk.zOffsetSupport  / dd4hep::mm
                                << "\t widthInnerSupport: " << disk.widthInnerSupport / dd4hep::mm
                                << "\t widthOuterSupport: " << disk.widthOuterSupport / dd4hep::mm
                                <<"\t lengthSupport: " << disk.lengthSupport / dd4hep::mm
                                << "\n distanceSensitive:  " << disk.distanceSensitive / dd4hep::mm
                                << "\t thicknessSensitive: " << disk.thicknessSensitive / dd4hep::mm
                                << "\n zOffsetSensitive: " << disk.zOffsetSensitive / dd4hep::mm
                                << "\t widthInnerSensitive: " << disk.widthInnerSensitive / dd4hep::mm
                                << "\t widthOuterSensitive: " << disk.widthOuterSensitive / dd4hep::mm
                                << "\t lengthSensitive: " << disk.lengthSensitive / dd4hep::mm << std::endl;  

        _endcap_z.push_back(disk.zPosition / dd4hep::mm);

      }

    }
    catch (std::exception &e) {
      if( _endcap_z.size() != endcapDets.size() ) _endcap_z.clear();
      streamlog_out( DEBUG0 ) << "Caught exception " << e.what() << std::endl;
    }

  }

}

void SkimSecondaryVertices::processRunHeader( LCRunHeader* ){
  ++_runNumber;
}

void SkimSecondaryVertices::processEvent( LCEvent* evt ){

  streamlog_out( MESSAGE ) << "Processing event " << _eventNumber << std::endl;

  // Get the collection of Vertices
  LCCollection *vtxCollection = 0;
  getCollection(vtxCollection, _inputBuildUpVertices, evt);
  if(vtxCollection == 0) return;

  int nVtx = vtxCollection->getNumberOfElements();
  _vtxN = nVtx;
  streamlog_out( DEBUG0 ) << nVtx << " vertices in event " << _eventNumber << std::endl;

  // Output vertex collection to fill in
  LCCollectionVec* vtxVec = new LCCollectionVec(LCIO::VERTEX);
  vtxVec->setSubset(true); // flag as subset
  int vtxskimN = 0;

  for(int i=0; i<nVtx; i++){
    Vertex *vtx = dynamic_cast<Vertex*> (vtxCollection->getElementAt(i));
    double vtxX = vtx->getPosition()[0];
    double vtxY = vtx->getPosition()[1];
    double vtxR = sqrt(vtxX*vtxX + vtxY*vtxY);
    double vtxZ = vtx->getPosition()[2];
    _vtxX.push_back(vtxX);
    _vtxY.push_back(vtxY);
    _vtxR.push_back(vtxR);
    _vtxZ.push_back(vtxZ);
    streamlog_out( DEBUG0 ) << "** Vertex " << i << ": x = " << vtxX << "\t y = " << vtxY << "\t z = " << vtxZ << std::endl;

    bool keep = true;
    for(unsigned int l=0; l<_barrel_radii.size(); l++){
      if(l%2==0){
        if(_vtxR.at(i) >= _barrel_radii.at(l) && _vtxR.at(i) < _barrel_radii.at(l+1) && abs(_vtxZ.at(i)) <= _barrel_halfLength ){
          keep = false;
          streamlog_out( DEBUG0 ) << "Vertex skimmed: r = " << _vtxR.at(i) << "\t z = " << _vtxZ.at(i) << std::endl;
        }
      }
    }
    for(unsigned int l=0; l<_endcap_z.size(); l++){
      if(l%2==0){
        if(abs(_vtxZ.at(i)) <= _endcap_z.at(l) && abs(_vtxZ.at(i)) > _endcap_z.at(l+1)){
          keep = false;
          streamlog_out( DEBUG0 ) << "Vertex skimmed: z = " << _vtxZ.at(i) << std::endl;
        }
      }
    }

    if(keep){
      vtxVec->addElement(vtx);
      _vtxskimR.push_back(vtxR);
      _vtxskimZ.push_back(vtxZ);
      vtxskimN++;
    }
  }

  _vtxskimN = vtxskimN;

  evt->addCollection(vtxVec, _outputBuildUpVerticesSkimmed);

  _vertexTree->Fill();
  _vertexSkimmedTree->Fill();
  _eventNumber++;

  // Clear the trees
  clearTreeVar();
}

void SkimSecondaryVertices::check( LCEvent *) {
  //nothing to check at the moment
}

void SkimSecondaryVertices::end(){

  streamlog_out( MESSAGE ) << " end() " << name() << " processed " << _eventNumber << " events in " << _runNumber << " runs " << std::endl;

}

void SkimSecondaryVertices::getCollection(LCCollection*& collection, std::string collectionName, LCEvent* evt){

  try{
    collection = evt->getCollection(collectionName);
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG5 ) << " - cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

void SkimSecondaryVertices::clearTreeVar(){

  _vtxX.clear();
  _vtxY.clear();
  _vtxR.clear();
  _vtxZ.clear();

  _vtxskimR.clear();
  _vtxskimZ.clear();
 
}
