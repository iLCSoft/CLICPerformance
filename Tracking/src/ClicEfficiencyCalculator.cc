/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "ClicEfficiencyCalculator.h"

#include "marlin/AIDAProcessor.h"

#include <EVENT/LCCollection.h>
#include <EVENT/SimTrackerHit.h>

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

ClicEfficiencyCalculator aClicEfficiencyCalculator;

/*
 
 Efficiency calculator for tracking. WARNING: Only works with 1 input collection at present
 
 */

ClicEfficiencyCalculator::ClicEfficiencyCalculator() : Processor("ClicEfficiencyCalculator") {
  
  // modify processor description
  _description = "ClicEfficiencyCalculator calculates the tracking efficiency and makes some performance plots" ;
  
  // Input collections
  registerInputCollection( LCIO::TRACK,
                          "TrackCollectionName",
                          "Track collection name",
                          m_inputTrackCollection,
                          std::string("CATracks"));
  
  registerInputCollection( LCIO::MCPARTICLE,
                          "MCParticleCollectionName",
                          "Name of the MCParticle input collection",
                          m_inputParticleCollection,
                          std::string("MCParticle"));

  registerInputCollection( LCIO::MCPARTICLE,
                           "MCPhysicsParticleCollectionName",
                           "Name of the MCPhysicsParticle input collection",
                           m_inputPhysicsParticleCollection,
                           std::string("MCPhysicsParticle"));
  
  // All tracker hit collections
  std::vector<std::string> inputTrackerHitCollections;
  inputTrackerHitCollections.push_back(std::string("VXDTrackerHits"));
  inputTrackerHitCollections.push_back(std::string("VXDEndcapTrackerHits"));
  inputTrackerHitCollections.push_back(std::string("ITrackerHits"));
  inputTrackerHitCollections.push_back(std::string("OTrackerHits"));
  inputTrackerHitCollections.push_back(std::string("ITrackerEndcapHits"));
  inputTrackerHitCollections.push_back(std::string("OTrackerEndcapHits"));
  
  registerInputCollections( LCIO::TRACKERHITPLANE,
                           "TrackerHitCollectionNames" ,
                           "Name of the TrackerHit input collections"  ,
                           m_inputTrackerHitCollections ,
                           inputTrackerHitCollections ) ;
  
  // All tracker hit relation collections
  std::vector<std::string> inputTrackerHitRelationCollections;
  inputTrackerHitRelationCollections.push_back(std::string("VXDTrackerHitRelations"));
  inputTrackerHitRelationCollections.push_back(std::string("VXDEndcapTrackerHitRelations"));
  inputTrackerHitRelationCollections.push_back(std::string("InnerTrackerBarrelHitsRelations"));
  inputTrackerHitRelationCollections.push_back(std::string("OuterTrackerBarrelHitsRelations"));
  inputTrackerHitRelationCollections.push_back(std::string("InnerTrackerEndcapHitsRelations"));
  inputTrackerHitRelationCollections.push_back(std::string("OuterTrackerEndcapHitsRelations"));
  
  registerInputCollections( LCIO::LCRELATION,
                           "TrackerHitRelCollectionNames",
                           "Name of TrackerHit relation collections",
                           m_inputTrackerHitRelationCollections,
                           inputTrackerHitRelationCollections );
  
  
  // Flag to decide which definition of reconstructable to use
  registerProcessorParameter( "reconstructableDefinition",
                             "Set of cuts to define 'reconstractable' particles for eff computation. The options are: NHits, NHitsVXD, ILDLike",
                             m_cuts,
                             std::string("NHits"));
  
  // Flag for ntuples to be written
  registerProcessorParameter( "debugPlots",
                             "If true additional plots (n of hits per subdetector per mc particle, mc theta, mc pt, info if the particle is decayed in the tracker) will be added to the Ntuple mctree",
                             m_fullOutput,
                             bool(true));

  // Simple ntuple (separate flag)
  registerProcessorParameter( "simpleOutput",
                             "A simplified tree can be output with the MC particle properties, including whether or not it was reconstructed",
                             m_simpleOutput,
                             bool(true));

  // Output collections - particles reconstructed and not reconstructed
  registerOutputCollection( LCIO::MCPARTICLE,
                           "EfficientMCParticleCollectionName",
                           "MC Particle Collection Name for particles reconstructed",
                           m_outputEfficientMCParticleCollection,
                           std::string("EfficientMCParticles"));
  
  registerOutputCollection( LCIO::MCPARTICLE,
                           "InefficientMCParticleCollectionName",
                           "MC Particle Collection Name for particles not reconstructed",
                           m_outputInefficientMCParticleCollection,
                           std::string("InefficientMCParticles"));
  
  // Output ntuple names
  registerProcessorParameter("efficiencyTreeName",
                             "Name of the efficiency tree",
                             m_efficiencyTreeName,
                             std::string("trackTree"));
  
  registerProcessorParameter("purityTreeName",
                             "Name of the purity tree",
                             m_purityTreeName,
                             std::string("purityTree"));
  
  registerProcessorParameter("mcTreeName",
                             "Name of the mc tree",
                             m_mcTreeName,
                             std::string("mctree"));
  
}


void ClicEfficiencyCalculator::init() {
  
  // Print the initial parameters
  printParameters() ;
  
  // Reset counters
  m_runNumber = 0 ;
  m_eventNumber = 0 ;
  
  // Get the magnetic field
  m_magneticField = MarlinUtil::getBzAtOrigin();
  
  // Initialise histograms
  AIDAProcessor::histogramFactory(this);

  
  // Set up trees if ntuple output enabled
  if(m_fullOutput){
    
    // Tree for track information
    m_trackTree = new TTree(m_efficiencyTreeName.c_str(),m_efficiencyTreeName.c_str());
    int bufsize_trk = 32000; //default buffer size 32KB
    m_trackTree->Branch("vx_reconstructable", "std::vector<double >",&m_vec_vx_reconstructable,bufsize_trk,0);
    m_trackTree->Branch("vy_reconstructable", "std::vector<double >",&m_vec_vy_reconstructable,bufsize_trk,0);
    m_trackTree->Branch("vz_reconstructable", "std::vector<double >",&m_vec_vz_reconstructable,bufsize_trk,0);
    m_trackTree->Branch("vr_reconstructable", "std::vector<double >",&m_vec_vr_reconstructable,bufsize_trk,0);
    m_trackTree->Branch("pt_reconstructable", "std::vector<double >",&m_vec_pt_reconstructable,bufsize_trk,0);
    m_trackTree->Branch("theta_reconstructable", "std::vector<double >",&m_vec_theta_reconstructable,bufsize_trk,0);
    m_trackTree->Branch("is_reconstructed", "std::vector<bool >",&m_vec_is_reconstructed,bufsize_trk,0);
    
    // Tree for purity information
    m_purityTree = new TTree(m_purityTreeName.c_str(),m_purityTreeName.c_str());
    int bufsize_purity = 32000; //default buffer size 32KB
    m_purityTree->Branch("trk_nhits_vtx", "std::vector<int >",&m_vec_nhits_vtx,bufsize_purity,0);
    m_purityTree->Branch("trk_nhits_trk", "std::vector<int >",&m_vec_nhits_trk,bufsize_purity,0);
    m_purityTree->Branch("trk_nhits", "std::vector<int >",&m_vec_nhits,bufsize_purity,0);
    m_purityTree->Branch("trk_purity", "std::vector<double >",&m_vec_purity,bufsize_purity,0);
    m_purityTree->Branch("mc_pdg", "std::vector<int >",&m_vec_pdg,bufsize_purity,0);
    m_purityTree->Branch("mc_theta", "std::vector<double >",&m_vec_theta,bufsize_purity,0);
    m_purityTree->Branch("mc_phi", "std::vector<double >",&m_vec_phi,bufsize_purity,0);
    m_purityTree->Branch("mc_p", "std::vector<double >",&m_vec_p,bufsize_purity,0);
    
    // Tree for MC information
    m_mctree = new TTree(m_mcTreeName.c_str(),m_mcTreeName.c_str());
    int bufsize = 32000; //default buffer size 32KB
    //mc categorization: 0 is charge but not reconstructable, 1 is reconstructable but not reconstructed, 2 is reconstructed
    m_mctree->Branch("mcCat", "std::vector<int >",&m_mcCat,bufsize,0);
    m_mctree->Branch("mcTheta", "std::vector<double >",&m_mcTheta,bufsize,0);
    m_mctree->Branch("mcPt", "std::vector<double >",&m_mcPt,bufsize,0);
    m_mctree->Branch("mcNHitsTot", "std::vector<int >",&m_mcNHitsTot,bufsize,0);
    m_mctree->Branch("mcNHitsVXD", "std::vector<int >",&m_mcNHitsVXD,bufsize,0);
    m_mctree->Branch("mcIsDecayedInTracker", "std::vector<int >",&m_mcIsDecayedInTracker,bufsize,0);
    m_mctree->Branch("mcNTrks", "std::vector<int >",&m_mcNTracks,bufsize,0);
    m_mctree->Branch("mcNTrkHits", "std::vector<int >",&m_mcNTrkHits,bufsize,0);
    m_mctree->Branch("mcThetaTrk", "std::vector<int >",&m_mcThetaTrk,bufsize,0);
    m_mctree->Branch("mcPtTrk", "std::vector<int >",&m_mcPtTrk,bufsize,0);
    m_mctree->Branch("mcPhiTrk", "std::vector<int >",&m_mcPhiTrk,bufsize,0);
    m_mctree->Branch("mcNTrksCone", "std::vector<int >",&m_mcNTracksCone,bufsize,0);
    
  }
  
  if(m_simpleOutput){
    m_simplifiedTree = new TTree("simplifiedEfficiencyTree","simplifiedEfficiencyTree");
    m_simplifiedTree->Branch("m_signal",&m_signal,"m_signal/O");
    m_simplifiedTree->Branch("m_type", &m_type, "m_type/D");
    m_simplifiedTree->Branch("m_pt", &m_pt, "m_pt/D");
    m_simplifiedTree->Branch("m_theta", &m_theta, "m_theta/D");
    m_simplifiedTree->Branch("m_phi", &m_phi, "m_phi/D");
    m_simplifiedTree->Branch("m_vertexR", &m_vertexR, "m_vertexR/D");
    m_simplifiedTree->Branch("m_distClosestMCPart",&m_distClosestMCPart,"m_distClosestMCPart/D");
    m_simplifiedTree->Branch("m_closeTracks", &m_closeTracks, "m_closeTracks/D");
    m_simplifiedTree->Branch("m_reconstructed", &m_reconstructed, "m_reconstructed/O");
    m_simplifiedTree->Branch("m_nHits", &m_nHits, "m_nHits/I");
    m_simplifiedTree->Branch("m_nHitsMC", &m_nHitsMC, "m_nHitsMC/I");
    m_simplifiedTree->Branch("m_purity",&m_purity,"m_purity/D");
    m_simplifiedTree->Branch("m_eventNumber", &m_eventNumber, "m_eventNumber/I");
  }
  
  // Set up histograms
  Double_t ptBins[14] = {0.1,0.25,0.5,0.75,1,2.5,5,7.5,10,25,50,75,100,250};
  m_thetaPtMCParticle = new TH2F("m_thetaPtMCParticle","m_thetaPtMCParticle",360,0,180,13,ptBins);
  m_thetaPtMCParticleReconstructed = new TH2F("m_thetaPtMCParticleReconstructed","m_thetaPtMCParticleReconstructed",360,0,180,13,ptBins);
  m_thetaPtMCParticleEfficiency = new TH2F("m_thetaPtMCParticleEfficiency","m_thetaPtMCParticleEfficiency",360,0,180,13,ptBins);
  
}


void ClicEfficiencyCalculator::processRunHeader( LCRunHeader* ) {
  ++m_runNumber ;
}

void ClicEfficiencyCalculator::processEvent( LCEvent* evt ) {
  
  streamlog_out( DEBUG8 ) << "Processing event " << m_eventNumber << std::endl;
  
  // First pick up all of the collections that will be used - tracks, MCparticles, hits from relevent subdetectors - and their relations
  
  // Initialise CELLID encoder to get the subdetector ID from a hit
  UTIL::BitField64 m_encoder( lcio::LCTrackerCellID::encoding_string() ) ;
  std::vector<int> m_collections;
  
  // Get the collection of tracks
  LCCollection* trackCollection = 0 ;
  getCollection(trackCollection, m_inputTrackCollection, evt); if(trackCollection == 0) return;
  
  // Get the collection of MC particles
  LCCollection* particleCollection = 0 ;
  getCollection(particleCollection, m_inputParticleCollection, evt); if(particleCollection == 0) return;

  // Get the collection of MC physics particles
  LCCollection* physicsParticleCollection = 0;
  getCollection(physicsParticleCollection, m_inputPhysicsParticleCollection, evt); 
  if(physicsParticleCollection == 0){
    getCollection(physicsParticleCollection, m_inputParticleCollection, evt);  
  }
  
  // Make the output collections
  LCCollectionVec* efficientParticleCollection = new LCCollectionVec( LCIO::MCPARTICLE )  ;
  efficientParticleCollection->setSubset(true);
  LCCollectionVec* inefficientParticleCollection = new LCCollectionVec( LCIO::MCPARTICLE )  ;
  inefficientParticleCollection->setSubset(true);
  
  // Make objects to hold all of the tracker hit and relation collections
  std::map<int,LCCollection*> trackerHitCollections;
  std::map<int,LCCollection*> trackerHitRelationCollections;
  std::map<int,std::shared_ptr<LCRelationNavigator> > relations;
  std::map<int,bool> activeDetectors;
  
  // Loop over each input collection and get the data
  for(unsigned int collection=0; collection<m_inputTrackerHitCollections.size();collection++){
    
    // Get the collection of tracker hits
    LCCollection* trackerHitCollection = 0 ;
    getCollection(trackerHitCollection, m_inputTrackerHitCollections[collection], evt); if(trackerHitCollection == 0) continue;
    
    // If there are no hits in the hit collection, continue
    if(trackerHitCollection->getNumberOfElements() == 0) continue;
    
    // Get the collection of tracker hit relations
    LCCollection* trackerHitRelationCollection = 0 ;
    getCollection(trackerHitRelationCollection, m_inputTrackerHitRelationCollections[collection], evt);
    
    // Create the relations navigator
    std::shared_ptr<LCRelationNavigator> relation =
      std::make_shared<LCRelationNavigator>( LCRelationNavigator( trackerHitRelationCollection ) );
    
    // Get the subdetector for this set of collections
    int subdetector = getSubdetector(trackerHitCollection,m_encoder);
    activeDetectors[subdetector] = true;
    
    // Save all of the created collections
    trackerHitCollections[subdetector] = trackerHitCollection;
    trackerHitRelationCollections[subdetector] = trackerHitRelationCollection;
    relations[subdetector] = relation;
    m_collections.push_back(subdetector);
    
  }
  
  // Store all physics MCParticles in a std::set
  std::set<MCParticle*> physicsParticles;

  for(int itPhys=0; itPhys<physicsParticleCollection->getNumberOfElements(); itPhys++){
    MCParticle *signal = static_cast<MCParticle*> (physicsParticleCollection->getElementAt(itPhys));
    physicsParticles.insert(signal);
  }

  std::map<MCParticle*,int> particleTracks;
  std::map<MCParticle*,int> particleTrackHits;
  std::list<std::pair<MCParticle*,double>> purityTrack;
  std::list<std::pair<MCParticle*,double>>::const_iterator iter;

  /*
   Look at all tracks that were reconstructed and calculate their purity. This can be used to point to all
	 	reconstructed MC particles. Then loop over MC particles and for each that was not reconstructed check the
	  particle properties. Use these to define efficiency etc. For this we also need the list of hits belonging
	  to each particle.
   */
  
  // Loop over all tracks
  int nTracks = trackCollection->getNumberOfElements();
  for(int itTrack=0;itTrack<nTracks;itTrack++){
    
    // Get the track
    Track* track = static_cast<Track*>( trackCollection->getElementAt(itTrack) ) ;
    
    // Get the hits
    const TrackerHitVec& hitVector = track->getTrackerHits();
    
    // If no hits then ignore this track
    if(hitVector.size() == 0){
      streamlog_out(MESSAGE)<<"Track found with no hits on it"<<std::endl;
      continue;
    }
    
    // Some storage to keep track of which particles have hits on this track
    std::vector<MCParticle*> trackParticles;
    std::map<MCParticle*,int> trackParticleHits;
    
    // Loop over hits and check the MC particle that they correspond to, and which subdetector they are on
    int nHitsVertex(0), nHitsTracker(0), nExcluded(0);
    int nHits = hitVector.size();
    for(int itHit=0;itHit<nHits;itHit++){
      // Get the tracker hit
      TrackerHitPlane* hit = static_cast<TrackerHitPlane*>(hitVector.at(itHit));
      // Get the subdetector number
      int subdetector = getSubdetector(hit,m_encoder);
      // Check if this subdetector is to be included in the efficiency calculation
      if( activeDetectors.count(subdetector) == 0 ){nExcluded++; continue;}
      // FIXME - I am hard coded and horrible!
      if (subdetector==1 || subdetector==2) nHitsVertex++;
      else nHitsTracker++;
      // Get the simulated hit
      const LCObjectVec& simHitVector = relations[subdetector]->getRelatedToObjects( hit );
      // Take the first hit only (this should be changed? Yes - loop over all related simHits and add an entry for each mcparticle so that this hit is in each fit)
      SimTrackerHit* simHit = static_cast<SimTrackerHit*>(simHitVector.at(0));
      // Get the particle belonging to that hit
      MCParticle* particle = simHit->getMCParticle();
      // Check if this particle already has hits on this track
      if( std::find(trackParticles.begin(), trackParticles.end(), particle) == trackParticles.end() ) trackParticles.push_back(particle);
      // Increment the number of counts
      trackParticleHits[particle]++;
    }
    
    // Check the purity of the track, treat ghosts
    int nTrackParticles = trackParticles.size(); int maxHits=0;
    MCParticle* associatedParticle = 0;
    for(int itParticle=0; itParticle<nTrackParticles; itParticle++){
      int nParticleHits = trackParticleHits[trackParticles[itParticle]];
      if(nParticleHits>maxHits){ maxHits=nParticleHits; associatedParticle = trackParticles[itParticle]; }
    }
    double purity = (double)maxHits/(double)(nHits-nExcluded);
    
    purityTrack.push_back(make_pair(associatedParticle, purity));

    // Save additional information
    if(m_fullOutput){
      
      // Number of hits of each type
      m_vec_nhits_vtx.push_back(nHitsVertex);
      m_vec_nhits_trk.push_back(nHitsTracker);
      m_vec_nhits.push_back(nHits-nExcluded);
      m_vec_purity.push_back(purity);
      
      // If no particle associated to this track
      if (associatedParticle == 0){
        m_vec_theta.push_back(-1.);
        m_vec_phi.push_back(0.);
        m_vec_p.push_back(0.);
        m_vec_pdg.push_back(-1);
      }else {
        
        // Store the MC particle information
        TLorentzVector mc_helper;
        mc_helper.SetPxPyPzE(associatedParticle->getMomentum()[0],associatedParticle->getMomentum()[1],associatedParticle->getMomentum()[2],associatedParticle->getEnergy());
        m_vec_theta.push_back(mc_helper.Theta());
        m_vec_phi.push_back(mc_helper.Phi());
        m_vec_p.push_back(mc_helper.P());
        m_vec_pdg.push_back(fabs(associatedParticle->getPDG()));
      }
    }
    
    // Now have a track which is associated to a particle
    particleTracks[associatedParticle]++;
    if(particleTrackHits.count(associatedParticle) == 0 || nHits > particleTrackHits[associatedParticle]) particleTrackHits[associatedParticle] = nHits;
    
  }

  /*
   Now for each MC particle we want the list of hits belonging to it. The most
   efficient way is to loop over all hits once, and store the pointers in a
   map, with the key a pointer to the MC particle. We can then loop over each
   MC particle at the end and get all of the hits, before making a track.
   */
  
  // Loop over all input collections
  for(unsigned int itCollection=0; itCollection<m_collections.size();itCollection++){
    
    // Get the collection ID (subdetector number)
    int collection = m_collections[itCollection];
    
    // Loop over tracker hits
    int nHits = trackerHitCollections[collection]->getNumberOfElements();
    for(int itHit=0;itHit<nHits;itHit++){
      
      // Get the hit
      TrackerHitPlane* hit = static_cast<TrackerHitPlane*>( trackerHitCollections[collection]->getElementAt(itHit) ) ;
      
      // Get the related simulated hit(s)
      const LCObjectVec& simHitVector = relations[collection]->getRelatedToObjects( hit );
      
      // Take the first hit only (this should be changed? Yes - loop over all related simHits and add an entry for each mcparticle so that this hit is in each fit)
      SimTrackerHit* simHit = static_cast<SimTrackerHit*>(simHitVector.at(0));
      
      // Get the particle belonging to that hit
      MCParticle* particle = simHit->getMCParticle();

      // Push back the element into the container
      particleHits[particle].push_back(hit);
      
    }
  }
  
  /*
   Now there is a map containing the number of tracks associated to each MC particle, and a
   map of all hits on each particle. Loop over all MC particles, decide if they are
   reconstructable or not, and check if they were reconstructed more than once (clones).
   */
  
  // Loop over particles
  int nReconstructed(0), nReconstructable(0), nChargePart(0), nCloseTrk(0);
  int nParticles = particleCollection->getNumberOfElements();
  for(int itParticle=0;itParticle<nParticles;itParticle++){
    
    // Get the particle
    nCloseTrk=0;
    MCParticle* particle = static_cast<MCParticle*>( particleCollection->getElementAt(itParticle) ) ;
    
    bool isSignal = physicsParticles.count(particle) == 1;

    // Calculate the particle properties
    TLorentzVector particleVector;
    particleVector.SetPxPyPzE(particle->getMomentum()[0],particle->getMomentum()[1],particle->getMomentum()[2],particle->getEnergy());
    double mcTheta=particleVector.Theta()*180./M_PI;
    double mcPhi=particleVector.Phi()*180./M_PI;
    double mcPt=particleVector.Pt();
    double mcVertexX=particle->getVertex()[0];
    double mcVertexY=particle->getVertex()[1];
    double mcVertexZ=particle->getVertex()[2];
    double mcVertexR=sqrt(pow(mcVertexX,2)+pow(mcVertexY,2));
    double mcCharge=fabs(particle->getCharge());
    bool mcDecayedInTracker=particle->isDecayedInTracker();
    
    // Now decide on criteria for different particle types/classifications
    m_particles["all"]++; // all particles
    
    // No hits in the input collections
    if(particleHits.count(particle) == 0) continue;
    
    // Cut on stable particles
    if(particle->getGeneratorStatus() != 1) continue;
    
    // Now decide on criteria for different particle types/classifications
    m_particles["all"]++; // all particles
    
    // Is this particle reconstructable?
    if (!isReconstructable(particle,m_cuts,m_encoder)) continue;
    m_particles["reconstructable"]++; // reconstructable particles
    nReconstructable++;
    m_thetaPtMCParticle->Fill(mcTheta,mcPt);
    std::vector<TrackerHit*> trackHits = particleHits[particle];
    int uniqueHits = getUniqueHits(trackHits,m_encoder);
    //int nHitsOnTrack(0);

    // Check if it was reconstructed
    bool reconstructed=false;
    if(particleTracks.count(particle)){
      reconstructed=true;
      for(iter = purityTrack.begin(); iter != purityTrack.end(); iter++){
        if((*iter).first == particle){
          m_purity = (*iter).second;
        }
      }
      nReconstructed++;
      m_particles["reconstructed"]++;
      m_thetaPtMCParticleReconstructed->Fill(mcTheta,mcPt);
      efficientParticleCollection->addElement(particle);
      streamlog_out( DEBUG8 ) << "Reconstructed particle with pt: " << mcPt << std::endl;
    }else{
      inefficientParticleCollection->addElement(particle);
        streamlog_out( DEBUG8 ) << "Failed to reconstruct particle with pt: " << mcPt
                                 << " and pdg value " << particle->getPDG()
                                 <<std::endl;
    }
    streamlog_out( DEBUG8 ) << "- particle produced at r = " << mcVertexR << std::endl;

    // Check for particles close to each other (needs cleaning up FIXME)
    double minDR = DBL_MAX;
    for(int j=0; j<nParticles; j++){
      if (itParticle!=j){
        MCParticle* particle2 = static_cast<MCParticle*>( particleCollection->getElementAt(j) );
        bool part2IsCharge = fabs(particle2->getCharge()) > 0.5;
        bool part2IsStable = particle2->getGeneratorStatus() == 1 ;
        if ( part2IsCharge && part2IsStable ) {
          TLorentzVector particleVector2;
          particleVector2.SetPxPyPzE(particle2->getMomentum()[0],particle2->getMomentum()[1],particle2->getMomentum()[2],particle2->getEnergy());
          double DR = particleVector.DeltaR(particleVector2);
          minDR = std::min(DR,minDR);
          if (DR<0.4) nCloseTrk++;
        }
      }
    }

    // Store data for trees
    if (m_fullOutput) {
      m_vec_vx_reconstructable.push_back(mcVertexX);
      m_vec_vy_reconstructable.push_back(mcVertexY);
      m_vec_vz_reconstructable.push_back(mcVertexZ);
      m_vec_vr_reconstructable.push_back(mcVertexR);
      m_vec_pt_reconstructable.push_back(mcPt);
      m_vec_theta_reconstructable.push_back(mcTheta);
      if(reconstructed){
        m_vec_is_reconstructed.push_back(true);
        // m_mcCat can be emtpy at this point
        // m_mcCat.pop_back();
        m_mcCat.push_back(2);
      }else{
        // m_mcCat.pop_back();
        m_mcCat.push_back(1);
      }
  
      // // information for unreconstructed particles
      // if (!reconstructed) {
      //   m_mcNTracks.push_back(nChargePart);
      //   m_mcNTrkHits.push_back(m_mcNHitsTot.back());
      //   m_mcThetaTrk.push_back(m_mcTheta.back());
      //   m_mcPtTrk.push_back(m_mcPt.back());
      //   m_mcNTracksCone.push_back(nCloseTrk);
      // }

      // Only look at charged stable particles. Some descriptive comment should be added
      if (mcCharge>0.5 && particle->getGeneratorStatus()==1){
        nChargePart++;
        m_mcCat.push_back(0);
        m_mcTheta.push_back(mcTheta);
        m_mcPt.push_back(mcPt);
        m_mcIsDecayedInTracker.push_back(mcDecayedInTracker);
        std::vector<TrackerHit*> trackHits_helper = particleHits[particle];
        int nhitsVXD = 0;
        for (size_t ihit=0; ihit<trackHits_helper.size(); ihit++){
          int subdetector = getSubdetector(trackHits_helper.at(ihit), m_encoder);
          if (subdetector==1 || subdetector==2) nhitsVXD++;
        }
        m_mcNHitsVXD.push_back(nhitsVXD);
        m_mcNHitsTot.push_back(trackHits_helper.size());
      }
    }
    
    if(m_simpleOutput){
      m_signal = isSignal;
      m_type = particle->getPDG();
      m_pt = mcPt;
      m_theta = mcTheta;
      m_phi = mcPhi;
      m_vertexR = mcVertexR;
      m_reconstructed = reconstructed;
      m_nHits = particleTrackHits[particle];
      m_nHitsMC = uniqueHits;
      m_distClosestMCPart = minDR;
      m_closeTracks = nCloseTrk;
      m_simplifiedTree->Fill();
    }
  }
  
  // Fill the output trees
  if (m_fullOutput) {
    m_mctree->Fill();
    m_trackTree->Fill();
    m_purityTree->Fill();
  }
  
  // Store the reconstructed and unreconstructed particle collections
  evt->addCollection( efficientParticleCollection , m_outputEfficientMCParticleCollection ) ;
  evt->addCollection( inefficientParticleCollection , m_outputInefficientMCParticleCollection ) ;
  
  // Increment the event number
  streamlog_out( DEBUG9 ) <<"For this event reconstructed "
                          << 100.*(double)nReconstructed/(double)nReconstructable
                          <<" % (" << nReconstructed << "/" << nReconstructable <<")"
                          << std::endl;
  m_eventNumber++ ;

  // Clean the trees
  clearTreeVar();

}

void ClicEfficiencyCalculator::check( LCEvent *  ) {
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ClicEfficiencyCalculator::end(){
  
  streamlog_out(MESSAGE) << " end()  " << name()
  << " processed " << m_eventNumber << " events in " << m_runNumber << " runs "
  << std::endl ;
  
  // Calculate efficiency results
  streamlog_out(MESSAGE)<<std::fixed<<std::setprecision(2)<<"Reconstructable particle efficiency: "<<100.*m_particles["reconstructed"]/m_particles["reconstructable"]<<" % ("<<std::setprecision(0)<<m_particles["reconstructed"]<<"/"<<m_particles["reconstructable"]<<")"<<std::endl;
  
  
}

void ClicEfficiencyCalculator::getCollection(LCCollection* &collection, std::string collectionName, LCEvent* evt){
  try{
    collection = evt->getCollection( collectionName ) ;
  }
  catch(DataNotAvailableException &e){
    streamlog_out( DEBUG5 )<<"- cannot get collection. Collection " << collectionName.c_str() << " is unavailable" << std::endl;
    return;
  }
  return;
}

int ClicEfficiencyCalculator::getSubdetector(LCCollection* collection, UTIL::BitField64 &encoder){
  TrackerHitPlane* hit = static_cast<TrackerHitPlane*>( collection->getElementAt(0) ) ;
  const int celId = hit->getCellID0() ;
  encoder.setValue(celId) ;
  int subdet = encoder[lcio::LCTrackerCellID::subdet()];
  return subdet;
}


int ClicEfficiencyCalculator::getSubdetector(TrackerHit* hit, UTIL::BitField64 &encoder){
  const int celId = hit->getCellID0() ;
  encoder.setValue(celId) ;
  int subdet = encoder[lcio::LCTrackerCellID::subdet()];
  return subdet;
}

int ClicEfficiencyCalculator::getLayer(TrackerHit* hit, UTIL::BitField64 &encoder){
  const int celId = hit->getCellID0() ;
  encoder.setValue(celId) ;
  int layer = encoder[lcio::LCTrackerCellID::layer()];
  return layer;
}



bool ClicEfficiencyCalculator::isReconstructable(MCParticle*& particle, std::string cut, UTIL::BitField64 &m_encoder){
  
  
  if (cut=="NHits") {
    
    // Only make tracks with 6 or more hits
    std::vector<TrackerHit*> trackHits = particleHits[particle];
    int uniqueHits = getUniqueHits(trackHits,m_encoder);
    if(uniqueHits >= 4) return true;
    
  } else if (cut=="NHitsVXD") {
    
    // Only make tracks with 4 or more hits in the vertex detector
    std::vector<TrackerHit*> trackHits = particleHits[particle];
    int nVXDHits = 0;
    UTIL::BitField64 encoder( lcio::LCTrackerCellID::encoding_string() ) ;
    for (size_t ihit=0; ihit<trackHits.size(); ihit++){
      int subdetector = getSubdetector(trackHits.at(ihit), encoder);
      if (subdetector==1 || subdetector==2) nVXDHits++;
    }
    if (nVXDHits >= 4) return true;
    
  } else if (cut=="ILDLike") {
    
    
    // Only consider particles: charged, stable, pT>0.1GeV, cosTheta<0.89, nHits>=4, IP in 10 cm
    //(a.t.m. cut in cosTheta is up to 0.89 instead of the usal 0.99 for cutting also particles in the vertex endcap region)
    //(a.t.m. no same cut on reco pT/theta but no bias because in that case also the mc particle is kept)
    
    //add condition for decay outside vertex/tracker?
    //bool isCharge = false;
    bool isStable = false;
    bool passPt = false;
    bool passTheta = false;
    bool passNHits = false;
    //bool passIP = false;
    //bool passEndPoint = false;
    // bool isOverlay = particle->isOverlay();
    // if (isOverlay) std::cout<<"----- ehi there are overlay particles in this mc collection"<<std::endl;

    // //isCharge Cut
    // double charge = fabs(particle->getCharge());
    // if (charge>0.5) isCharge = true;
    //else std::cout<<"----- mc not charged"<<std::endl;
    int genStatus = particle->getGeneratorStatus();
    //int nDaughters = particle->getDaughters().size();
    if (genStatus == 1 ) isStable = true;
    //else std::cout<<"----- mc not stable in generator"<<std::endl;
    
    
    TLorentzVector p;
    p.SetPxPyPzE(particle->getMomentum()[0], particle->getMomentum()[1], particle->getMomentum()[2], particle->getEnergy());//in GeV
    if ( p.Pt()>=0.1 ) passPt = true;
    if ( fabs(cos(p.Theta()))<0.99 ) passTheta = true;
    //if ( fabs(cos(p.Theta()))>0.86 && fabs(cos(p.Theta()))<0.99 ) passTheta = true;
    //if ( fabs(cos(p.Theta()))<0.89 ) passTheta = true;
    //else std::cout<<"----- mc does not pass theta acceptance cut"<<std::endl;
    std::vector<TrackerHit*> trackHits = particleHits[particle];
    int uniqueHits = getUniqueHits(trackHits,m_encoder);
    if(uniqueHits >= 4) passNHits = true;
    //else std::cout<<"----- mc has not 4 hits associated"<<std::endl;
    //if(trackHits.size() >= 6) passNHits = true;
    
    // int nVXDHits = 0;
    // UTIL::BitField64 encoder( lcio::LCTrackerCellID::encoding_string() ) ;
    // for (size_t ihit=0; ihit<trackHits.size(); ihit++){
    //   int subdetector = getSubdetector(trackHits.at(ihit), encoder);
    //   if (subdetector == m_vertexBarrelID) nVXDHits++;
    // }
    // if (nVXDHits >= 4) return true;
    
    // //passIP Cut
    // double dist = sqrt( pow(particle->getVertex()[0],2) + pow(particle->getVertex()[1],2) );
    // if (dist<100.) passIP = true;

    //else std::cout<<"----- mc has dist from IP > 100mm"<<std::endl;
    //if (dist<30.) passIP = true;

    // //passEndpoint Cut
    // double e = sqrt( pow(particle->getEndpoint()[0],2) + pow(particle->getEndpoint()[1],2) );
    // if (e==0. || e>40.) passEndPoint=true;

    //else std::cout<<"----- mc has endpoint < 40mm"<<std::endl;
    
    
    //bool keepParticle = isCharge && isStable && passPt && passTheta && passNHits && passIP && passEndPoint;
    //bool keepParticle = isCharge && passTheta && passNHits;
    //bool keepParticle = passTheta && passNHits && passPt && passIP && passEndPoint && isStable; //81.98%
    //bool keepParticle = passTheta && passNHits && passPt && passIP && isStable; //81.98%
    //bool keepParticle = passTheta && passNHits && passPt && isStable; //81.98%
    //bool keepParticle = passTheta && passNHits && passPt; //59.09
    //bool keepParticle = passTheta && passNHits && passPt && passIP && passEndPoint; //71.65%
    
    bool keepParticle = passTheta && passNHits && passPt && isStable; // && !isOverlay;
    
    
    if (keepParticle) return true;
    
  }  else if (cut=="SingleMu") {
    
    // Only consider particles: muons (|PDG|=13) with at least 4 hits pass
    
    bool isMu = false;
    bool passNHits = false;
    bool passPt = false;
    bool passTheta = false;
    bool passIP = false;
    bool passEndPoint = false;
    //bool passVertexHits = false;
    bool isNotLoop = true;
    
    double pdg = fabs(particle->getPDG());
    if (pdg==13) isMu = true;
    
    std::vector<TrackerHit*> trackHits = particleHits[particle];
    if(trackHits.size() >= 4) passNHits = true;
    
    int nVXDHits = 0;
    UTIL::BitField64 encoder( lcio::LCTrackerCellID::encoding_string() ) ;
    std::vector<int > vec_hit_subdet;
    std::vector<int > vec_hit_layer;
    for (size_t ihit=0; ihit<trackHits.size(); ihit++){
      int subdetector = getSubdetector(trackHits.at(ihit), encoder);
      int layer = getLayer(trackHits.at(ihit), encoder);
      if (subdetector==1 || subdetector==2) nVXDHits++;
      vec_hit_subdet.push_back(subdetector);
      vec_hit_layer.push_back(layer);
    }
    //if (nVXDHits >= 6) passVertexHits = true;
    
    for (size_t j = 0; j < vec_hit_subdet.size()-1; j++){
      for (size_t k = j+1; k < vec_hit_subdet.size(); k++){
        if ( vec_hit_subdet.at(j) == vec_hit_subdet.at(k) ) {
          if ( vec_hit_layer.at(j) == vec_hit_layer.at(k) ) {
            isNotLoop = false;
            streamlog_out( MESSAGE ) << "=====> loop track " << std::endl;
          }
        }
      }
    }
    
    vec_hit_subdet.clear();
    vec_hit_layer.clear();
    
    double dist = sqrt( pow(particle->getVertex()[0],2) + pow(particle->getVertex()[1],2) );
    streamlog_out( DEBUG8 ) << "----- for mc_pdg = " << pdg << "   dist from IP = " << dist << std::endl;
    if (dist<100.) passIP = true;
    
    double e = sqrt( pow(particle->getEndpoint()[0],2) + pow(particle->getEndpoint()[1],2) );
    if (e==0. || e>40.) passEndPoint=true;
    
    TLorentzVector p;
    p.SetPxPyPzE(particle->getMomentum()[0], particle->getMomentum()[1], particle->getMomentum()[2], particle->getEnergy());//in GeV
    if ( p.Pt()>=0.1 ) passPt = true;
    if ( fabs(cos(p.Theta()))<0.99 ) passTheta = true;
    
    //bool keepParticle = isMu && passNHits && passTheta && passPt && passIP && passEndPoint && passVertexHits && isNotLoop;
    bool keepParticle = isMu && passNHits && passTheta && passPt && passIP && passEndPoint && isNotLoop;
    if (keepParticle) return true;
    
  } else {
    streamlog_out( ERROR )<<"Set of cuts " << cut.c_str() << " not defined" << std::endl;
    return false;
  }
  
  return false;
  
}

int ClicEfficiencyCalculator::getUniqueHits(std::vector<TrackerHit*> trackHits, UTIL::BitField64 &encoder){
  
  std::vector<int> uniqueIds;
  
  // Loop over each hit
  for(size_t iHit=0;iHit<trackHits.size();iHit++){
    
    // Get the hit and the cell ID
    TrackerHit* hit = trackHits[iHit];
    const int celId = hit->getCellID0();
    
    // Get the subdetector information
    encoder.setValue(celId) ;
    int subdet = encoder[lcio::LCTrackerCellID::subdet()];
    int side = encoder[lcio::LCTrackerCellID::side()];
    int layer = encoder[lcio::LCTrackerCellID::layer()];
    
    // Make a unique id
    int id = ((subdet & 0xFF) << 16) + ((side & 0xFF) << 8) + (layer & 0xFF);
    if( std::find(uniqueIds.begin(), uniqueIds.end(), id) == uniqueIds.end()) uniqueIds.push_back(id);
  }
  
  return uniqueIds.size();
  
}


void ClicEfficiencyCalculator::clearTreeVar(){
  
  m_mcCat.clear();
  m_mcPt.clear();
  m_mcTheta.clear();
  m_mcIsDecayedInTracker.clear();
  m_mcNHitsTot.clear();
  m_mcNHitsVXD.clear();
  
  m_mcNTracks.clear();
  m_mcNTrkHits.clear();
  m_mcThetaTrk.clear();
  m_mcPtTrk.clear();
  m_mcPhiTrk.clear();
  m_mcNTracksCone.clear();
  
  m_vec_vx_reconstructable.clear();
  m_vec_vy_reconstructable.clear(); 
  m_vec_vz_reconstructable.clear(); 
  m_vec_vr_reconstructable.clear(); 
  m_vec_pt_reconstructable.clear(); 
  m_vec_theta_reconstructable.clear(); 
  m_vec_is_reconstructed.clear(); 
  
  m_vec_nhits_vtx.clear();
  m_vec_nhits_trk.clear();
  m_vec_nhits.clear();
  m_vec_purity.clear();
  m_vec_pdg.clear();
  m_vec_theta.clear();
  m_vec_phi.clear();
  m_vec_p.clear();
  
  particleHits.clear();

}  





