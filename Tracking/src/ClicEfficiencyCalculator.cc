/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
#include "ClicEfficiencyCalculator.h"
#include "DDRec/API/IDDecoder.h"

#include "MarlinTrk/MarlinTrkUtils.h"
#include "MarlinTrk/HelixTrack.h"
#include "MarlinTrk/HelixFit.h"
#include "MarlinTrk/IMarlinTrack.h"
#include "MarlinTrk/Factory.h"
#include "MarlinTrk/MarlinTrkDiagnostics.h"
#include "MarlinTrk/IMarlinTrkSystem.h"

#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <IMPL/TrackImpl.h>

#include <UTIL/CellIDEncoder.h>
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>
#include <UTIL/LCRelationNavigator.h>

#include "DD4hep/LCDD.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/SurfaceManager.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "marlin/ProcessorEventSeeder.h"
#include "marlin/Global.h"

#include "CLHEP/Vector/TwoVector.h"

#include <AIDA/IAnalysisFactory.h>
#include <AIDA/IHistogramFactory.h>

#include <TLorentzVector.h>
#include "TFile.h"
#include "TTree.h"
#include "LinkDef.h"

#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <climits>
#include <cfloat>
#include "sys/types.h"
#include "sys/sysinfo.h"


using namespace lcio ;
using namespace marlin ;
using namespace std ;
using namespace DD4hep ;
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
	
//	registerInputCollection( LCIO::LCRELATION,
//													"TrackRelationCollectionName",
//													"Track relation collection name",
//													m_inputTrackRelationCollection,
//													std::string("SiTrackRelations"));
	
	registerInputCollection( LCIO::MCPARTICLE,
													"MCParticleCollectionName",
													"Name of the MCParticle input collection",
													m_inputParticleCollection,
													std::string("MCParticle"));
	
	// All tracker hit collections
	std::vector<std::string> inputTrackerHitCollections;
	inputTrackerHitCollections.push_back(std::string("VXDTrackerHits"));
	
	registerInputCollections( LCIO::TRACKERHITPLANE,
													 "TrackerHitCollectionNames" ,
													 "Name of the TrackerHit input collections"  ,
													 m_inputTrackerHitCollections ,
													 inputTrackerHitCollections ) ;
	
	// All tracker hit relation collections
	std::vector<std::string> inputTrackerHitRelationCollections;
	inputTrackerHitRelationCollections.push_back(std::string("VXDTrackerHitRelations"));
	
	registerInputCollections( LCIO::LCRELATION,
													 "TrackerHitRelCollectionNames",
													 "Name of TrackerHit relation collections",
													 m_inputTrackerHitRelationCollections,
													 inputTrackerHitRelationCollections );


  registerProcessorParameter( "reconstructableDefinition",
                              "Set of cuts to define 'reconstractable' particles for eff computation. The options are: NHits, NHitsVXD, ILDLike",
                              m_cuts,
                              std::string("NHits"));

  registerProcessorParameter( "vertexBarrelID",
                              "Detector element ID for the vertex Barrel",
                              m_vertexBarrelID,
                              int(1));

  registerProcessorParameter( "morePlots",
                              "If true additional plots (n of hits per subdetector per mc particle, mc theta, mc pt, info if the particle is decayed in the tracker) will be added to the Ntuple mctree",
                              m_morePlots,
                              bool(false));


  
  registerOutputCollection( LCIO::MCPARTICLE,
                            "MCParticleNotReco" , 
                            "Name of the output collection of the not reconstructed charged MCParticle"  ,
                            m_notRecoMCColName ,
                            std::string("MCParticleNotReco") ) ;



  // registerProcessorParameter("outFileName",
  //                            "Name of the output root file",
  //                            _outFileName,
  //                            std::string("out_eff.root")
  //                            );

  registerProcessorParameter("efficiencyTreeName",
                             "Name of the efficiency tree",
                             _effTreeName,
                             std::string("trktree")
                             );

  registerProcessorParameter("purityTreeName",
                             "Name of the purity tree",
                             _purityTreeName,
                             std::string("puritytree")
                             );


  registerProcessorParameter("mcTreeName",
                             "Name of the efficiency tree",
                             _mcTreeName,
                             std::string("mctree")
                             );
	
}


void ClicEfficiencyCalculator::init() {
	
	// Print the initial parameters
	printParameters() ;
	
	// Reset counters
	m_runNumber = 0 ;
	m_eventNumber = 0 ;
	
	// Get the magnetic field
	DD4hep::Geometry::LCDD& lcdd = DD4hep::Geometry::LCDD::getInstance();
	const double position[3]={0,0,0}; // position to calculate magnetic field at (the origin in this case)
	double magneticFieldVector[3]={0,0,0}; // initialise object to hold magnetic field
	lcdd.field().magneticField(position,magneticFieldVector); // get the magnetic field vector from DD4hep
	m_magneticField = magneticFieldVector[2]/dd4hep::tesla; // z component at (0,0,0)
	
	// Define the purity cut - should be made configurable
	m_purity = 0.5;
	
	// Register this process
	Global::EVENTSEEDER->registerProcessor(this);

	


  // _outRootFile = new TFile(_outFileName.c_str(),"RECREATE");
  


  // // Plots

  // eff_vs_theta = new TCanvas("eff_vs_theta","Trk Eff vs Theta",800,800);
  // g_eff_vs_theta = new TGraphAsymmErrors() ;

  // // const int nbins_theta = 13;
  // // double theta_edges[nbins_theta+1] = { 0.0, 10., 20., 30., 45., 60., 70, 87., 93., 110., 145., 160., 170., 180. } ;
  // const int nbins_theta = 18;
  // double theta_edges[nbins_theta+1] = { 0.0, 7., 15., 23., 30., 45., 60, 75, 88., 90., 92., 105., 120., 135., 150., 157., 165., 173., 180.} ;

  // h_theta_reconstructed  = new TH1F( "h_theta_reconstructed", "Theta distributions of reconstructed tracks passing purity criteria", nbins_theta , theta_edges ) ;
  // h_theta_reconstructable  = new TH1F( "h_theta_reconstructable", "Theta distribution of reconstructable tracks", nbins_theta , theta_edges ) ;



  // eff_vs_pt = new TCanvas("eff_vs_pt","Trk Eff vs Pt",800,800);
  // g_eff_vs_pt = new TGraphAsymmErrors() ;

  // // const int nbins_pt = 12;
  // // double pt_edges[nbins_pt+1] = { 0.01, 0.1, 0.5, 1., 2., 5., 10., 20., 50., 100., 200., 500., 1000. } ;
  // const int nbins_pt = 16;
  // double pt_edges[nbins_pt+1] = { 0.01, 0.1, 0.2, 0.5, 1., 2., 4., 8., 10., 20., 40., 80., 100., 150., 200., 500., 1000. } ;

  // h_pt_reconstructed  = new TH1F( "h_pt_reconstructed", "Pt distributions of reconstructed tracks passing purity criteria", nbins_pt , pt_edges ) ;
  // h_pt_reconstructable  = new TH1F( "h_pt_reconstructable", "Pt distribution of reconstructable tracks", nbins_pt , pt_edges ) ;


   
  trktree = new TTree(_effTreeName.c_str(),_effTreeName.c_str());
  int bufsize_trk = 32000; //default buffer size 32KB
  trktree->Branch("vx_reconstructable", "std::vector<double >",&m_vec_vx_reconstructable,bufsize_trk,0); 
  trktree->Branch("vy_reconstructable", "std::vector<double >",&m_vec_vy_reconstructable,bufsize_trk,0); 
  trktree->Branch("vz_reconstructable", "std::vector<double >",&m_vec_vz_reconstructable,bufsize_trk,0); 
  trktree->Branch("vr_reconstructable", "std::vector<double >",&m_vec_vr_reconstructable,bufsize_trk,0); 
  trktree->Branch("pt_reconstructable", "std::vector<double >",&m_vec_pt_reconstructable,bufsize_trk,0); 
  trktree->Branch("theta_reconstructable", "std::vector<double >",&m_vec_theta_reconstructable,bufsize_trk,0); 
  trktree->Branch("is_reconstructed", "std::vector<bool >",&m_vec_is_reconstructed,bufsize_trk,0); 

   
  puritytree = new TTree(_purityTreeName.c_str(),_purityTreeName.c_str());
  int bufsize_purity = 32000; //default buffer size 32KB
  puritytree->Branch("trk_nhits_vtx", "std::vector<int >",&m_vec_nhits_vtx,bufsize_purity,0); 
  puritytree->Branch("trk_nhits_trk", "std::vector<int >",&m_vec_nhits_trk,bufsize_purity,0); 
  puritytree->Branch("trk_nhits", "std::vector<int >",&m_vec_nhits,bufsize_purity,0); 
  puritytree->Branch("trk_purity", "std::vector<double >",&m_vec_purity,bufsize_purity,0); 
  puritytree->Branch("mc_pdg", "std::vector<int >",&m_vec_pdg,bufsize_purity,0); 
  puritytree->Branch("mc_theta", "std::vector<double >",&m_vec_theta,bufsize_purity,0); 
  puritytree->Branch("mc_phi", "std::vector<double >",&m_vec_phi,bufsize_purity,0); 
  puritytree->Branch("mc_p", "std::vector<double >",&m_vec_p,bufsize_purity,0); 



  if (m_morePlots){
    // Tree
      
    mctree = new TTree(_mcTreeName.c_str(),_mcTreeName.c_str());
    int bufsize = 32000; //default buffer size 32KB

    mctree->Branch("mcCat", "std::vector<int >",&m_mcCat,bufsize,0); //mc particle categorization: 0 is charge but not reconstructable, 1 is reconstructable but not reconstructed, 2 is reconstructed
    mctree->Branch("mcTheta", "std::vector<double >",&m_mcTheta,bufsize,0); 
    mctree->Branch("mcPt", "std::vector<double >",&m_mcPt,bufsize,0); 
    //mctree->Branch("mcNHits", "std::vector<std::vector<int > >",&m_mcNHits,bufsize,0); 
    mctree->Branch("mcNHitsTot", "std::vector<int >",&m_mcNHitsTot,bufsize,0); 
    mctree->Branch("mcNHitsVXD", "std::vector<int >",&m_mcNHitsVXD,bufsize,0); 
    mctree->Branch("mcIsDecayedInTracker", "std::vector<int >",&m_mcIsDecayedInTracker,bufsize,0); 

    // not reconstructed tracks
    mctree->Branch("mcNTrks", "std::vector<int >",&m_mcNTracks,bufsize,0); 
    mctree->Branch("mcNTrkHits", "std::vector<int >",&m_mcNTrkHits,bufsize,0); 
    mctree->Branch("mcThetaTrk", "std::vector<int >",&m_mcThetaTrk,bufsize,0); 
    mctree->Branch("mcPtTrk", "std::vector<int >",&m_mcPtTrk,bufsize,0); 
    mctree->Branch("mcPhiTrk", "std::vector<int >",&m_mcPhiTrk,bufsize,0); 
    mctree->Branch("mcNTrksCone", "std::vector<int >",&m_mcNTracksCone,bufsize,0); 

  }




}


void ClicEfficiencyCalculator::processRunHeader( LCRunHeader* run) {
	++m_runNumber ;
}

void ClicEfficiencyCalculator::processEvent( LCEvent* evt ) {

  // if( isFirstEvent() ) { 
  // }

  
	std::cout<<"Processing event "<<m_eventNumber<<std::endl;
	// First pick up all of the collections that will be used - tracks, MCparticles, hits from relevent subdetectors - and their relations
	
  clearTreeVar();

  //LCCollectionVec* skimVec = new LCCollectionVec( LCIO::MCPARTICLE )  ;
  //skimVec->setSubset(true) ; 


	// Initialise CELLID encoder to get the subdetector ID from a hit
	UTIL::BitField64 m_encoder( lcio::ILDCellID0::encoder_string ) ;
	std::vector<int> m_collections;

	// Get the collection of tracks
	LCCollection* trackCollection = 0 ;
	getCollection(trackCollection, m_inputTrackCollection, evt); if(trackCollection == 0) return;
	
	// Get the collection of MC particles
	LCCollection* particleCollection = 0 ;
	getCollection(particleCollection, m_inputParticleCollection, evt); if(particleCollection == 0) return;
	
	// Get the collection of track relations to MC particles
  //	LCCollection* trackRelationCollection = 0 ;
  //	getCollection(trackRelationCollection, m_inputTrackRelationCollection, evt); if(trackRelationCollection == 0) return;
	
	// Create the relations navigator
  //	LCRelationNavigator* trackRelation = new LCRelationNavigator( trackRelationCollection );
	
	// Make objects to hold all of the tracker hit and relation collections
	std::map<int,LCCollection*> trackerHitCollections;
	std::map<int,LCCollection*> trackerHitRelationCollections;
	std::map<int,LCRelationNavigator*> relations;
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
		LCRelationNavigator* relation = new LCRelationNavigator( trackerHitRelationCollection );
		
		// Get the subdetector for this set of collections
		int subdetector = getSubdetector(trackerHitCollection,m_encoder);
		activeDetectors[subdetector] = true;

		// Save all of the created collections
		trackerHitCollections[subdetector] = trackerHitCollection;
		trackerHitRelationCollections[subdetector] = trackerHitRelationCollection;
		relations[subdetector] = relation;
		m_collections.push_back(subdetector);
  
	}



	std::map<MCParticle*,int> particleTracks;
	
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
		Track* track = dynamic_cast<Track*>( trackCollection->getElementAt(itTrack) ) ;
		
		// Get the hits
		const TrackerHitVec& hitVector = track->getTrackerHits();
		
		// Some storage to keep track of which particles have hits on this track
		std::vector<MCParticle*> trackParticles;
		std::map<MCParticle*,int> trackParticleHits;
		
		// Loop over hits and check the MC particle that they correspond to, and which subdetector they are on
    int nHitsVertex = 0;
    int nHitsTracker = 0;
		int nHits = hitVector.size();
		int nExcluded = 0;
		for(int itHit=0;itHit<nHits;itHit++){
			// Get the tracker hit
			TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>(hitVector.at(itHit));
			// Get the subdetector number
			int subdetector = getSubdetector(hit,m_encoder);
			// Check if this subdetector is to be included in the efficiency calculation
			if( activeDetectors.count(subdetector) == 0 ){nExcluded++; continue;}
      if (subdetector==1 || subdetector==2) nHitsVertex++;
      else nHitsTracker++;
			// Get the simulated hit
			const LCObjectVec& simHitVector = relations[subdetector]->getRelatedToObjects( hit );
			// Take the first hit only (this should be changed? Yes - loop over all related simHits and add an entry for each mcparticle so that this hit is in each fit)
			SimTrackerHit* simHit = dynamic_cast<SimTrackerHit*>(simHitVector.at(0));
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

    /*
    bool found = false;
    //check if MC particle associated to the track (from the simhit) is from te physics signal mc collection
    for(int imcp=0; imcp<particleCollection->getNumberOfElements(); imcp++){
      MCParticle* mcp = dynamic_cast<MCParticle*>(particleCollection->getElementAt(imcp));
      if (associatedParticle == mcp){
        found = true;
        break;
      }
    } // end loop on signal physics mc particle collection
    if (!found) {
      //streamlog_out(MESSAGE) << " I AM DOING SOMETHING ABOUT MC PARTICLE NOT FROM THE PHYSICS COLLECTION " << std::endl;
      associatedParticle = 0;
    }
    */


		double purity = (double)maxHits/(double)(nHits-nExcluded);
    streamlog_out(MESSAGE) << " purity = " << purity << std::endl;
    
    m_vec_nhits_vtx.push_back(nHitsVertex);
    m_vec_nhits_trk.push_back(nHitsTracker);
    m_vec_nhits.push_back(nHits-nExcluded);
    m_vec_purity.push_back(purity);

    if (associatedParticle == 0){

      m_vec_theta.push_back(-1.);
      m_vec_phi.push_back(0.);
      m_vec_p.push_back(0.);
      m_vec_pdg.push_back(-1);

    }else {

      TLorentzVector mc_helper;
      mc_helper.SetPxPyPzE(associatedParticle->getMomentum()[0],associatedParticle->getMomentum()[1],associatedParticle->getMomentum()[2],associatedParticle->getEnergy());

      m_vec_theta.push_back(mc_helper.Theta());
      m_vec_phi.push_back(mc_helper.Phi());
      m_vec_p.push_back(mc_helper.P());
      m_vec_pdg.push_back(fabs(associatedParticle->getPDG()));


    }


		if(purity < m_purity) continue; // do something with ghosts here

		// Now have a track which is associated to a particle
		particleTracks[associatedParticle]++;


    // // Fill the histogram of the theta distribution for the reconstructed tracks (numerator for the eff)
    // double tanLambda = track->getTanLambda();
    // double theta = M_PI/2 - atan(tanLambda);
    // h_theta_reconstructed->Fill(theta*180./M_PI);

	}
	
	
	/*
	 Now for each MC particle we want the list of hits belonging to it. The most
	 efficient way is to loop over all hits once, and store the pointers in a
	 map, with the key a pointer to the MC particle. We can then loop over each
	 MC particle at the end and get all of the hits, before making a track.
  */
	
	// Make the container
	//std::map<MCParticle*, std::vector<TrackerHit*> > particleHits;
  particleHits.clear();

	// Loop over all input collections
	for(unsigned int itCollection=0; itCollection<m_collections.size();itCollection++){
	
		// Get the collection ID (subdetector number)
		int collection = m_collections[itCollection];

		// Loop over tracker hits
		int nHits = trackerHitCollections[collection]->getNumberOfElements();
		for(int itHit=0;itHit<nHits;itHit++){
			
			// Get the hit
			TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>( trackerHitCollections[collection]->getElementAt(itHit) ) ;
			
			// Get the related simulated hit(s)
			const LCObjectVec& simHitVector = relations[collection]->getRelatedToObjects( hit );
			
			// Take the first hit only (this should be changed? Yes - loop over all related simHits and add an entry for each mcparticle so that this hit is in each fit)
			SimTrackerHit* simHit = dynamic_cast<SimTrackerHit*>(simHitVector.at(0));
			
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
	
	int nReconstructed=0;
	int nReconstructable=0;
	// Loop over particles
  int nChargePart=0;
  int nCloseTrk=0;
	int nParticles = particleCollection->getNumberOfElements();
	for(int itParticle=0;itParticle<nParticles;itParticle++){
		 nCloseTrk=0;
		// Get the particle
		MCParticle* particle = dynamic_cast<MCParticle*>( particleCollection->getElementAt(itParticle) ) ;
    // get theta angle for plot: eff vs theta
    TLorentzVector tv_mcp;
    tv_mcp.SetPxPyPzE(particle->getMomentum()[0],particle->getMomentum()[1],particle->getMomentum()[2],particle->getEnergy());
    double theta_mcp=tv_mcp.Theta()*180./M_PI;
    double pt_mcp=tv_mcp.Pt();
    double vx=particle->getVertex()[0];
    double vy=particle->getVertex()[1];
    double vz=particle->getVertex()[2];
    double vr=sqrt(pow(vx,2)+pow(vy,2));

    if (m_morePlots){
      double charge_mcp=fabs(particle->getCharge());
      bool decay_mcp=particle->isDecayedInTracker();

      if (charge_mcp>0.5 && particle->getGeneratorStatus()==1){
        nChargePart++;

        m_mcCat.push_back(0);
        m_mcTheta.push_back(theta_mcp);
        m_mcPt.push_back(pt_mcp);
        m_mcIsDecayedInTracker.push_back(decay_mcp);  
        std::vector<TrackerHit*> trackHits_helper = particleHits[particle];
        int nhitsVXD = 0;
        for (size_t ihit=0; ihit<trackHits_helper.size(); ihit++){
          int subdetector = getSubdetector(trackHits_helper.at(ihit), m_encoder);
          if (subdetector == m_vertexBarrelID) nhitsVXD++;
          //m_mcNHits_helper.push_back(subdetector);
        }
        m_mcNHitsVXD.push_back(nhitsVXD);
        m_mcNHitsTot.push_back(trackHits_helper.size());
        //m_mcNHits.push_back(m_mcNHits_helper);
      }    
    }
		
    bool isReco=false;
		// Check if it was reconstructed
		if(particleTracks.count(particle)){
      
			// Assumption: if particle was reconstructed then it is reconstructable!
			m_particles["all"]++; m_particles["reconstructable"]++; nReconstructable++;
			m_reconstructedParticles["all"]++; nReconstructed++;
      // std::cout<< "--- nReconstructed = " << nReconstructed << std::endl;
      isReco=true;
			// Check if clones were produced (1 particle, more than 1 track)
			if(particleTracks[particle] > 1) m_reconstructedParticles["clones"]+=(particleTracks[particle]-1);

      // //fill theta distribution for reconstructable mc particle (denominator of the eff) and for the reco one (numerator of the eff)
      // h_theta_reconstructable -> Fill(theta_mcp);
      // h_theta_reconstructed -> Fill(theta_mcp);

      // h_pt_reconstructable -> Fill(pt_mcp);
      // h_pt_reconstructed -> Fill(pt_mcp);


      //fill also vector for tree variable - it would allow to make plots with different binning directly offline 

      m_vec_vx_reconstructable.push_back(vx);
      m_vec_vy_reconstructable.push_back(vy);
      m_vec_vz_reconstructable.push_back(vz);
      m_vec_vr_reconstructable.push_back(vr);
      m_vec_pt_reconstructable.push_back(pt_mcp);
      m_vec_theta_reconstructable.push_back(theta_mcp);
      m_vec_is_reconstructed.push_back(true);
  
      // m_vec_theta_reconstructed.push_back(theta_mcp);
      // m_vec_pt_reconstructed.push_back(pt_mcp);


      if (m_morePlots) {
        m_mcCat.pop_back();
        m_mcCat.push_back(2);
      }
			continue;
		}
		
		// Now decide on criteria for different particle types/classifications



		m_particles["all"]++; // all particles
		
		// No hits in the input collections
		if(particleHits.count(particle) == 0) continue; // No hits in the input collections


		// // Exclude particles which are not "real"
		// if(particle->vertexIsNotEndpointOfParent()) continue;

		
    if (!isReconstructable(particle,m_cuts)) continue;


		m_particles["reconstructable"]++; // reconstructable particles
		nReconstructable++;

    

    // //fill theta distribution for reconstructable mc particle (denominator of the eff)
    // h_theta_reconstructable -> Fill(theta_mcp);

    // h_pt_reconstructable -> Fill(pt_mcp);

    m_vec_vx_reconstructable.push_back(vx);
    m_vec_vy_reconstructable.push_back(vy);
    m_vec_vz_reconstructable.push_back(vz);
    m_vec_vr_reconstructable.push_back(vr);
    m_vec_pt_reconstructable.push_back(pt_mcp);
    m_vec_theta_reconstructable.push_back(theta_mcp);
    m_vec_is_reconstructed.push_back(false);


    if (m_morePlots) {
      m_mcCat.pop_back();
      m_mcCat.push_back(1);

      //for(int j=itParticle+1; j<nParticles; j++){
      for(int j=0; j<nParticles; j++){
        if (itParticle!=j){ 
          MCParticle* particle2 = dynamic_cast<MCParticle*>( particleCollection->getElementAt(j) );
          bool part2IsCharge = fabs(particle2->getCharge()) > 0.5;
          bool part2IsStable = particle2->getGeneratorStatus() == 1 ;
          if ( part2IsCharge && part2IsStable ) {            
            TLorentzVector tv_mcp2;
            tv_mcp2.SetPxPyPzE(particle2->getMomentum()[0],particle2->getMomentum()[1],particle2->getMomentum()[2],particle2->getEnergy());
            double DR = tv_mcp.DeltaR(tv_mcp2);
            if (DR<0.4) nCloseTrk++; 
          }
        }
      }

      if (!isReco) {
        m_mcNTracks.push_back(nChargePart);
        m_mcNTrkHits.push_back(m_mcNHitsTot.back());
        m_mcThetaTrk.push_back(m_mcTheta.back());
        m_mcPtTrk.push_back(m_mcPt.back());         
        m_mcNTracksCone.push_back(nCloseTrk);

        //skimVec->addElement( particle );

      }
    }

	}

  if (m_morePlots) {
    mctree->Fill();
	}

  trktree->Fill();
  puritytree->Fill();

  //evt->addCollection(  skimVec , m_notRecoMCColName ) ;


  struct sysinfo memInfo;

  sysinfo (&memInfo);

  unsigned long physMemUsed = memInfo.totalram - memInfo.freeram;
  //Multiply in next statement to avoid int overflow on right hand side...                                         
                                                                                                     
  physMemUsed *= memInfo.mem_unit;


  streamlog_out(MESSAGE) << " Processed event  "<<m_eventNumber
                         << " Physical memory in use: "<<physMemUsed
                         << std::endl ;

	// Increment the event number
	m_eventNumber++ ;
	std::cout<<"For this event reconstructed "<<100.*(double)nReconstructed/(double)nReconstructable<<" % ("<<nReconstructed<<"/"<<nReconstructable<<")"<<std::endl;
	
}

void ClicEfficiencyCalculator::check( LCEvent * evt ) {
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}


void ClicEfficiencyCalculator::end(){

	streamlog_out(MESSAGE) << " end()  " << name()
	<< " processed " << m_eventNumber << " events in " << m_runNumber << " runs "
	<< std::endl ;
	
	// Calculate efficiency results
	streamlog_out(MESSAGE)<<std::fixed<<std::setprecision(2)<<"Reconstructable particle efficiency: "<<100.*m_reconstructedParticles["all"]/m_particles["reconstructable"]<<" % ("<<std::setprecision(0)<<m_reconstructedParticles["all"]<<"/"<<m_particles["reconstructable"]<<")"<<std::endl;
	

  // // Plots

  // eff_vs_theta->cd();

  // g_eff_vs_theta->Divide( h_theta_reconstructed , h_theta_reconstructable , "v" ) ;
  // g_eff_vs_theta->SetMarkerColor(kRed) ;
  // g_eff_vs_theta->SetLineColor(kRed) ;
  // g_eff_vs_theta->GetYaxis()->SetTitle( "#epsilon_{trk}" );
  // g_eff_vs_theta->GetXaxis()->SetTitle( "#theta [rad]" );
  // g_eff_vs_theta->Draw("AP");
  // g_eff_vs_theta->Write();
  // eff_vs_theta->Write();

  // h_theta_reconstructed->Write();
  // h_theta_reconstructable->Write();


  // eff_vs_pt->cd();
  // gPad->SetLogx();

  // g_eff_vs_pt->Divide( h_pt_reconstructed , h_pt_reconstructable , "v" ) ;
  // g_eff_vs_pt->SetMarkerColor(kRed) ;
  // g_eff_vs_pt->SetLineColor(kRed) ;
  // g_eff_vs_pt->GetYaxis()->SetTitle( "#epsilon_{trk}" );
  // g_eff_vs_pt->GetXaxis()->SetTitle( "p_{T} [GeV]" );
  // g_eff_vs_pt->Draw("AP");
  // g_eff_vs_pt->Write();
  // eff_vs_pt->Write();

  // h_pt_reconstructed->Write();
  // h_pt_reconstructable->Write();

  trktree->Write();
  puritytree->Write();
  if (m_morePlots) mctree->Write();
  // _outRootFile->Write();
  // _outRootFile->Close();



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
	TrackerHitPlane* hit = dynamic_cast<TrackerHitPlane*>( collection->getElementAt(0) ) ;
	const int celId = hit->getCellID0() ;
	encoder.setValue(celId) ;
	int subdet = encoder[lcio::ILDCellID0::subdet];
	return subdet;
}


int ClicEfficiencyCalculator::getSubdetector(TrackerHit* hit, UTIL::BitField64 &encoder){
	const int celId = hit->getCellID0() ;
	encoder.setValue(celId) ;
	int subdet = encoder[lcio::ILDCellID0::subdet];
	return subdet;
}

int ClicEfficiencyCalculator::getLayer(TrackerHit* hit, UTIL::BitField64 &encoder){
	const int celId = hit->getCellID0() ;
	encoder.setValue(celId) ;
	int layer = encoder[lcio::ILDCellID0::layer];
	return layer;
}



bool ClicEfficiencyCalculator::isReconstructable(MCParticle*& particle, std::string cut){


  if (cut=="NHits") {

    // Only make tracks with 6 or more hits
		std::vector<TrackerHit*> trackHits = particleHits[particle];
		if(trackHits.size() >= 6) return true;

  } else if (cut=="NHitsVXD") {

    // Only make tracks with 4 or more hits in the vertex detector
		std::vector<TrackerHit*> trackHits = particleHits[particle];
    int nVXDHits = 0;
    UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ;
    for (size_t ihit=0; ihit<trackHits.size(); ihit++){
      int subdetector = getSubdetector(trackHits.at(ihit), encoder);
      if (subdetector == m_vertexBarrelID) nVXDHits++;
    }
    if (nVXDHits >= 4) return true;

  } else if (cut=="ILDLike") {


    // Only consider particles: charged, stable, pT>0.1GeV, cosTheta<0.89, nHits>=4, IP in 10 cm 
    //(a.t.m. cut in cosTheta is up to 0.89 instead of the usal 0.99 for cutting also particles in the vertex endcap region)
    //(a.t.m. no same cut on reco pT/theta but no bias because in that case also the mc particle is kept)

    //add condition for decay outside vertex/tracker? 
    bool isCharge = false;
    bool isStable = false;
    bool passPt = false;
    bool passTheta = false;
    bool passNHits = false;
    bool passIP = false;
    bool passEndPoint = false;
    // bool isOverlay = particle->isOverlay();
    // if (isOverlay) std::cout<<"----- ehi there are overlay particles in this mc collection"<<std::endl;

    double charge = fabs(particle->getCharge());
    if (charge>0.5) isCharge = true;
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
		if(trackHits.size() >= 4) passNHits = true;
    //else std::cout<<"----- mc has not 4 hits associated"<<std::endl;
		//if(trackHits.size() >= 6) passNHits = true;

    // int nVXDHits = 0;
    // UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ;
    // for (size_t ihit=0; ihit<trackHits.size(); ihit++){
    //   int subdetector = getSubdetector(trackHits.at(ihit), encoder);
    //   if (subdetector == m_vertexBarrelID) nVXDHits++;
    // }
    // if (nVXDHits >= 4) return true;

    double dist = sqrt( pow(particle->getVertex()[0],2) + pow(particle->getVertex()[1],2) );
    if (dist<100.) passIP = true;
    //else std::cout<<"----- mc has dist from IP > 100mm"<<std::endl;
    //if (dist<30.) passIP = true;
    
    double e = sqrt( pow(particle->getEndpoint()[0],2) + pow(particle->getEndpoint()[1],2) );
    if (e==0. || e>40.) passEndPoint=true;
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
    UTIL::BitField64 encoder( lcio::ILDCellID0::encoder_string ) ;
    std::vector<int > vec_hit_subdet;
    std::vector<int > vec_hit_layer;
    for (size_t ihit=0; ihit<trackHits.size(); ihit++){
      int subdetector = getSubdetector(trackHits.at(ihit), encoder);
      int layer = getLayer(trackHits.at(ihit), encoder);
      if (subdetector == m_vertexBarrelID) nVXDHits++;
      vec_hit_subdet.push_back(subdetector);
      vec_hit_layer.push_back(layer);
    }
    //if (nVXDHits >= 6) passVertexHits = true;

    for (size_t j = 0; j < vec_hit_subdet.size()-1; j++){
      for (size_t k = j+1; k < vec_hit_subdet.size(); k++){
        if ( vec_hit_subdet.at(j) == vec_hit_subdet.at(k) ) {
          if ( vec_hit_layer.at(j) == vec_hit_layer.at(k) ) {
            isNotLoop = false;
            std::cout << "=====> loop track " << std::endl;
          }
        }
      }
    }

    vec_hit_subdet.clear();
    vec_hit_layer.clear();

    double dist = sqrt( pow(particle->getVertex()[0],2) + pow(particle->getVertex()[1],2) );
    std::cout<<"----- for mc_pdg = "<< pdg << "   dist from IP = " << dist <<std::endl;
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



void ClicEfficiencyCalculator::clearTreeVar(){

  m_mcCat.clear();
  m_mcTheta.clear();
  m_mcPt.clear();
  m_mcTheta.clear();
  m_mcIsDecayedInTracker.clear();
  //m_mcNHits_helper.clear();
  //m_mcNHits.clear();
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

  // m_vec_theta_reconstructed.clear(); 
  // m_vec_pt_reconstructed.clear(); 

  m_vec_nhits_vtx.clear();
  m_vec_nhits_trk.clear();
  m_vec_nhits.clear();
  m_vec_purity.clear();
  m_vec_pdg.clear();
  m_vec_theta.clear();
  m_vec_phi.clear();
  m_vec_p.clear();


}  
  
  
  
  
  
