# v02-04-01

* 2020-09-21 Andre Sailer ([PR#132](https://github.com/iLCSoft/ClicPerformance/pull/132))
  - CLIC Simulation: change the path to the particle.tbl file
  - FCC Simulation: change the path to the particle.tbl file
  - CLIC Reconstruction: change path of default detector and input file, fixes #131 
  - FCC Reconstruction: change path of default  input file, fixes #131

* 2020-06-03 Erica Brondolin ([PR#130](https://github.com/iLCSoft/ClicPerformance/pull/130))
  - Adapt `fcc-` and `clicReconstruction.xml` according to the [PR#45 in MarlinTrkProcessors](https://github.com/iLCSoft/MarlinTrkProcessors/pull/45)
  - `mergeSplitTracks` is set to true in the case of `clicReconstruction.xml`

* 2020-03-31 Emilia Leogrande ([PR#128](https://github.com/iLCSoft/ClicPerformance/pull/128))
  - New cuts for new definitions of reconstructability added to ClicEfficiencyCalculator
  - Cut "AnyGenStatus" is optimized for charginos (no requirement on number of hits and on stability)
  - Cut "AnyGenStatusLowPt" is optimized for soft pions from chargino decays (like "AnyGenStatus" but without requirement on pt)

* 2020-02-07 Andre Sailer ([PR#127](https://github.com/iLCSoft/ClicPerformance/pull/127))
  - ClicReconstruction: Reduce ConformalTracking.TooManyTracks to 100k from 500k
  - FccReconstruction: Reduce ConformalTracking.TooManyTracks to 100k from 500k

* 2020-02-05 Erica Brondolin ([PR#126](https://github.com/iLCSoft/ClicPerformance/pull/126))
  - Add track purity in perfTree for tracking validation (already present in purTree)

* 2020-01-24 Matthias Artur Weber ([PR#125](https://github.com/iLCSoft/ClicPerformance/pull/125))
  - 1. JetAnalyzer for Jet and missing (transverse) energy studies
  - 2. TrueMCintoRecoForJets MCParticle and ReconstructedParticle for jets selector, giving functionality to select different input particles
      a) options for MC particles: stable particles, including all neutrinos, including all neutrinos from hadronic decays
     b) remove reconstructed particles after angular matching to leptons from vector boson decays

* 2019-10-04 Erica Brondolin ([PR#124](https://github.com/iLCSoft/ClicPerformance/pull/124))
  - CLIC Reco: ConformalTracking: Set HighPTCut to 0
      - According to a study on the conformal tracking fit for the ExtendTracks step, the quadratic term in the fit is not needed for the low-pT tracks. Results for ttbar w/o and w/ overlay shows unchanged level of tracking efficiency, slightly increase of fakerate, and reduction CPU time of 13%. Full set of results [here](https://indico.cern.ch/event/848916/contributions/3567053/attachments/1909256/3154418/EricaBrondolin_20190917_LowPtStudies.pdf).

* 2019-10-02 Andre Sailer ([PR#123](https://github.com/iLCSoft/ClicPerformance/pull/123))
  - FCC/CLICReco:  fix 'spelling' of TooManyTracks parameter value

* 2019-09-13 Andre Sailer ([PR#119](https://github.com/iLCSoft/ClicPerformance/pull/119))
  - clic_steer, fcc_steer: use g4units instead of dropped SystemOfUnits

* 2019-08-26 Andre Sailer ([PR#118](https://github.com/iLCSoft/ClicPerformance/pull/118))
  - DDSim Steer: fix spelling of ClassicalRK4 for magnetic field stepper, this is also the geant4 default

* 2019-08-26 Erica Brondolin ([PR#115](https://github.com/iLCSoft/ClicPerformance/pull/115))
  - Update xml file with z cut on ConformalTracking (https://github.com/iLCSoft/ConformalTracking/pull/53)

* 2019-08-19 Emilia Leogrande ([PR#117](https://github.com/iLCSoft/ClicPerformance/pull/117))
  - ClicEfficiencyCalculator: fixed m_vec_is_reconstructed branch for not reconstructed tracks in trackTree, fixes #116

* 2019-07-08 Andre Sailer ([PR#114](https://github.com/iLCSoft/ClicPerformance/pull/114))
  - FCC/ClicReco: add optional parameters for Output_REC, Output_DST so that they can be overwritten via the command line

# v02-03

* 2019-02-20 Oleksandr Viazlo ([PR#106](https://github.com/ilcsoft/ClicPerformance/pull/106))
  - add PandoraSettingsPhotonTraining.xml needed for calorimeter calibration service

* 2019-02-07 Emilia Leogrande ([PR#105](https://github.com/ilcsoft/ClicPerformance/pull/105))
  - fccReconstruction:: JetClusteringAndRefiner: B field changed from 4 T to 2 T

* 2018-11-19 Oleksandr Viazlo ([PR#104](https://github.com/ilcsoft/ClicPerformance/pull/104))
  - new calorimeter calibration for CLD detector with 400ns integration time window; 
  - refactoring of fccReconstruction.xml to enable possibility to switch between settings corresponding to 10ns or 400ns integration time window

# v02-02

* 2018-11-14 Marko Petric ([PR#103](https://github.com/ilcsoft/ClicPerformance/pull/103))
  - Changes to colour and layout of the visualization model

* 2018-11-06 Marko Petric ([PR#102](https://github.com/ilcsoft/ClicPerformance/pull/102))
  - Enable energy coloring in visualization in CLIC visualization template and use new colors for the event display

* 2018-11-01 Emilia Leogrande ([PR#101](https://github.com/ilcsoft/ClicPerformance/pull/101))
  - CLIC: beam sizes from CDR
  - FCCee: beam sizes from FCC Week 2018

* 2018-11-01 Andre Sailer ([PR#100](https://github.com/ilcsoft/ClicPerformance/pull/100))
  - ClicReconstruction::Vertexing add JetClustering.PrimaryVertexCollectionName parameter
  - CI: simulate and reconstruct some real events
  - CI: move truth tracking reconstruction to TRUTH configuration
  - CI: remove overlay test as we don't have an overlay file at the moment

* 2018-10-31 Emilia Leogrande ([PR#99](https://github.com/ilcsoft/ClicPerformance/pull/99))
  * ClicReconstruction: Setup for running the vertexing and jet clustering
    - Prior to vertexing, FastJet is run in case of overlay. Without overlay, dummy MergeCollections is used to shallow copy the PandoraPFOs collection and have the same output collection name as in case of overlay
    - VertexFinder and JetClusteringAndRefiner are run
    - if VertexUnconstrainedOn, also unconstrained vertex finder is run (for vertex resolutions)

* 2018-10-31 Andre Sailer ([PR#98](https://github.com/ilcsoft/ClicPerformance/pull/98))
  -  ClicRecoConfig: add Not + option with negative condition,  'Config.OverlayNotFalse' which would be true if the choice is not False (3TeV, etc

* 2018-10-31 Andre Sailer ([PR#95](https://github.com/ilcsoft/ClicPerformance/pull/95))
  - CI: added -Werror flag to CI compilation

* 2018-10-26 Emilia Leogrande ([PR#97](https://github.com/ilcsoft/ClicPerformance/pull/97))
  - Two vertexFinder processors defined in a group, one for flavour tagging (beamSpotConstraint = 1), one for vertex resolution (beamSpotConstraint = 0)
  - Option 'VertexUnconstrained' in CLICRecoConfig, to enable the vertexFinderUnconstrained and the JetClustering. Default is OFF. The vertexFinder constrained is always enabled
  - Same updates in clicReconstruction.xml and fccReconstruction.xml
  - Changed LumiCalReco verbosity (from DEBUG0 to WARNING)
  - Removed unused variable uniqueHits in ClicEfficiencyCalculator.cc
  - Added 'deg' units in clic_steer.py and fcc_steer.py
  - Added folder vtxprob/ in Flavour_tagging. Needed to run the makeNtuple

* 2018-10-25 Oleksandr Viazlo ([PR#96](https://github.com/ilcsoft/ClicPerformance/pull/96))
  - FCCee reconstruction: update track selection cuts (d0, z0, hit radius) to 200mm

* 2018-10-23 Andre Sailer ([PR#93](https://github.com/ilcsoft/ClicPerformance/pull/93))
  - ClicReco: PandoraPFA: update z0, d0, MaxBarrelTrackerInnerRDistance track cuts to 200 mm

* 2018-10-17 Emilia Leogrande ([PR#92](https://github.com/ilcsoft/ClicPerformance/pull/92))
  - Added VertexFinder to default clicReconstruction steering file. Relates to #68

* 2018-10-17 Emilia Leogrande ([PR#91](https://github.com/ilcsoft/ClicPerformance/pull/91))
  - Updated instructions on vertex reconstruction for flavour tagging

# v02-01

* 2018-10-05 Erica Brondolin ([PR#90](https://github.com/ilcsoft/ClicPerformance/pull/90))
  - In the tree the nHitsMC must be the number of hits that are reconstructed and that belong to the same particle.

* 2018-08-23 Oleksandr Viazlo ([PR#89](https://github.com/ilcsoft/ClicPerformance/pull/89))
  - Update FCCee_o1_v03_CED
    - add magnetic field
    - change color of Yoke

* 2018-08-22 Marko Petric ([PR#88](https://github.com/ilcsoft/ClicPerformance/pull/88))
  - Add `FCCee_o1_v03_CED` model into Visualisation folder. Model only suitable for event display not for simulation!
    - example usage: `ced2go -t ced2go-template-DD4.xml -d FCCee_o1_v03_CED.xml data/REC/00010521/000/Z_uds_rec_10521_105.slcio`
  - Final colour scheme still has to be decided by FCCee people

* 2018-08-09 Emilia Leogrande ([PR#87](https://github.com/ilcsoft/ClicPerformance/pull/87))
  - ConformalTracking: Set pT threshold to run extendHighPt as introduced in iLCSoft/ConformalTracking#42
    - Added to both clic and fccee steering files

* 2018-07-06 Erica Brondolin ([PR#86](https://github.com/ilcsoft/ClicPerformance/pull/86))
  - Clic Reconstruction: Update the Low Energy timing cuts: the variables "HcalBarrelLoose(/Tight)TimingCut" are now set consistently with the loose (tight) timing cuts imposed on neutral hadrons.

* 2018-06-29 Emilia Leogrande ([PR#85](https://github.com/ilcsoft/ClicPerformance/pull/85))
  - clicReconstruction.xml and fccReconstruction.xml: setup the configuration for ConformalTrackingV2
  - parameters match the configuration of old ConformalTracking
  - NB: parameters for clic and fcc are supposed to be different

* 2018-06-29 Andre Sailer ([PR#84](https://github.com/ilcsoft/ClicPerformance/pull/84))
  - Splitting fcc and clic configuration files into separate folders: clicConfig and fcceeConfig

* 2018-06-27 Oleksandr Viazlo ([PR#83](https://github.com/ilcsoft/ClicPerformance/pull/83))
  - enable Software Compensation for CLD
  - CLD calorimeter calibration constants for new default settings
  - correct CoilCorrectionMinInnerRadius value for CLD

* 2018-06-25 Emilia Leogrande ([PR#81](https://github.com/ilcsoft/ClicPerformance/pull/81))
  - Updates in the fccReconstruction.xml:
  - MCPhysicsParticles (instead of MCParticle) given as input to the RecoMCTruthLinker to make the MCParticlesSkimmed
  - SiTracks_Refitted added to the output collections in the DST
  - fixed setup for Overlay group

* 2018-05-18 Emilia Leogrande ([PR#80](https://github.com/ilcsoft/ClicPerformance/pull/80))
  - clicReconstruction: in ConformalTracking, added flag to enable tight cuts in the vertex combined (barrel+endcap) hits reconstruction
  - done similarly to fccReconstruction
  - to be upgraded soon with a more elegant solution

* 2018-05-16 Emilia Leogrande ([PR#79](https://github.com/ilcsoft/ClicPerformance/pull/79))
  - fccReconstruction: in ConformalTracking settings, added processor parameter to enable/disable the tight cuts for the combined vertex barrel + endcap reconstruction 
  - default = true (enabled)
  - for FCCee, it needs to be disabled because of 10degrees tracks (2 hits in the vertex barrel + 6 in the endcap). Looser cuts allow to pick the vertex (b+e) hits in one track, otherwise tight cuts allow for a track made with a subset of hits, while the leftover hits are used for a second track
  - the external processor parameter is a temporary solution

* 2018-04-26 Emilia Leogrande ([PR#77](https://github.com/ilcsoft/ClicPerformance/pull/77))
  - Tracking/include(src)/ClicTrackingValidate.h(cc)
  - Processor to validate the tracking using MCParticle information at particle and hit level
  - Checks the misassociated hits (from wrong MCParticle) and the missed hits

* 2018-04-26 Andre Sailer ([PR#67](https://github.com/ilcsoft/ClicPerformance/pull/67))
  - Reconstruction: add all parameters for DDMarlinPandora, fixes #66 
  - Reconstruction: DDMarlinPandora increase D0 and Z0 cut offs
  - Reconstruction: Add MCPhysicsParticles to DST output file

* 2018-04-24 Emilia Leogrande ([PR#76](https://github.com/ilcsoft/ClicPerformance/pull/76))
  - clicReconstruction.xml: Config.TrackingConformal runs the ConformalTracking + ClonesAndSplitTracksFinder, BUT with the option of mergeSplitTracks = false
  - mergeSplitTracks = true would enable also the merging of close tracks, which is still work in progress
  - same implementation updated also for fccReconstruction.xml

* 2018-04-19 Emilia Leogrande ([PR#75](https://github.com/ilcsoft/ClicPerformance/pull/75))
  - clicReconstruction.xml : edited AIDAprocessor root output file name
  - fccReconstruction.xml: edited AIDAprocessor root output file name, set new detector model and conformal tracking as default

* 2018-04-19 Andre Sailer ([PR#74](https://github.com/ilcsoft/ClicPerformance/pull/74))
  - CLICReco: ConformalTracking: add the new parameters for collection names in VXDBarrel/Endcap and MainTracker
  - FCCReco: ConformalTracking: add the new parameters for collection names in VXDBarrel/Endcap and MainTracker

* 2018-04-17 Andre Sailer ([PR#73](https://github.com/ilcsoft/ClicPerformance/pull/73))
  - ClicReconstruction: drop ConformalTracking+Extrapolator option for Tracking
  - FCCReconstruction: drop ConformalTracking+Extrapolator option for Tracking
  - ClicRecoConfig: drop alternative and broken Config names
  - FCCReconstruction: use correct TrackingChoices

* 2018-03-23 Andre Sailer ([PR#72](https://github.com/ilcsoft/ClicPerformance/pull/72))
  -  ClicReconstruction::Overlay: correct the eventsPerBX numbers for L*=6m overlay instances

* 2018-03-09 Nacho Garcia ([PR#71](https://github.com/ilcsoft/ClicPerformance/pull/71))
  - Code to run the flavour-tagging training and produce performance plots.

* 2018-03-05 Andre Sailer ([PR#70](https://github.com/ilcsoft/ClicPerformance/pull/70))
  - Reconstruction: add overlay for 350GeV, 380GeV, and 3TeV L*=6m layout and L*=4.3m, keep CDR values in separate processors

* 2018-03-02 Emilia Leogrande ([PR#69](https://github.com/ilcsoft/ClicPerformance/pull/69))
  - examples/clicReconstruction.xml: added configuration for ClonesAndSplitTracksFinder
  - the processor is included in the Config.TrackingConformal condition, to run automatically after the conformal tracking
  - the output track collection from ClonesAndSplitTracksFinder is given as input to RefitFinal
  - analogous for fccReconstruction.xml

* 2018-02-07 Emilia Leogrande ([PR#61](https://github.com/ilcsoft/ClicPerformance/pull/61))
  - Tests/CMakeLists.txt: Added 3 tests for FCCee_o1_v02:
  - sim
  - reco_truth
  - reco_conformal

* 2018-02-06 Emilia Leogrande ([PR#65](https://github.com/ilcsoft/ClicPerformance/pull/65))
  - **fcc_steer.py**: Lorentz boost set to 0.015 rad [NB: to be used on pairs files only if those have been produced without crossing angle. If this is not the case, set the Lorentz boost to 0.0]
  - **fccReconstruction.xml**: adjusted Overlay parameters to FCCee (number of BX, Delta_t, integration time windows per subdetector, number of events to overlay = 1 -- not Poissonian -- given that the gen files from GuineaPig contain 1 evt per file) and removed BeamCalReco (not available for FCCee yet)
  - **CLICRecoConfig.h**: added two options for FCCee background overlay (91GeV, 365GeV)

* 2018-02-06 Andre Sailer ([PR#64](https://github.com/ilcsoft/ClicPerformance/pull/64))
  - Simulation: add Lorentz boost for the crossing angle to clic_steer.py

* 2018-02-02 Andre Sailer ([PR#63](https://github.com/ilcsoft/ClicPerformance/pull/63))
  - Reconstruction: add subClusterEnergyID to LumiCalReco and BeamCalReco

* 2018-01-29 Andre Sailer ([PR#62](https://github.com/ilcsoft/ClicPerformance/pull/62))
  - Use BeamCalReco for LumiCal reconstruction
  - Merge Cluster collections and feed all to RecoMCTruthLiner
  - Merge ReconstructedParticle collections and feed all o RecoMCTruthLinker
  - Only gives physics particles to RecoMCTruthLinker
  - Switch Tests to CLIC_o3_v14

* 2018-01-18 Andre Sailer ([PR#60](https://github.com/ilcsoft/ClicPerformance/pull/60))
  - Reconstruction: DST Files: add low energy Selected PFO collection
  - Fix typo in LE_TightSelectedPandoraPFOs

* 2017-12-21 Emilia Leogrande ([PR#59](https://github.com/ilcsoft/ClicPerformance/pull/59))
  - Updated steering files clicReconstruction.xml and fccReconstruction.xml
  - Single point resolutions in the tracker (inner and outer, barrel and endcap): 7um x 90um
  - all layers except 1st inner tracker endcap: still 5um x 5um

* 2017-12-12 Andre Sailer ([PR#57](https://github.com/ilcsoft/ClicPerformance/pull/57))
  - Refactoring Config Processor: allow extending options via parameters, easier to add new options
  - Add 380GeV BeamCalBackground file based on 380GeV L*=6m (http://clic-beam-beam.web.cern.ch/clic-beam-beam/380gev_l6_bx8mm.html)
  - Waiting for ilcsoft/Marlin#28 to be deployed to optimize selection of processors, which makes the Config processor obsolete and command line arguments would would be `--constant....` instead of `--Config...`

* 2017-12-07 Matthias ([PR#56](https://github.com/ilcsoft/ClicPerformance/pull/56))
  - new calibration constants for detector model CLIC_o3_v13
  - new photonlikelihood file, now derived in 12 instead of 9 energy bins
  - CLIC specific compensation weights, extend range and tuning in energy and energy density. Weights are derived from a sum of single particle neutron and K0L events (ratio 1:1), energy points from 2 GeV up to 1000 GeV

* 2017-12-07 Andre Sailer ([PR#55](https://github.com/ilcsoft/ClicPerformance/pull/55))
  - CLIC Reconstruction:
    - Add SiTracks_Refitted to DST Output, add LumiCal and BeamCal Cluster and RecoParticle to DST Output
    - Note that BeamCalReconstruction output collections have changed names to be consistent with LumiCalOutput
    - Add new options to ConformalTracking (iLCSoft/ConformalTracking#28, iLCSoft/ConformalTracking#29)
    - add call for running vtune on reconstruction (only works in CERN intranet)
       `ctest  -C PROF -R t_test_CLIC_o3_v13_reco_conformal_zuds_overlay_profile -V`

* 2017-12-06 Oleksandr Viazlo ([PR#54](https://github.com/ilcsoft/ClicPerformance/pull/54))
  - Pandora settings for FCC-ee detector (no software compensation)

* 2017-12-05 Emilia Leogrande ([PR#50](https://github.com/ilcsoft/ClicPerformance/pull/50))
  - `clicReconstruction.xml`: removed unused processor (`DDCellsAutomatonMV`)
  - `fccReconstruction.xml`: new steering file for reconstruction with the FCCee detector (conformal tracking parameters are different than for CLIC)

* 2017-12-04 Marko Petric ([PR#53](https://github.com/ilcsoft/ClicPerformance/pull/53))
  - Update field stepping to Geant4 defaults

* 2017-11-27 Emilia Leogrande ([PR#51](https://github.com/ilcsoft/ClicPerformance/pull/51))
  - ClicEfficiencyCalculator: `perftree` copied from `TrackChecker` 
  - `trueVertexR` info added
  - innermost hit Radius added

* 2017-11-08 Emilia Leogrande ([PR#49](https://github.com/ilcsoft/ClicPerformance/pull/49))
  - ClicEfficiencyCalculator: fixed storing of track purity: now done with a std::list<std::pair<MCParticle*,double>>, which stores the MCParticle associated to the track and the track purity

* 2017-11-08 Emilia Leogrande ([PR#48](https://github.com/ilcsoft/ClicPerformance/pull/48))
  - ClicEfficiencyCalculator: added 'isSignal' flag for MCPhysicsParticles
  - MCPhysicsParticles collection is the output of the Overlay Processor
  - in case of no overlay, the MCParticle collection is used instead
  - clicReconstruction.xml: MCPhysicsParticles collection as input

* 2017-10-20 Emilia Leogrande ([PR#47](https://github.com/ilcsoft/ClicPerformance/pull/47))
  - ClicEfficiencyCalculator: enabled simpleOutput tree (very important for efficiency plots)

* 2017-10-18 Emilia Leogrande ([PR#45](https://github.com/ilcsoft/ClicPerformance/pull/45))
  - ClicEfficiencyCalculator: updated ILDLike cuts to accept unique hits only, i.e. hits from different subdetector layers, as in the default NHits cuts

* 2017-10-18 Andre Sailer ([PR#44](https://github.com/ilcsoft/ClicPerformance/pull/44))
  - Updated Tests to CLIC_o3_v13
  - Added tests running Z_uds events if EOS is available
  - Added more tests running valgrind (memcheck and massif), only if enabled

* 2017-10-13 Emilia Leogrande ([PR#43](https://github.com/ilcsoft/ClicPerformance/pull/43))
  Updated parameters for Conformal Tracking in clicReconstruction.xml
  - ThetaRange: from 0.03 to 0.05
  - MinClustersOnTrack: from 6 to 4
  - MaxChi2: from 50 to 100

* 2017-10-06 Andre Sailer ([PR#42](https://github.com/ilcsoft/ClicPerformance/pull/42))
  - Drop unused and no longer existing header includes AidaSoft/DD4hep#241

* 2017-10-02 Emilia Leogrande ([PR#41](https://github.com/ilcsoft/ClicPerformance/pull/41))
  - ClicEfficiencyCalculator: added min distance between MCParticles in the simplifiedEfficiencyTree output
  - The distance is computed as the angular separation DeltaR = sqrt( DeltaEta^2 + DeltaPhi^2 )

* 2017-09-15 Emilia Leogrande ([PR#40](https://github.com/ilcsoft/ClicPerformance/pull/40))
  - Removed purity selection on reconstructed tracks to compute efficiency. This makes tracking efficiency independent on purity.
  - Added purity info in simplifiedEfficiencyTree: can be used to cut offline on the purity.

* 2017-09-14 Emilia Leogrande ([PR#39](https://github.com/ilcsoft/ClicPerformance/pull/39))
  - TrackChecker: added MCParticle PDG info in checktree

* 2017-08-31 Andre Sailer ([PR#38](https://github.com/ilcsoft/ClicPerformance/pull/38))
  - CLICEfficiency/HitResiduals/TrackChecker: use streamlog instead of cout

* 2017-08-25 Andre Sailer ([PR#37](https://github.com/ilcsoft/ClicPerformance/pull/37))
  - HitResiduals: protect against missing trackstate
  - HitResiduals: use AtIP instead of AtFirstHit
  - Reconstruction: fix typo in parameter value for RefitFinal (inconsequential as using default value anyway)
  - Reconstruction: use SiTracks_Refitted for HitResiduals and CLICEfficiencyCalculator

* 2017-08-24 Emilia Leogrande ([PR#35](https://github.com/ilcsoft/ClicPerformance/pull/35))
  * CLICRecoConfig: Updated Config for Tracking: Truth, ConformalPlusExtrapolator, Conformal

* 2017-08-23 Matthias ([PR#36](https://github.com/ilcsoft/ClicPerformance/pull/36))
  - Update in PandoraSettingsDefault, enable SoftwareCompensation, including new cleaning of clusters
  - change calibration constants in clicReconstruction to reflect changes in new Pandora default settings as well as switch to model CLIC_o3_v13, effectively removing the MaxHCalHitHadronicEnergy cut, used to scale hot hadrons with the previous PandoraSettings
  - newly adapted PhotonLikelihoodFile

* 2017-08-22 Andre Sailer ([PR#34](https://github.com/ilcsoft/ClicPerformance/pull/34))
  - Reco: add refit processor after the track reconstruction

* 2017-08-21 Andre Sailer ([PR#33](https://github.com/ilcsoft/ClicPerformance/pull/33))
  - Overlay processors now in a group to avoid duplicating parameters
  - Add overlay parameters for 350, 420, 500, and 1.4TeV

* 2017-07-24 Andre Sailer ([PR#32](https://github.com/ilcsoft/ClicPerformance/pull/32))
  - Add processor printing out current event number
  - overlay: add parameters for 380GeV
  - overlay: change Number of BX to overlay to 30 fitting the time interval

# v02-00-00

* 2017-04-21 Andre Sailer ([PR#10](https://github.com/ilcsoft/CLICPerformance/pull/10))
  - Reconstruction::MuonReconstruction change number of required layers for muon identification to reflect smaller number of layers in the muon system. Increases Muon ID efficiency

* 2017-04-24 Andre Sailer ([PR#11](https://github.com/ilcsoft/CLICPerformance/pull/11))
  - Reconstruction::MuonReco: adapt more parameters to new geometry, as recommended by @bonoxu

* 2017-05-24 Daniel Hynds ([PR#12](https://github.com/ilcsoft/CLICPerformance/pull/12))
  - ClicEfficiencyCalculator: Added number of hits on track to the simplified tree output
  - Reconstruction: Updated parameters for the ConformalTracking

* 2017-05-24 Matthias ([PR#13](https://github.com/ilcsoft/CLICPerformance/pull/13))
  - Reconstruction: new calibration constants for calorimeters and MUON system, using model CLIC_o3_v10. Photon likelihood file modified to adapt to those numbers.

* 2017-05-26 Andre Sailer ([PR#14](https://github.com/ilcsoft/CLICPerformance/pull/14))
  - Add PFOSelector processors, copied from CDR reconstruction

* 2017-05-30 Andre Sailer ([PR#17](https://github.com/ilcsoft/CLICPerformance/pull/17))
  - Add CLICRecoConfig processor to chose which configuration to use for reconstruction: truth conformal tracking, if and which overlay

* 2017-05-30 Andre Sailer ([PR#19](https://github.com/ilcsoft/CLICPerformance/pull/19))
  - ClicEfficiencyCalculator: Fix #18 , comment code storing information when particle not reconstructed, however some of that information is not aggregated beforehand so there was a crash

* 2017-06-19 Emilia Leogrande ([PR#20](https://github.com/ilcsoft/CLICPerformance/pull/20))
  * HitResiduals: replaced 'subdet' with 'system' in the call to the encoding string (otherwise the processor won't work) and added the 'side' information needed to propagateToLayer when hits are in the endcaps.

* 2017-06-20 Andre Sailer ([PR#21](https://github.com/ilcsoft/CLICPerformance/pull/21))
  - Adapt to the rest of namespace changes in DD4hep

* 2017-06-22 Andre Sailer ([PR#23](https://github.com/ilcsoft/CLICPerformance/pull/23))
  - tests running valgrind:
      By default uses valgrind from /cvmfs/sft
        to run:
          `ctest -C valgrind -R _testName_`
       See test names via `ctest -C valgrind -N`
  - Added option to run tests for CLIC_o2_v04 model for faster debugging 
     E.g.: `ctest -C o2_v04 -R o2_v04 -E valgrind`
    needs special steering file without GlobalTrackerReadoutID specified
  
  - Add suppressions file for valgrind to use in addition to root suppressions file for uninitialised member in a TObject

* 2017-06-28 Andre Sailer ([PR#24](https://github.com/ilcsoft/CLICPerformance/pull/24))
  * TrackChecker: 
    * use c++ function instead of string for TF1, reduced memory allocation by 50MB
    * cleanup included header files
    * move creation of histograms and tree to init, register with AIDAProcessor to create own folder
  
  * HitResiduals:
    * cleanup; use AIDArocessor to store output
    * clean marlin_trk object
    * use numeric indices to acces subdet, side and layer
  
  * ClicEfficiencyCalculator:
    * comment unused variables for selection cuts
    * static_cast --> dynamic cast
    * move clearing of tree variables to end of event
    * cleanup header files, drop registration with eventseeder
    * wrap relations in shared_ptr to clean them automatically

* 2017-06-28 Andre Sailer ([PR#25](https://github.com/ilcsoft/CLICPerformance/pull/25))
  - ReconstructionSteering: remove duplicate setting of DebugPlots for ConformalTracking

* 2017-06-30 Andre Sailer ([PR#26](https://github.com/ilcsoft/CLICPerformance/pull/26))
  - Ensure that histograms and trees are written only to their folder in the RootFile. AIDAProcessor already writes the outputs so we don't have to as well
  - TrackChecker: beautify fit output by naming parameters. Fix issue with unnormalised gauss function

* 2017-07-05 Andre Sailer ([PR#27](https://github.com/ilcsoft/CLICPerformance/pull/27))
  - Add suppressions for root things, lcio buffers etc.

* 2017-07-10 Andre Sailer ([PR#28](https://github.com/ilcsoft/CLICPerformance/pull/28))
  - Add suppressions for 4-8 byte definite leaks in root compiled with llvm39
  - Add suppressions for more still reachable leaks in root (based on gcc)

* 2017-07-11 Andre Sailer ([PR#29](https://github.com/ilcsoft/CLICPerformance/pull/29))
  - Reconstruction: Add REC/DST OutputProcessors
  - Reconstruction: Change default name of the output file to Output_REC.slcio and Output_DST.slcio instead of sitracks.slcio
  - Reconstruction: Change name of `SelectedPandoraPFO` collections. Instead of `SelectedPandoraPFANewPFO` collections

* 2017-07-11 Andre Sailer ([PR#30](https://github.com/ilcsoft/CLICPerformance/pull/30))
  - Add LumiCalReco processor
  - Add BeamCalReco processor

* 2017-07-11 Daniel Hynds ([PR#31](https://github.com/ilcsoft/CLICPerformance/pull/31))
  - updated chi2 for conformal tracking to account for new track fit
  - removed tracker hit collections from steering file. If included, the conformal tracking will extend tracks through them and the extrapolator should not be used. Still under testing

# v01-00
*  General: At present the package contains only processors for the Tracking systems. Further dedicated folders for Calorimetry performance etc. should be added in future. The configuration folder contains steering files for the reconstruction in Marlin, and a test script (testClicBuild.sh) which allows a fast check of new relseases/geometries etc.
*  ClicEfficiencyCalculator: Processor to calculate the tracking efficiency. Only the input hit collections are considered, allowing the efficiency to be restricted to the given subdetectors.
*  TrackChecker: Processor to evaluate the pull and residual distributions of the track fitting. 

