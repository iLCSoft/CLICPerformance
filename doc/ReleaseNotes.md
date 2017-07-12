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

