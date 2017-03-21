# v01-00
*  General: At present the package contains only processors for the Tracking systems. Further dedicated folders for Calorimetry performance etc. should be added in future. The configuration folder contains steering files for the reconstruction in Marlin, and a test script (testClicBuild.sh) which allows a fast check of new relseases/geometries etc.
*  ClicEfficiencyCalculator: Processor to calculate the tracking efficiency. Only the input hit collections are considered, allowing the efficiency to be restricted to the given subdetectors.
*  TrackChecker: Processor to evaluate the pull and residual distributions of the track fitting. 

