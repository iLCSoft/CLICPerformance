#!/bin/bash

# check if the detector model has been specified
ARGS=1
if [ $# -ne "$ARGS" ]
then 
  echo -e "\nPlease give the detector model to run over, eg. CLIC_o2_v01\n"
  exit
fi

# set up the environment
export MARLIN_DLL=${MARLIN_DLL}../lib/libClicPerformance.so:
export PYTHONPATH=${LCIO}/src/python:${ROOTSYS}/lib
source $DD4HEP/bin/thisdd4hep.sh
 
# find the detector model specified
DETECTOR="$(find $ILCSOFT/lcgeo -name ${1}.xml)" 
echo -e "\nPicking up detector from: $DETECTOR"

# generate 1000 muons in random directions
echo "Generating particles"
python lcio_particle_gun.py &> particle.log

# run geant4 with the detector model specified
echo "Running GEANT4"
python lcio_geant4.py $DETECTOR &> geant4.log

# run the CLIC reconstruction chain includin all tracking detectors
Marlin clicReconstruction.xml  --global.LCIOInputFiles=simulatedMuons.slcio

