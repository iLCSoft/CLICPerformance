#!/bin/bash

# check if the detector model has been specified
ARGS=1
if [ $# -ne "$ARGS" ]
then
  echo "Please give the detector model to run over"
  exit
fi

# set up the environment
unset MARLIN_DLL
source ~/clic/setup.sh
export PYTHONPATH=${LCIO}/src/python:${ROOTSYS}/lib
source $DD4HEP/bin/thisdd4hep.sh
 
# find the detector model specified
DETECTOR="$(find $ILCSOFT/lcgeo -name ${1}.xml)" 
echo Detector is $DETECTOR

# generate 1000 muons in random directions
#python lcio_particle_gun.py

# run geant4 with the detector model specified
#python lcio_geant4.py $DETECTOR

# run the CLIC reconstruction chain includin all tracking detectors
#echo $MARLIN_DLL
Marlin clicReconstruction.xml

