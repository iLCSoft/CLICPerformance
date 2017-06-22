#!/bin/bash
#
# Simple script to run marlin tests with valgrind on top
#   calls the command (given as first argument)
#   with all following arguments
#

#----- initialize environment for this package - including DD4hep

#----- parse command line - first argument is the
#      test to run
## Drop existing ClicPerformance from MARLIN_DLL
MARLIN_DLL=`echo $MARLIN_DLL | sed -e 's/:\?[-a-z0-9A-Z_/.]*libClicPerformance.so[.0-9]*:\?/:/'`
export MARLIN_DLL=$1:$MARLIN_DLL
export VALGRIND_LIB=$2
command=$3
theargs=""

##Drop three
shift
shift
shift
for i in "$@" ; do
    theargs="${theargs} $i"
done

echo " #### MARLIN_DLL = :  ${MARLIN_DLL}"

echo " ### running test :  '${command} ${theargs}'"
exec ${command} ${theargs}
