#!/bin/bash
#
# Simple script to run marlin tests with profiler on top
#   calls the command (given as second argument)
#   with all following arguments

#----- initialize profiler environment
source /cvmfs/projects.cern.ch/intelsw/psxe/linux/all-setup.sh
#----- parse command line - first argument is the
#      test to run
## Drop existing ClicPerformance from MARLIN_DLL
MARLIN_DLL=`echo $MARLIN_DLL | sed -E 's/(:?)[^:]*libClicPerformance\.so[.\d]*:?/\1/'`
export MARLIN_DLL=$1:$MARLIN_DLL
command=$2
theargs=""

##Drop two
shift
shift
for i in "$@" ; do
    theargs="${theargs} $i"
done

echo " #### MARLIN_DLL = :  ${MARLIN_DLL}"

echo " ### running test :  '${command} ${theargs}'"
exec ${command} ${theargs}
