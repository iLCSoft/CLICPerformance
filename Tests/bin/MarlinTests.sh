#!/bin/bash
#
# Simple script to run marlin tests
#   calls the command (given as first argument)
#   with all following arguments
#

#----- initialize environment for this package - including DD4hep

#----- parse command line - first argument is the
#      test to run
## Drop existing ClicPerformance from MARLIN_DLL
MARLIN_DLL=`echo $MARLIN_DLL | sed -E 's/(:?)[^:]*libClicPerformance\.so[.\d]*:?/\1/'`
export MARLIN_DLL=$1:$MARLIN_DLL
command=$2
theargs=""

shift
shift
for i in "$@" ; do
    theargs="${theargs} $i"
done

echo " #### MARLIN_DLL = :  ${MARLIN_DLL}"

echo " ### running test :  '${command} ${theargs}'"
exec ${command} ${theargs}
