


SET( test_name "test_CLIC_o2_v04" )
ADD_TEST( t_${test_name} "$ENV{lcgeo_DIR}/bin/run_test_lcgeo.sh"
  ddsim --compactFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml --runType=batch -G -N=1 --outputFile=testCLIC_o2_v04.slcio --steeringFile=${CMAKE_SOURCE_DIR}/examples/clic_steer.py )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )
