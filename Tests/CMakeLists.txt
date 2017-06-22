


SET( test_name "test_CLIC_o2_v04_sim" )
ADD_TEST( NAME t_${test_name}
  COMMAND "$ENV{lcgeo_DIR}/bin/run_test_lcgeo.sh"
  ddsim --compactFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml --runType=batch -G -N=100 --outputFile=testCLIC_o2_v04.slcio --gun.isotrop=true --steeringFile=${CMAKE_SOURCE_DIR}/examples/clic_steer.py
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ CONFIGURATIONS o2_v04
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_CLIC_o2_v04_reco_truth" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_SOURCE_DIR}/Tests/bin/MarlinTests.sh
  ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so  Marlin --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml --global.LCIOInputFiles=testCLIC_o2_v04.slcio clicReconstruction_o2.xml --Config.Tracking=Truth --global.MaxRecordNumber=3 --MyLCIOOutputProcessor.LCIOOutputFile=reco_truth.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ CONFIGURATIONS o2_v04
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )


SET( test_name "test_CLIC_o2_v04_reco_conformal" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_SOURCE_DIR}/Tests/bin/MarlinTests.sh
  ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so  Marlin --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml --global.LCIOInputFiles=testCLIC_o2_v04.slcio clicReconstruction_o2.xml --Config.Tracking=Conformal --global.MaxRecordNumber=3  --MyLCIOOutputProcessor.LCIOOutputFile=reco_conformal.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ CONFIGURATIONS o2_v04
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )

SET( test_name "test_CLIC_o2_v04_reco_overlay" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_SOURCE_DIR}/Tests/bin/MarlinTests.sh
  ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so  Marlin --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml --global.LCIOInputFiles=testCLIC_o2_v04.slcio clicReconstruction_o2.xml --Config.Overlay=3TeV --Overlay3TeV.BackgroundFileNames=testCLIC_o2_v04.slcio --global.MaxRecordNumber=3  --MyLCIOOutputProcessor.LCIOOutputFile=reco_overlay.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ CONFIGURATIONS o2_v04
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )


SET( test_name "test_CLIC_o2_v04_reco_conformal_overlay_valgrind" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_SOURCE_DIR}/Tests/bin/MarlinTests.sh
  ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so  ${VALGRIND} --leak-check=full --show-leak-kinds=all --suppressions=${ROOT_DIR}/../etc/valgrind-root.supp Marlin --InitDD4hep.DD4hepXMLFile=CLIC_o2_v04/CLIC_o2_v04.xml --global.LCIOInputFiles=testCLIC_o2_v04.slcio clicReconstruction_o2.xml --Config.Overlay=3TeV --Config.Tracking=Conformal --Overlay3TeV.BackgroundFileNames=testCLIC_o2_v04.slcio --global.MaxRecordNumber=10  --MyLCIOOutputProcessor.LCIOOutputFile=reco_overlay_truth_v04.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ CONFIGURATIONS o2_v04
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )




SET( test_name "test_CLIC_o3_v11_sim" )
ADD_TEST( NAME t_${test_name}
  COMMAND "$ENV{lcgeo_DIR}/bin/run_test_lcgeo.sh"
  ddsim --compactFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v11/CLIC_o3_v11.xml --runType=batch -G -N=20 --outputFile=testCLIC_o3_v11.slcio --gun.isotrop=true --steeringFile=${CMAKE_SOURCE_DIR}/examples/clic_steer.py
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_CLIC_o3_v11_reco_truth" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_SOURCE_DIR}/Tests/bin/MarlinTests.sh
  ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so  Marlin --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v11/CLIC_o3_v11.xml --global.LCIOInputFiles=testCLIC_o3_v11.slcio clicReconstruction.xml --Config.Tracking=Truth
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )


SET( test_name "test_CLIC_o3_v11_reco_conformal" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_SOURCE_DIR}/Tests/bin/MarlinTests.sh
  ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so  Marlin --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v11/CLIC_o3_v11.xml --global.LCIOInputFiles=testCLIC_o3_v11.slcio clicReconstruction.xml --Config.Tracking=Conformal
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )

SET( test_name "test_CLIC_o3_v11_reco_overlay" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_SOURCE_DIR}/Tests/bin/MarlinTests.sh
  ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so  Marlin --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v11/CLIC_o3_v11.xml --global.LCIOInputFiles=testCLIC_o3_v11.slcio clicReconstruction.xml --Config.Overlay=3TeV --Overlay3TeV.BackgroundFileNames=testCLIC_o3_v11.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )
