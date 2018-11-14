SET( test_name "test_CLIC_o2_v04_sim" )
ADD_TEST( NAME t_${test_name}
  COMMAND "run_test_lcgeo.sh"
  ddsim --compactFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml --runType=batch -G -N=100 --outputFile=testCLIC_o2_v04.slcio --gun.isotrop=true --steeringFile=${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/clic_steer.py
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/ CONFIGURATIONS o2_v04
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_CLIC_o2_v04_reco_truth" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/MarlinTests.sh
  ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml
  --global.LCIOInputFiles=testCLIC_o2_v04.slcio clicReconstruction_o2.xml
  --Config.Tracking=Truth
  --global.MaxRecordNumber=3
  --Output_REC.LCIOOutputFile=reco_truth_rec.slcio
  --Output_DST.LCIOOutputFile=reco_truth_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/ CONFIGURATIONS o2_v04
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )


SET( test_name "test_CLIC_o2_v04_reco_conformal" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/MarlinTests.sh
  ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml
  --global.LCIOInputFiles=testCLIC_o2_v04.slcio clicReconstruction_o2.xml
  --Config.Tracking=Conformal
  --global.MaxRecordNumber=3
  --Output_REC.LCIOOutputFile=reco_conformal_rec.slcio
  --Output_DST.LCIOOutputFile=reco_conformal_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/ CONFIGURATIONS o2_v04
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )

SET( test_name "test_CLIC_o2_v04_reco_overlay" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/MarlinTests.sh
  ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml
  --global.LCIOInputFiles=testCLIC_o2_v04.slcio clicReconstruction_o2.xml
  --Config.Overlay=3TeV
  --Overlay3TeV.BackgroundFileNames=testCLIC_o2_v04.slcio
  --global.MaxRecordNumber=3
  --Output_REC.LCIOOutputFile=reco_overlay_rec.slcio
  --Output_DST.LCIOOutputFile=reco_overlay_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/ CONFIGURATIONS o2_v04
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )


SET( test_name "test_CLIC_o2_v04_reco_conformal_overlay_valgrind" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/ValgrindTests.sh ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so ${VALGRIND_LIB}
  ${VALGRIND}
  --leak-check=full
  --show-leak-kinds=all
  --suppressions=${ROOT_DIR}/../etc/valgrind-root.supp Marlin
  --InitDD4hep.DD4hepXMLFile=CLIC_o2_v04/CLIC_o2_v04.xml
  --global.LCIOInputFiles=testCLIC_o2_v04.slcio clicReconstruction_o2.xml
  --Config.Overlay=3TeV
  --Config.Tracking=Conformal
  --Overlay3TeV.BackgroundFileNames=testCLIC_o2_v04.slcio
  --global.MaxRecordNumber=10
  --Output_REC.LCIOOutputFile=reco_overlay_conform_valgrind_v04_rec.slcio
  --Output_DST.LCIOOutputFile=reco_overlay_conform_valgrind_v04_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/ CONFIGURATIONS o2_v04
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )



SET(EOS_BACKGROUND_FILE /eos/experiment/clicdp/grid/ilc/prod/clic/3tev/gghad/CLIC_o3_v14/SIM/00009245/000/gghad_sim_9245_37.slcio)
SET(CLIC_DETECTOR_MODEL CLIC_o3_v14)
SET( test_name "test_${CLIC_DETECTOR_MODEL}_sim" )
ADD_TEST( NAME t_${test_name}
  COMMAND "run_test_lcgeo.sh"
  ddsim
  --compactFile=$ENV{lcgeo_DIR}/CLIC/compact/${CLIC_DETECTOR_MODEL}/${CLIC_DETECTOR_MODEL}.xml
  --runType=batch --inputFile ../Tests/yyxyev_000.stdhep -N 2
  --outputFile=test${CLIC_DETECTOR_MODEL}.slcio
  --gun.isotrop=true
  --steeringFile=${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/clic_steer.py
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_${CLIC_DETECTOR_MODEL}_reco_truth" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/MarlinTests.sh
  ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/${CLIC_DETECTOR_MODEL}/${CLIC_DETECTOR_MODEL}.xml
  --global.LCIOInputFiles=test${CLIC_DETECTOR_MODEL}.slcio clicReconstruction.xml
  --Config.Tracking=Truth
  --global.MaxRecordNumber=2
  --Output_REC.LCIOOutputFile=reco_truth_${CLIC_DETECTOR_MODEL}_rec.slcio
  --Output_DST.LCIOOutputFile=reco_truth_${CLIC_DETECTOR_MODEL}_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/
  CONFIGURATIONS TRUTH
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )


SET( test_name "test_${CLIC_DETECTOR_MODEL}_reco_conformal" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/MarlinTests.sh
  ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/${CLIC_DETECTOR_MODEL}/${CLIC_DETECTOR_MODEL}.xml
  --global.LCIOInputFiles=test${CLIC_DETECTOR_MODEL}.slcio clicReconstruction.xml
  --Config.Tracking=Conformal
  --global.MaxRecordNumber=2
  --Output_REC.LCIOOutputFile=reco_conformal_${CLIC_DETECTOR_MODEL}_rec.slcio
  --Output_DST.LCIOOutputFile=reco_conformal_${CLIC_DETECTOR_MODEL}_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )

SET( test_name "test_${CLIC_DETECTOR_MODEL}_reco_real_overlay" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/MarlinTests.sh
  ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/${CLIC_DETECTOR_MODEL}/${CLIC_DETECTOR_MODEL}.xml
  --global.LCIOInputFiles=test${CLIC_DETECTOR_MODEL}.slcio clicReconstruction.xml
  --Config.Tracking=Truth
  --Config.Overlay=3TeV
  --Overlay3TeV.BackgroundFileNames=${EOS_BACKGROUND_FILE}
  --global.MaxRecordNumber=4
  --Output_REC.LCIOOutputFile=reco_overlay_${CLIC_DETECTOR_MODEL}_rec.slcio
  --Output_DST.LCIOOutputFile=reco_overlay_${CLIC_DETECTOR_MODEL}_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/
  CONFIGURATIONS EOS
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  REQUIRED_FILES ${EOS_BACKGROUND_FILE}
)


SET( test_name "test_${CLIC_DETECTOR_MODEL}_reco_overlay_valgrind" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/ValgrindTests.sh ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so ${VALGRIND_LIB}
  ${VALGRIND}
  --leak-check=full
  --show-leak-kinds=all
  --suppressions=${ROOT_DIR}/../etc/valgrind-root.supp
  --suppressions=${CMAKE_CURRENT_SOURCE_DIR}/../etc/extra-root.supp
  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/${CLIC_DETECTOR_MODEL}/${CLIC_DETECTOR_MODEL}.xml
  --global.LCIOInputFiles=test${CLIC_DETECTOR_MODEL}.slcio clicReconstruction.xml
  --Config.Tracking=Truth
  --Config.Overlay=3TeV
  --Overlay3TeV.BackgroundFileNames=test${CLIC_DETECTOR_MODEL}.slcio
  --global.MaxRecordNumber=3
  --Output_REC.LCIOOutputFile=reco_overlay_valgrind_${CLIC_DETECTOR_MODEL}_rec.slcio
  --Output_DST.LCIOOutputFile=reco_overlay_valgrind_${CLIC_DETECTOR_MODEL}_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/ CONFIGURATIONS valgrind
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )

SET( test_name "test_${CLIC_DETECTOR_MODEL}_reco_conformal_valgrind" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/ValgrindTests.sh ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so ${VALGRIND_LIB}
  ${VALGRIND}
  --leak-check=full
  --show-leak-kinds=all
  --suppressions=${ROOT_DIR}/../etc/valgrind-root.supp
  --suppressions=${CMAKE_CURRENT_SOURCE_DIR}/../etc/extra-root.supp
  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/${CLIC_DETECTOR_MODEL}/${CLIC_DETECTOR_MODEL}.xml
  --global.LCIOInputFiles=test${CLIC_DETECTOR_MODEL}.slcio clicReconstruction.xml
  --Config.Tracking=Conformal
  --global.MaxRecordNumber=3
  --Output_REC.LCIOOutputFile=reco_conformal_valgrind_${CLIC_DETECTOR_MODEL}_rec.slcio
  --Output_DST.LCIOOutputFile=reco_conformal_valgrind_${CLIC_DETECTOR_MODEL}_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/ CONFIGURATIONS valgrind
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  TIMEOUT 300000
  )


SET( test_name "test_${CLIC_DETECTOR_MODEL}_reco_conformal_overlay_valgrind" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/ValgrindTests.sh ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so ${VALGRIND_LIB}
  ${VALGRIND}
  --leak-check=full
  --show-leak-kinds=all
  --suppressions=${ROOT_DIR}/../etc/valgrind-root.supp
  --suppressions=${CMAKE_CURRENT_SOURCE_DIR}/../etc/extra-root.supp
  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/${CLIC_DETECTOR_MODEL}/${CLIC_DETECTOR_MODEL}.xml
  --global.LCIOInputFiles=test${CLIC_DETECTOR_MODEL}.slcio clicReconstruction.xml
  --Config.Tracking=Conformal
  --Config.Overlay=3TeV
  --Overlay3TeV.BackgroundFileNames=test${CLIC_DETECTOR_MODEL}.slcio
  --global.MaxRecordNumber=3
  --Output_REC.LCIOOutputFile=reco_overlay_conformal_valgrind_${CLIC_DETECTOR_MODEL}_rec.slcio
  --Output_DST.LCIOOutputFile=reco_overlay_conformal_valgrind_${CLIC_DETECTOR_MODEL}_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/ CONFIGURATIONS valgrind
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  TIMEOUT 300000
  )


SET( test_name "test_${CLIC_DETECTOR_MODEL}_sim_zuds" )
ADD_TEST( NAME t_${test_name}
  COMMAND "run_test_lcgeo.sh"
  ddsim
  --compactFile=$ENV{lcgeo_DIR}/CLIC/compact/${CLIC_DETECTOR_MODEL}/${CLIC_DETECTOR_MODEL}.xml
  --runType=batch -N=4
  --outputFile=test${CLIC_DETECTOR_MODEL}_zuds.slcio
  --inputFiles /eos/experiment/clicdp/grid/ilc/user/w/webermat/stdhep/Zuds3000/Zuds3000__01.stdhep
  --steeringFile=${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/clic_steer.py
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/
  CONFIGURATIONS EOS
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error"
  REQUIRED_FILES /eos/experiment/clicdp/grid/ilc/user/w/webermat/stdhep/Zuds3000/Zuds3000__01.stdhep
  TIMEOUT 300000
  )

SET( test_name "test_${CLIC_DETECTOR_MODEL}_reco_conformal_zuds" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/MarlinTests.sh
  ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/${CLIC_DETECTOR_MODEL}/${CLIC_DETECTOR_MODEL}.xml
  --global.LCIOInputFiles=test${CLIC_DETECTOR_MODEL}_zuds.slcio clicReconstruction.xml
  --Config.Tracking=Conformal
  --global.MaxRecordNumber=4
  --Output_REC.LCIOOutputFile=reco_conformal_${CLIC_DETECTOR_MODEL}_zuds_rec.slcio
  --Output_DST.LCIOOutputFile=reco_conformal_${CLIC_DETECTOR_MODEL}_zuds_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/
  CONFIGURATIONS EOS
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  DEPENDS test_${CLIC_DETECTOR_MODEL}_sim_zuds
  TIMEOUT 300000
  )

SET( test_name "test_${CLIC_DETECTOR_MODEL}_reco_conformal_zuds_overlay" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/MarlinTests.sh
  ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/${CLIC_DETECTOR_MODEL}/${CLIC_DETECTOR_MODEL}.xml
  --global.LCIOInputFiles=test${CLIC_DETECTOR_MODEL}_zuds.slcio clicReconstruction.xml
  --Config.Tracking=Conformal
  --global.MaxRecordNumber=4
  --Config.Overlay=3TeV
  --Overlay3TeV.BackgroundFileNames=${EOS_BACKGROUND_FILE}
  --Output_REC.LCIOOutputFile=reco_conformal_${CLIC_DETECTOR_MODEL}_zuds_rec.slcio
  --Output_DST.LCIOOutputFile=reco_conformal_${CLIC_DETECTOR_MODEL}_zuds_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/
  CONFIGURATIONS EOS
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  DEPENDS test_${CLIC_DETECTOR_MODEL}_sim_zuds
  REQUIRED_FILES "${EOS_BACKGROUND_FILE}"
  TIMEOUT 300000
  )



SET( test_name "test_${CLIC_DETECTOR_MODEL}_reco_conformal_zuds_valgrind" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/ValgrindTests.sh ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so ${VALGRIND_LIB}
  ${VALGRIND}
  --leak-check=full
  --show-leak-kinds=all
  --suppressions=${ROOT_DIR}/../etc/valgrind-root.supp
  --suppressions=${CMAKE_CURRENT_SOURCE_DIR}/../etc/extra-root.supp
  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/${CLIC_DETECTOR_MODEL}/${CLIC_DETECTOR_MODEL}.xml
  --global.LCIOInputFiles=test${CLIC_DETECTOR_MODEL}_zuds.slcio clicReconstruction.xml
  --Config.Tracking=Conformal
  --global.MaxRecordNumber=4
  --Output_REC.LCIOOutputFile=reco_conformal_valgrind_${CLIC_DETECTOR_MODEL}_zuds_rec.slcio
  --Output_DST.LCIOOutputFile=reco_conformal_valgrind_${CLIC_DETECTOR_MODEL}_zuds_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/ CONFIGURATIONS valgrind
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  DEPENDS test_${CLIC_DETECTOR_MODEL}_sim_zuds
  TIMEOUT 300000
  )

SET( test_name "test_${CLIC_DETECTOR_MODEL}_reco_conformal_overlay_zuds_valgrind" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/ValgrindTests.sh ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so ${VALGRIND_LIB}
  ${VALGRIND}
  --leak-check=full
  --show-leak-kinds=all
  --suppressions=${ROOT_DIR}/../etc/valgrind-root.supp
  --suppressions=${CMAKE_CURRENT_SOURCE_DIR}/../etc/extra-root.supp
  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/${CLIC_DETECTOR_MODEL}/${CLIC_DETECTOR_MODEL}.xml
  --global.LCIOInputFiles=test${CLIC_DETECTOR_MODEL}_zuds.slcio clicReconstruction.xml
  --Config.Tracking=Conformal
  --Config.Overlay=3TeV
  --Overlay3TeV.BackgroundFileNames=${EOS_BACKGROUND_FILE}
  --global.MaxRecordNumber=4
  --Output_REC.LCIOOutputFile=reco_overlay_conformal_valgrind_${CLIC_DETECTOR_MODEL}_zuds_rec.slcio
  --Output_DST.LCIOOutputFile=reco_overlay_conformal_valgrind_${CLIC_DETECTOR_MODEL}_zuds_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/ CONFIGURATIONS valgrind
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  DEPENDS test_${CLIC_DETECTOR_MODEL}_sim_zuds
  REQUIRED_FILES ${EOS_BACKGROUND_FILE}
  TIMEOUT 300000
  )



SET( test_name "test_${CLIC_DETECTOR_MODEL}_reco_conformal_zuds_massif" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/ValgrindTests.sh ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so ${VALGRIND_LIB}
  ${VALGRIND}
  --tool=massif
  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/${CLIC_DETECTOR_MODEL}/${CLIC_DETECTOR_MODEL}.xml
  --global.LCIOInputFiles=test${CLIC_DETECTOR_MODEL}_zuds.slcio clicReconstruction.xml
  --Config.Tracking=Conformal
  --global.MaxRecordNumber=4
  --Output_REC.LCIOOutputFile=reco_conformal_msassif_${CLIC_DETECTOR_MODEL}_zuds_rec.slcio
  --Output_DST.LCIOOutputFile=reco_conformal_msassif_${CLIC_DETECTOR_MODEL}_zuds_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/ CONFIGURATIONS valgrind
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  DEPENDS test_${CLIC_DETECTOR_MODEL}_sim_zuds
  TIMEOUT 300000
  )


SET( test_name "test_${CLIC_DETECTOR_MODEL}_reco_conformal_overlay_zuds_massif" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/ValgrindTests.sh ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so ${VALGRIND_LIB}
  ${VALGRIND}
  --tool=massif
  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/${CLIC_DETECTOR_MODEL}/${CLIC_DETECTOR_MODEL}.xml
  --global.LCIOInputFiles=test${CLIC_DETECTOR_MODEL}_zuds.slcio clicReconstruction.xml
  --Config.Tracking=Conformal
  --Config.Overlay=3TeV
  --Overlay3TeV.BackgroundFileNames=${EOS_BACKGROUND_FILE}
  --global.MaxRecordNumber=4
  --Output_REC.LCIOOutputFile=reco_overlay_conformal_${CLIC_DETECTOR_MODEL}_zuds_rec.slcio
  --Output_DST.LCIOOutputFile=reco_overlay_conformal_${CLIC_DETECTOR_MODEL}_zuds_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/ CONFIGURATIONS valgrind
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  DEPENDS test_${CLIC_DETECTOR_MODEL}_sim_zuds
  REQUIRED_FILES ${EOS_BACKGROUND_FILE}
  TIMEOUT 300000
  )


SET( test_name "test_${CLIC_DETECTOR_MODEL}_reco_conformal_zuds_profile" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/Profiler.sh ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so
  amplxe-cl -collect hotspots -r ../ClicReco/r@@@{at}
  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/${CLIC_DETECTOR_MODEL}/${CLIC_DETECTOR_MODEL}.xml
  --global.LCIOInputFiles=test${CLIC_DETECTOR_MODEL}_zuds.slcio clicReconstruction.xml
  --Config.Tracking=Conformal
  --global.MaxRecordNumber=10
  --Output_REC.LCIOOutputFile=reco_conformal_${CLIC_DETECTOR_MODEL}_zuds_rec_prof.slcio
  --Output_DST.LCIOOutputFile=reco_conformal_${CLIC_DETECTOR_MODEL}_zuds_dst_prof.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/
  CONFIGURATIONS PROF
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  DEPENDS test_${CLIC_DETECTOR_MODEL}_sim_zuds
  TIMEOUT 300000
  )


SET( test_name "test_${CLIC_DETECTOR_MODEL}_reco_conformal_zuds_overlay_profile" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/Profiler.sh ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so
  amplxe-cl -collect hotspots -r ../ClicReco/r@@@{at}_over
  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/${CLIC_DETECTOR_MODEL}/${CLIC_DETECTOR_MODEL}.xml
  --global.LCIOInputFiles=test${CLIC_DETECTOR_MODEL}_zuds.slcio clicReconstruction.xml
  --Config.Tracking=Conformal
  --Config.Overlay=3TeV
  --Overlay3TeV.BackgroundFileNames=${EOS_BACKGROUND_FILE}
  --global.MaxRecordNumber=1
  --Output_REC.LCIOOutputFile=reco_overlay_conformal_${CLIC_DETECTOR_MODEL}_zuds_rec_over_prof.slcio
  --Output_DST.LCIOOutputFile=reco_overlay_conformal_${CLIC_DETECTOR_MODEL}_zuds_dst_over_prof.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../clicConfig/
  CONFIGURATIONS PROF
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  DEPENDS test_${CLIC_DETECTOR_MODEL}_sim_zuds
  REQUIRED_FILES ${EOS_BACKGROUND_FILE}
  TIMEOUT 300000
  )

SET(FCCee_DETECTOR_MODEL FCCee_o1_v04)
SET( test_name "test_${FCCee_DETECTOR_MODEL}_sim" )
ADD_TEST( NAME t_${test_name}
  COMMAND "run_test_lcgeo.sh"
  ddsim --compactFile=$ENV{lcgeo_DIR}/FCCee/compact/${FCCee_DETECTOR_MODEL}/${FCCee_DETECTOR_MODEL}.xml
  --runType=batch --inputFiles ../Tests/yyxyev_000.stdhep -N 2
  --outputFile=test${FCCee_DETECTOR_MODEL}.slcio --gun.isotrop=true --steeringFile=${CMAKE_CURRENT_SOURCE_DIR}/../fcceeConfig/fcc_steer.py
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../fcceeConfig/
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )
 
SET( test_name "test_${FCCee_DETECTOR_MODEL}_reco_truth" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/MarlinTests.sh
  ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/FCCee/compact/${FCCee_DETECTOR_MODEL}/${FCCee_DETECTOR_MODEL}.xml
  --global.LCIOInputFiles=test${FCCee_DETECTOR_MODEL}.slcio fccReconstruction.xml
  --Config.Tracking=Truth
  --global.MaxRecordNumber=3
  --Output_REC.LCIOOutputFile=reco_truth_rec.slcio
  --Output_DST.LCIOOutputFile=reco_truth_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../fcceeConfig/
  CONFIGURATIONS TRUTH
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )
 
SET( test_name "test_${FCCee_DETECTOR_MODEL}_reco_conformal" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../Tests/bin/MarlinTests.sh
  ${CMAKE_INSTALL_PREFIX}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/FCCee/compact/${FCCee_DETECTOR_MODEL}/${FCCee_DETECTOR_MODEL}.xml
  --global.LCIOInputFiles=test${FCCee_DETECTOR_MODEL}.slcio fccReconstruction.xml
  --Config.Tracking=Conformal
  --global.MaxRecordNumber=3
  --Output_REC.LCIOOutputFile=reco_conformal_rec.slcio
  --Output_DST.LCIOOutputFile=reco_conformal_dst.slcio
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../fcceeConfig/
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )
