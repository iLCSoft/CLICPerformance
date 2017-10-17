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
  ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml
  --global.LCIOInputFiles=testCLIC_o2_v04.slcio clicReconstruction_o2.xml
  --Config.Tracking=Truth
  --global.MaxRecordNumber=3
  --Output_REC.LCIOOutputFile=reco_truth_rec.slcio
  --Output_DST.LCIOOutputFile=reco_truth_dst.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ CONFIGURATIONS o2_v04
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )


SET( test_name "test_CLIC_o2_v04_reco_conformal" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_SOURCE_DIR}/Tests/bin/MarlinTests.sh
  ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml
  --global.LCIOInputFiles=testCLIC_o2_v04.slcio clicReconstruction_o2.xml
  --Config.Tracking=Conformal
  --global.MaxRecordNumber=3
  --Output_REC.LCIOOutputFile=reco_conformal_rec.slcio
  --Output_DST.LCIOOutputFile=reco_conformal_dst.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ CONFIGURATIONS o2_v04
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )

SET( test_name "test_CLIC_o2_v04_reco_overlay" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_SOURCE_DIR}/Tests/bin/MarlinTests.sh
  ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o2_v04/CLIC_o2_v04.xml
  --global.LCIOInputFiles=testCLIC_o2_v04.slcio clicReconstruction_o2.xml
  --Config.Overlay=3TeV
  --Overlay3TeV.BackgroundFileNames=testCLIC_o2_v04.slcio
  --global.MaxRecordNumber=3
  --Output_REC.LCIOOutputFile=reco_overlay_rec.slcio
  --Output_DST.LCIOOutputFile=reco_overlay_dst.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ CONFIGURATIONS o2_v04
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )


SET( test_name "test_CLIC_o2_v04_reco_conformal_overlay_valgrind" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_SOURCE_DIR}/Tests/bin/ValgrindTests.sh ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so ${VALGRIND_LIB}
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
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ CONFIGURATIONS o2_v04
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )




SET( test_name "test_CLIC_o3_v13_sim" )
ADD_TEST( NAME t_${test_name}
  COMMAND "$ENV{lcgeo_DIR}/bin/run_test_lcgeo.sh"
  ddsim
  --compactFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v13/CLIC_o3_v13.xml
  --runType=batch -G -N=100
  --outputFile=testCLIC_o3_v13.slcio
  --gun.isotrop=true
  --steeringFile=${CMAKE_SOURCE_DIR}/examples/clic_steer.py
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error" )

SET( test_name "test_CLIC_o3_v13_reco_truth" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_SOURCE_DIR}/Tests/bin/MarlinTests.sh
  ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v13/CLIC_o3_v13.xml
  --global.LCIOInputFiles=testCLIC_o3_v13.slcio clicReconstruction.xml
  --Config.Tracking=Truth
  --global.MaxRecordNumber=4
  --Output_REC.LCIOOutputFile=reco_truth_v13_rec.slcio
  --Output_DST.LCIOOutputFile=reco_truth_v13_dst.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )


SET( test_name "test_CLIC_o3_v13_reco_conformal" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_SOURCE_DIR}/Tests/bin/MarlinTests.sh
  ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v13/CLIC_o3_v13.xml
  --global.LCIOInputFiles=testCLIC_o3_v13.slcio clicReconstruction.xml
  --Config.Tracking=Conformal
  --global.MaxRecordNumber=20
  --Output_REC.LCIOOutputFile=reco_conformal_v13_rec.slcio
  --Output_DST.LCIOOutputFile=reco_conformal_v13_dst.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )

SET( test_name "test_CLIC_o3_v13_reco_overlay" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_SOURCE_DIR}/Tests/bin/MarlinTests.sh
  ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v13/CLIC_o3_v13.xml
  --global.LCIOInputFiles=testCLIC_o3_v13.slcio clicReconstruction.xml
  --Config.Overlay=3TeV
  --Overlay3TeV.BackgroundFileNames=testCLIC_o3_v13.slcio
  --global.MaxRecordNumber=20
  --Output_REC.LCIOOutputFile=reco_overlay_v13_rec.slcio
  --Output_DST.LCIOOutputFile=reco_overlay_v13_dst.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )

SET( test_name "test_CLIC_o3_v13_reco_overlay_valgrind" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_SOURCE_DIR}/Tests/bin/ValgrindTests.sh ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so ${VALGRIND_LIB}
  ${VALGRIND}
  --leak-check=full
  --show-leak-kinds=all
  --suppressions=${ROOT_DIR}/../etc/valgrind-root.supp
  --suppressions=${CMAKE_SOURCE_DIR}/etc/extra-root.supp
  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v13/CLIC_o3_v13.xml
  --global.LCIOInputFiles=testCLIC_o3_v13.slcio clicReconstruction.xml
  --Config.Overlay=3TeV
  --Overlay3TeV.BackgroundFileNames=testCLIC_o3_v13.slcio
  --global.MaxRecordNumber=3
  --Output_REC.LCIOOutputFile=reco_overlay_valgrind_v13_rec.slcio
  --Output_DST.LCIOOutputFile=reco_overlay_valgrind_v13_dst.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ CONFIGURATIONS valgrind
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error" )

SET( test_name "test_CLIC_o3_v13_reco_conformal_valgrind" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_SOURCE_DIR}/Tests/bin/ValgrindTests.sh ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so ${VALGRIND_LIB}
  ${VALGRIND}
  --leak-check=full
  --show-leak-kinds=all
  --suppressions=${ROOT_DIR}/../etc/valgrind-root.supp
  --suppressions=${CMAKE_SOURCE_DIR}/etc/extra-root.supp
  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v13/CLIC_o3_v13.xml
  --global.LCIOInputFiles=testCLIC_o3_v13.slcio clicReconstruction.xml
  --Config.Tracking=Conformal
  --global.MaxRecordNumber=3
  --Output_REC.LCIOOutputFile=reco_conformal_valgrind_v13_rec.slcio
  --Output_DST.LCIOOutputFile=reco_conformal_valgrind_v13_dst.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ CONFIGURATIONS valgrind
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  TIMEOUT 300000
  )


SET( test_name "test_CLIC_o3_v13_reco_conformal_overlay_valgrind" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_SOURCE_DIR}/Tests/bin/ValgrindTests.sh ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so ${VALGRIND_LIB}
  ${VALGRIND}
  --leak-check=full
  --show-leak-kinds=all
  --suppressions=${ROOT_DIR}/../etc/valgrind-root.supp
  --suppressions=${CMAKE_SOURCE_DIR}/etc/extra-root.supp
  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v13/CLIC_o3_v13.xml
  --global.LCIOInputFiles=testCLIC_o3_v13.slcio clicReconstruction.xml
  --Config.Overlay=3TeV
  --Config.Tracking=Conformal
  --Overlay3TeV.BackgroundFileNames=testCLIC_o3_v13.slcio
  --global.MaxRecordNumber=3
  --Output_REC.LCIOOutputFile=reco_overlay_conformal_valgrind_v13_rec.slcio
  --Output_DST.LCIOOutputFile=reco_overlay_conformal_valgrind_v13_dst.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ CONFIGURATIONS valgrind
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  TIMEOUT 300000
  )


SET( test_name "test_CLIC_o3_v13_sim_zuds" )
ADD_TEST( NAME t_${test_name}
  COMMAND "$ENV{lcgeo_DIR}/bin/run_test_lcgeo.sh"
  ddsim
  --compactFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v13/CLIC_o3_v13.xml
  --runType=batch -N=4
  --outputFile=testCLIC_o3_v13_zuds.slcio
  --inputFiles /eos/experiment/clicdp/grid/ilc/user/w/webermat/stdhep/Zuds3000/Zuds3000__01.stdhep
  --steeringFile=${CMAKE_SOURCE_DIR}/examples/clic_steer.py
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/
  CONFIGURATIONS EOS
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;ERROR;Error"
  REQUIRED_FILES /eos/experiment/clicdp/grid/ilc/user/w/webermat/stdhep/Zuds3000/Zuds3000__01.stdhep
  TIMEOUT 300000
  )

SET( test_name "test_CLIC_o3_v13_reco_conformal_zuds" )
ADD_TEST( NAME t_${test_name}
  COMMAND ${CMAKE_SOURCE_DIR}/Tests/bin/MarlinTests.sh
  ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v13/CLIC_o3_v13.xml
  --global.LCIOInputFiles=testCLIC_o3_v13_zuds.slcio clicReconstruction.xml
  --Config.Tracking=Conformal
  --global.MaxRecordNumber=4
  --Output_REC.LCIOOutputFile=reco_conformal_v13_zuds_rec.slcio
  --Output_DST.LCIOOutputFile=reco_conformal_v13_zuds_dst.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/
  CONFIGURATIONS EOS
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  DEPENDS test_CLIC_o3_v13_sim_zuds
  TIMEOUT 300000
  )


SET( test_name "test_CLIC_o3_v13_reco_conformal_zuds_valgrind" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_SOURCE_DIR}/Tests/bin/ValgrindTests.sh ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so ${VALGRIND_LIB}
  ${VALGRIND}
  --leak-check=full
  --show-leak-kinds=all
  --suppressions=${ROOT_DIR}/../etc/valgrind-root.supp
  --suppressions=${CMAKE_SOURCE_DIR}/etc/extra-root.supp
  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v13/CLIC_o3_v13.xml
  --global.LCIOInputFiles=testCLIC_o3_v13_zuds.slcio clicReconstruction.xml
  --Config.Tracking=Conformal
  --global.MaxRecordNumber=4
  --Output_REC.LCIOOutputFile=reco_conformal_valgrind_v13_zuds_rec.slcio
  --Output_DST.LCIOOutputFile=reco_conformal_valgrind_v13_zuds_dst.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ CONFIGURATIONS valgrind
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  DEPENDS test_CLIC_o3_v13_sim_zuds
  TIMEOUT 300000
  )

SET( test_name "test_CLIC_o3_v13_reco_conformal_overlay_zuds_valgrind" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_SOURCE_DIR}/Tests/bin/ValgrindTests.sh ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so ${VALGRIND_LIB}
  ${VALGRIND}
  --leak-check=full
  --show-leak-kinds=all
  --suppressions=${ROOT_DIR}/../etc/valgrind-root.supp
  --suppressions=${CMAKE_SOURCE_DIR}/etc/extra-root.supp
  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v13/CLIC_o3_v13.xml
  --global.LCIOInputFiles=testCLIC_o3_v13_zuds.slcio clicReconstruction.xml
  --Config.Overlay=3TeV
  --Config.Tracking=Conformal
  --Overlay3TeV.BackgroundFileNames="/eos/experiment/clicdp/grid/ilc/prod/clic/3tev/gghad/CLIC_o3_v13/SIM/00008258/000/gghad_sim_8258_37.slcio /eos/experiment/clicdp/grid/ilc/prod/clic/3tev/gghad/CLIC_o3_v13/SIM/00008258/000/gghad_sim_8258_38.slcio"
  --global.MaxRecordNumber=4
  --Output_REC.LCIOOutputFile=reco_overlay_conformal_valgrind_v13_zuds_rec.slcio
  --Output_DST.LCIOOutputFile=reco_overlay_conformal_valgrind_v13_zuds_dst.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ CONFIGURATIONS valgrind
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  DEPENDS test_CLIC_o3_v13_sim_zuds
  REQUIRED_FILES "/eos/experiment/clicdp/grid/ilc/prod/clic/3tev/gghad/CLIC_o3_v13/SIM/00008258/000/gghad_sim_8258_37.slcio /eos/experiment/clicdp/grid/ilc/prod/clic/3tev/gghad/CLIC_o3_v13/SIM/00008258/000/gghad_sim_8258_38.slcio"
  TIMEOUT 300000
  )



SET( test_name "test_CLIC_o3_v13_reco_conformal_zuds_massif" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_SOURCE_DIR}/Tests/bin/ValgrindTests.sh ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so ${VALGRIND_LIB}
  ${VALGRIND}
  --tool=massif
  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v13/CLIC_o3_v13.xml
  --global.LCIOInputFiles=testCLIC_o3_v13_zuds.slcio clicReconstruction.xml
  --Config.Tracking=Conformal
  --global.MaxRecordNumber=4
  --Output_REC.LCIOOutputFile=reco_conformal_msassif_v13_zuds_rec.slcio
  --Output_DST.LCIOOutputFile=reco_conformal_msassif_v13_zuds_dst.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ CONFIGURATIONS valgrind
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  DEPENDS test_CLIC_o3_v13_sim_zuds
  TIMEOUT 300000
  )


SET( test_name "test_CLIC_o3_v13_reco_conformal_overlay_zuds_massif" )
ADD_TEST( NAME t_${test_name}
  COMMAND
  ${CMAKE_SOURCE_DIR}/Tests/bin/ValgrindTests.sh ${CMAKE_SOURCE_DIR}/lib/libClicPerformance.so ${VALGRIND_LIB}
  ${VALGRIND}
  --tool=massif
  Marlin
  --InitDD4hep.DD4hepXMLFile=$ENV{lcgeo_DIR}/CLIC/compact/CLIC_o3_v13/CLIC_o3_v13.xml
  --global.LCIOInputFiles=testCLIC_o3_v13_zuds.slcio clicReconstruction.xml
  --Config.Overlay=3TeV
  --Config.Tracking=Conformal
  --Overlay3TeV.BackgroundFileNames="/eos/experiment/clicdp/grid/ilc/prod/clic/3tev/gghad/CLIC_o3_v13/SIM/00008258/000/gghad_sim_8258_37.slcio /eos/experiment/clicdp/grid/ilc/prod/clic/3tev/gghad/CLIC_o3_v13/SIM/00008258/000/gghad_sim_8258_38.slcio"
  --global.MaxRecordNumber=4
  --Output_REC.LCIOOutputFile=reco_overlay_conformal_v13_zuds_rec.slcio
  --Output_DST.LCIOOutputFile=reco_overlay_conformal_v13_zuds_dst.slcio
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/ CONFIGURATIONS valgrind
  )
SET_TESTS_PROPERTIES( t_${test_name} PROPERTIES
  FAIL_REGULAR_EXPRESSION  "Exception;EXCEPTION;\\[ ERROR;Error"
  DEPENDS test_CLIC_o3_v13_sim_zuds
  REQUIRED_FILES "/eos/experiment/clicdp/grid/ilc/prod/clic/3tev/gghad/CLIC_o3_v13/SIM/00008258/000/gghad_sim_8258_37.slcio /eos/experiment/clicdp/grid/ilc/prod/clic/3tev/gghad/CLIC_o3_v13/SIM/00008258/000/gghad_sim_8258_38.slcio"
  TIMEOUT 300000
  )
