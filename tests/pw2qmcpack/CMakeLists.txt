#############################################################
# Add tests to ctest
#############################################################

message("Adding tests for pw.x, pw2qmcpack.x and QMCPACK workflow")

IF (NOT TEST_MAX_PROCS )
SET ( TEST_MAX_PROCS 16 )
message ("TEST_MAX_PROCS was unset. Set to 16")
ENDIF ()

set(QMCAPP_ERR "")
if((NOT ${TEST_MAX_PROCS}) OR (${TEST_MAX_PROCS} STREQUAL ""))
    set(QMCAPP_ERR "${QMCAPP_ERR}  TEST_MAX_PROCS not set: '${TEST_MAX_PROCS}'\n")
endif()
#
#
if(NOT ${QMCAPP_ERR} STREQUAL "")
  message("${QMCAPP_ERR}  skipping tests")
else()
  INCLUDE( "${qmcpack_SOURCE_DIR}/CMake/run_qe.cmake" )
  SET(LAST_TEST_NAME "NONE")

  RUN_QE_TEST(qe-LiH-unpolarized ${CMAKE_SOURCE_DIR}/examples/solids/dft-inputs
              16 16 1 16 1 1 LiH )

  SOFTLINK_H5( qe-LiH-unpolarized-np-16-16-1-nk-16-1-1
               LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 LiH LiH-gamma.pwscf.h5 LAST_TEST_NAME )

  # These values are slightly different from those in tests/systems
  # because the h5 is generated in different ways.
  LIST(APPEND LIH_GAMMA_SCALARS "kinetic" "7.46326 0.024")
  LIST(APPEND LIH_GAMMA_SCALARS "localecp" "-11.78311 0.025")
  LIST(APPEND LIH_GAMMA_SCALARS "samples" "1600000 0.0")
  LIST(APPEND LIH_GAMMA_SCALARS "flux" "0.0 0.05")

  QMC_RUN_AND_CHECK(qe-LiH-unpolarized-np-16-16-1-nk-16-1-1-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-gamma
                    hf_vmc_LiH-gamma.xml
                    1 16
                    LIH_GAMMA_SCALARS
                    0 # VMC
                    TRUE)

  SET_TESTS_PROPERTIES( qe-LiH-unpolarized-np-16-16-1-nk-16-1-1-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 
                        PROPERTIES DEPENDS ${LAST_TEST_NAME} )

  SOFTLINK_H5( qe-LiH-unpolarized-np-16-16-1-nk-16-1-1
               LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 LiH LiH-x.pwscf.h5 LAST_TEST_NAME )

  # These values are slightly different from those in tests/systems
  # because the h5 is generated in different ways.
  LIST(APPEND LIH_X_SCALARS "kinetic" "7.33899 0.024")
  LIST(APPEND LIH_X_SCALARS "localecp" "-11.37256 0.025")
  LIST(APPEND LIH_X_SCALARS "samples" "1600000 0.0")
  LIST(APPEND LIH_X_SCALARS "flux" "0.0 0.05")

  QMC_RUN_AND_CHECK(qe-LiH-unpolarized-np-16-16-1-nk-16-1-1-LiH_solid_1x1x1_pp-x-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-x
                    hf_vmc_LiH-x.xml
                    1 16
                    LIH_X_SCALARS
                    0 # VMC
                    TRUE)

  SET_TESTS_PROPERTIES( qe-LiH-unpolarized-np-16-16-1-nk-16-1-1-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 
                        PROPERTIES DEPENDS ${LAST_TEST_NAME} )

  RUN_QE_TEST(qe-LiH-unpolarized ${CMAKE_SOURCE_DIR}/examples/solids/dft-inputs
              16 16 2 16 2 2 LiH )

  SOFTLINK_H5( qe-LiH-unpolarized-np-16-16-2-nk-16-2-2
               LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 LiH LiH-gamma.pwscf.h5 LAST_TEST_NAME )

  QMC_RUN_AND_CHECK(qe-LiH-unpolarized-np-16-16-2-nk-16-2-2-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-gamma
                    hf_vmc_LiH-gamma.xml
                    1 16
                    LIH_GAMMA_SCALARS
                    0 # VMC
                    TRUE)

  SET_TESTS_PROPERTIES( qe-LiH-unpolarized-np-16-16-2-nk-16-2-2-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 
                        PROPERTIES DEPENDS ${LAST_TEST_NAME} )

  SOFTLINK_H5( qe-LiH-unpolarized-np-16-16-2-nk-16-2-2
               LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 LiH LiH-x.pwscf.h5 LAST_TEST_NAME )

  QMC_RUN_AND_CHECK(qe-LiH-unpolarized-np-16-16-2-nk-16-2-2-LiH_solid_1x1x1_pp-x-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-x
                    hf_vmc_LiH-x.xml
                    1 16
                    LIH_X_SCALARS
                    0 # VMC
                    TRUE)

  SET_TESTS_PROPERTIES( qe-LiH-unpolarized-np-16-16-2-nk-16-2-2-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 
                        PROPERTIES DEPENDS ${LAST_TEST_NAME} )

  # test case for no wf_collect
  RUN_QE_TEST(qe-LiH-unpolarized-no-collect ${CMAKE_SOURCE_DIR}/examples/solids/dft-inputs-no-collect
              16 16 16 16 1 1 LiH )

  SOFTLINK_H5( qe-LiH-unpolarized-no-collect-np-16-16-16-nk-16-1-1
               LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 LiH LiH-gamma.pwscf.h5 LAST_TEST_NAME )

  QMC_RUN_AND_CHECK(qe-LiH-unpolarized-no-collect-np-16-16-16-nk-16-1-1-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-gamma
                    hf_vmc_LiH-gamma.xml
                    1 16
                    LIH_GAMMA_SCALARS
                    0 # VMC
                    TRUE)

  SET_TESTS_PROPERTIES( qe-LiH-unpolarized-no-collect-np-16-16-16-nk-16-1-1-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 
                        PROPERTIES DEPENDS ${LAST_TEST_NAME} )

  SOFTLINK_H5( qe-LiH-unpolarized-no-collect-np-16-16-16-nk-16-1-1
               LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 LiH LiH-x.pwscf.h5 LAST_TEST_NAME )

  QMC_RUN_AND_CHECK(qe-LiH-unpolarized-no-collect-np-16-16-16-nk-16-1-1-LiH_solid_1x1x1_pp-x-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-x
                    hf_vmc_LiH-x.xml
                    1 16
                    LIH_X_SCALARS
                    0 # VMC
                    TRUE)

  SET_TESTS_PROPERTIES( qe-LiH-unpolarized-no-collect-np-16-16-16-nk-16-1-1-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 
                        PROPERTIES DEPENDS ${LAST_TEST_NAME} )

  RUN_QE_TEST(qe-LiH-unpolarized-no-collect ${CMAKE_SOURCE_DIR}/examples/solids/dft-inputs-no-collect
              16 16 16 16 2 2 LiH )

  SOFTLINK_H5( qe-LiH-unpolarized-no-collect-np-16-16-16-nk-16-2-2
               LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 LiH LiH-gamma.pwscf.h5 LAST_TEST_NAME )

  QMC_RUN_AND_CHECK(qe-LiH-unpolarized-no-collect-np-16-16-16-nk-16-2-2-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-gamma
                    hf_vmc_LiH-gamma.xml
                    1 16
                    LIH_GAMMA_SCALARS
                    0 # VMC
                    TRUE)

  SET_TESTS_PROPERTIES( qe-LiH-unpolarized-no-collect-np-16-16-16-nk-16-2-2-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 
                        PROPERTIES DEPENDS ${LAST_TEST_NAME} )

  SOFTLINK_H5( qe-LiH-unpolarized-no-collect-np-16-16-16-nk-16-2-2
               LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 LiH LiH-x.pwscf.h5 LAST_TEST_NAME )

  QMC_RUN_AND_CHECK(qe-LiH-unpolarized-no-collect-np-16-16-16-nk-16-2-2-LiH_solid_1x1x1_pp-x-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-x
                    hf_vmc_LiH-x.xml
                    1 16
                    LIH_X_SCALARS
                    0 # VMC
                    TRUE)

  SET_TESTS_PROPERTIES( qe-LiH-unpolarized-no-collect-np-16-16-16-nk-16-2-2-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 
                        PROPERTIES DEPENDS ${LAST_TEST_NAME} )

endif()
