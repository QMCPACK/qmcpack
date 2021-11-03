#############################################################
# Add tests to ctest
#############################################################

message("Adding RMG tests")

if(NOT TEST_MAX_PROCS)
  set(TEST_MAX_PROCS 16)
  message("TEST_MAX_PROCS was unset. Set to 16")
endif()

set(QMCAPP_ERR "")
if((NOT ${TEST_MAX_PROCS}) OR (${TEST_MAX_PROCS} STREQUAL ""))
  set(QMCAPP_ERR "${QMCAPP_ERR}  TEST_MAX_PROCS not set: '${TEST_MAX_PROCS}'\n")
endif()
#
#
if(NOT ${QMCAPP_ERR} STREQUAL "")
  message("${QMCAPP_ERR}  skipping tests")
else()
  include("${qmcpack_SOURCE_DIR}/CMake/run_rmg.cmake")
  set(LAST_TEST_NAME "NONE")
  set(RMG_TEST_NAME "NONE")

  set(RMG_SYSTEM "Diamond2-1x1x1-gamma-ccECP")
  run_rmg_test(
    rmg-${RMG_SYSTEM}
    ${qmcpack_SOURCE_DIR}/examples/solids/rmg-inputs/${RMG_SYSTEM}
    1
    12
    RMG_TEST_NAME)
  set(LAST_TEST_NAME "${RMG_TEST_NAME}")


  # RMG_TEST_NAME is now rmg-Diamond2-1x1x1-gamma-ccECP-np-1
  # workdir is CMAKE_CURRENT_BINARY_DIR/RMG_TEST_NAME
  MESSAGE("${RMG_TEST_NAME}")
  
  softlink_rmg_input( ${RMG_TEST_NAME} 1-16 rmg-${RMG_SYSTEM}  LAST_TEST_NAME)
  # TODO: do larger run to get better reference?
  # ref values:                                                        adjusted for 1.6E4 samples
  # TotalSamples          =         120000000 +/-          0
  # LocalEnergy           =         -10.19480 +/-          0.00019     0.017
  # Kinetic               =          11.78992 +/-          0.00049     0.042
  # LocalPotential        =         -21.98472 +/-          0.00047     0.041
  # LocalECP              =          -7.34003 +/-          0.00051     0.044
  # NonLocalECP           =           0.41632 +/-          0.00018     0.016
  # ElecElec              =         -2.235309 +/-          0.000149    0.013
  # MPC                   =          -1.93438 +/-          0.00016     0.014
  # Variance              =            4.4782 +/-          0.0063
  list(APPEND DIAMOND2_GAMMA_SCALARS "kinetic"         "11.78992 0.042")
  #list(APPEND DIAMOND2_GAMMA_SCALARS "localpotential" "-21.98472 0.041")
  list(APPEND DIAMOND2_GAMMA_SCALARS "localecp"        "-7.34003 0.044")
  list(APPEND DIAMOND2_GAMMA_SCALARS "nonlocalecp"      "0.41632 0.016")
  list(APPEND DIAMOND2_GAMMA_SCALARS "samples" "16000 0.0")
  
  #set(THIS_TEST_NAME "${RMG_TEST_NAME}-Diamond2_1x1x1-gamma-ccECP-vmc_noj")
  set(THIS_TEST_NAME "${RMG_TEST_NAME}")
    qmc_run_and_check(
      ${THIS_TEST_NAME}
      "${qmcpack_SOURCE_DIR}/examples/solids/rmg-inputs/${RMG_SYSTEM}"
      rmg-Diamond2-1x1x1-gamma-ccECP
      rmg-Diamond2-1x1x1-gamma-ccECP.vmc.xml
      1
      16
      TRUE
      0
      DIAMOND2_GAMMA_SCALARS # VMC
    )
  
    if(${THIS_TEST_NAME})
      set_tests_properties(${THIS_TEST_NAME}-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})
    endif()
  

endif()
