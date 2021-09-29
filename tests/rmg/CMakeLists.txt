#############################################################
# Add tests to ctest
#############################################################

message("Adding tests for rmg -> convert4qmc -> QMCPACK workflow")

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
  # TODO: do larger run to get better reference
  list(APPEND DIAMOND2_GAMMA_SCALARS "kinetic" "11.78992 0.05")
  list(APPEND DIAMOND2_GAMMA_SCALARS "localecp" "-7.34003 0.05")
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
