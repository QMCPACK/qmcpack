include(test_labels)

# Runs unit tests
function(ADD_UNIT_TEST TESTNAME PROCS THREADS TEST_BINARY)
  message(VERBOSE "Adding test ${TESTNAME}")
  math(EXPR TOT_PROCS "${PROCS} * ${THREADS}")
  if(HAVE_MPI)
    add_test(NAME ${TESTNAME} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS}
                                      ${TEST_BINARY} ${ARGN})
    set(TEST_ADDED TRUE)
  else()
    if((${PROCS} STREQUAL "1"))
      add_test(NAME ${TESTNAME} COMMAND ${TEST_BINARY} ${ARGN})
      set(TEST_ADDED TRUE)
    else()
      message(VERBOSE "Disabling test ${TESTNAME} (building without MPI)")
    endif()
  endif()

  if(TEST_ADDED)
    set_tests_properties(${TESTNAME} PROPERTIES PROCESSORS ${TOT_PROCS} ENVIRONMENT OMP_NUM_THREADS=${THREADS}
                                                PROCESSOR_AFFINITY TRUE)

    if(ENABLE_CUDA
       OR ENABLE_ROCM
       OR ENABLE_SYCL
       OR ENABLE_OFFLOAD)
      set_tests_properties(${TESTNAME} PROPERTIES RESOURCE_LOCK exclusively_owned_gpus)
    endif()

    if(ENABLE_OFFLOAD)
      set_property(
        TEST ${TESTNAME}
        APPEND
        PROPERTY ENVIRONMENT "OMP_TARGET_OFFLOAD=mandatory")
    endif()
  endif()

  set(TEST_LABELS_TEMP "")
  add_test_labels(${TESTNAME} TEST_LABELS_TEMP)
  set_property(TEST ${TESTNAME} APPEND PROPERTY LABELS "unit")
endfunction()

# Add a test to see if the target output exists in the desired location in the build directory.
function(add_test_target_in_output_location TARGET_NAME_TO_TEST EXE_DIR_RELATIVE_TO_BUILD)

  # obtain BASE_NAME
  get_target_property(BASE_NAME ${TARGET_NAME_TO_TEST} OUTPUT_NAME)
  if(NOT BASE_NAME)
    set(BASE_NAME ${TARGET_NAME_TO_TEST})
  endif()

  set(TESTNAME build_output_${TARGET_NAME_TO_TEST}_exists)
  add_test(NAME ${TESTNAME} COMMAND ls ${qmcpack_BINARY_DIR}/bin/${BASE_NAME})

  set_property(TEST ${TESTNAME} APPEND PROPERTY LABELS "unit;deterministic")
endfunction()
