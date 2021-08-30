include("${PROJECT_SOURCE_DIR}/CMake/test_labels.cmake")

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

    if(QMC_CUDA
       OR ENABLE_CUDA
       OR ENABLE_ROCM
       OR ENABLE_OFFLOAD)
      set_tests_properties(${TESTNAME} PROPERTIES RESOURCE_LOCK exclusively_owned_gpus)
    endif()
  endif()

  set(TEST_LABELS_TEMP "")
  add_test_labels(${TESTNAME} TEST_LABELS_TEMP)
  set_property(
    TEST ${TESTNAME}
    APPEND
    PROPERTY LABELS "unit")
endfunction()
