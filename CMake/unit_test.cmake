
INCLUDE("${PROJECT_SOURCE_DIR}/CMake/test_labels.cmake")

# Runs unit tests
#
# The NO_MPI_LINK option is to allow different MPIEXEC_PREFLAGS to be passed in the rare case
# such as on summit that special flags are required for serial tests (in the mpi sense)
# to run under the mpirunner.
#
# On summit this flag is '--smpiargs="-disable_gpu_hooks"'
#

FUNCTION( ADD_UNIT_TEST TESTNAME PROCS THREADS TEST_BINARY )
  set(options NO_MPI_LINK)
  set(oneValueArgs)
  set(multValueArgs)
  cmake_parse_arguments(ADD_UNIT_TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    message_verbose("Adding test ${TESTNAME}")
    math( EXPR TOT_PROCS "${PROCS} * ${THREADS}" )
    if ( HAVE_MPI )
      if (ADD_UNIT_TEST_NO_MPI_LINK)
        add_test(NAME ${TESTNAME} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_NO_MPI_PREFLAGS} ${TEST_BINARY} ${ADD_UNIT_TEST_UNPARSED_ARGUMENTS})
      else()
        add_test(NAME ${TESTNAME} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS} ${TEST_BINARY} ${ADD_UNIT_TEST_UNPARSED_ARGUMENTS})
      endif()
        set( TEST_ADDED TRUE )
    else()
        if ( ( ${PROCS} STREQUAL "1" ) )
            add_test(NAME ${TESTNAME} COMMAND ${TEST_BINARY} ${ADD_UNIT_TEST_UNPARSED_ARGUMENTS})
            set( TEST_ADDED TRUE )
        else()
            message_verbose("Disabling test ${TESTNAME} (building without MPI)")
        endif()
    endif()

    if (TEST_ADDED)
        set_tests_properties( ${TESTNAME} PROPERTIES
                              PROCESSORS ${TOT_PROCS}
                              ENVIRONMENT OMP_NUM_THREADS=${THREADS}
                              PROCESSOR_AFFINITY TRUE )

        if (QMC_CUDA OR ENABLE_CUDA OR ENABLE_OFFLOAD)
            set_tests_properties(${TESTNAME} PROPERTIES RESOURCE_LOCK exclusively_owned_gpus)
        endif()
    endif()

    set(TEST_LABELS_TEMP "")
    add_test_labels( ${TESTNAME} TEST_LABELS_TEMP )
    set_property(TEST ${TESTNAME} APPEND PROPERTY LABELS "unit")
ENDFUNCTION()

