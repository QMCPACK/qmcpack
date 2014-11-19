#############################################################
# Add tests to ctest
#############################################################
# Useful macros to compile and run an executable:
#   ADD_QMC_PROVISIONAL_TEST( EXECFILE )
#       Add a provisional test that will be compiled 
#       but not executed
#   ADD_QMC_TEST( EXECFILE ARGS )
#       Add a serial test passing the given args to the test
#   ADD_QMC_TEST_PARALLEL( EXECFILE PROCS ARGS )
#       Add a parallel test with the given number of
#       processors and arguments
#   ADD_QMC_TEST_1_2_4( EXECFILE args )
#       Add a test that will run on 1, 2, and 4 processors
#   ADD_QMC_WEEKLY_TEST( EXECFILE PROCS ARGS )
#       Add a "WEEKLY" parallel test with the given number 
#       of processors and arguments
#   ADD_QMC_TEST_THREAD_MPI( EXECFILE PROCS THREADS ARGS )
#       Add a parallel test with the given number of threads
#       per processes
# Useful macros to run an existing executable:
#   RUN_QMC_APP( TESTNAME SRC_DIR PROCS THREADS ARGS )
#############################################################

message("Attempting to add tests for QMCAPP")

set(QMCAPP_ERR "")
if((NOT ${TEST_MAX_PROCS}) OR (${TEST_MAX_PROCS} STREQUAL ""))
  set(QMCAPP_ERR "${QMCAPP_ERR}  TEST_MAX_PROCS not set: '${TEST_MAX_PROCS}'\n")
endif()


if(NOT ${QMCAPP_ERR} STREQUAL "")
  message("${QMCAPP_ERR}  skipping tests")
else()
  INCLUDE( "${qmcpack_SOURCE_DIR}/CMake/macros.cmake" )
  IF ( NOT QMC_DATA ) 
      MESSAGE( WARNING "QMC_DATA is not set, some tests will fail")
  ENDIF()
  
  # Add a simple hello world test
  ADD_QMC_TEST_1_2_4( HelloWorld )
  
  # Run the deterministic tests (sucess is running without errors)
  COPY_DIRECTORY( "${QMC_DATA}/dft" "${CMAKE_CURRENT_BINARY_DIR}/dft" )
  RUN_QMC_APP( deterministic-1-1 "${QMC_DATA}/deterministic" 1 1 qmc.in.xml )
  RUN_QMC_APP( deterministic-4-1 "${QMC_DATA}/deterministic" 4 1 qmc.in.xml )
  RUN_QMC_APP( deterministic-1-4 "${QMC_DATA}/deterministic" 1 4 qmc.in.xml )
  
  # Run the statistical tests
  RUN_QMC_APP( statistical_vmc-1-4 "${QMC_DATA}/statistical_vmc" 4 1 qmc.in.xml )
  RUN_QMC_APP( statistical_vmc-4-4 "${QMC_DATA}/statistical_vmc" 1 4 qmc.in.xml )
endif()
