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

message("Adding tests for QMCAPP")

set(QMCAPP_ERR "")
if((NOT ${TEST_MAX_PROCS}) OR (${TEST_MAX_PROCS} STREQUAL ""))
    set(QMCAPP_ERR "${QMCAPP_ERR}  TEST_MAX_PROCS not set: '${TEST_MAX_PROCS}'\n")
endif()


if(NOT ${QMCAPP_ERR} STREQUAL "")
    message("${QMCAPP_ERR}  skipping tests")
else()
    INCLUDE( "${qmcpack_SOURCE_DIR}/CMake/macros.cmake" )
#    IF ( NOT QMC_DATA ) 
#        MESSAGE( WARNING "QMC_DATA is not set, some tests will fail")
#    ENDIF()
#
#    # Add a simple hello world test
#    ADD_QMC_TEST_1_2_4( HelloWorld )
#
#    # Run the deterministic tests (sucess is running without errors)
#    COPY_DIRECTORY( "${QMC_DATA}/dft" "${CMAKE_CURRENT_BINARY_DIR}/dft" )
#    RUN_QMC_APP( deterministic-1-1 "${QMC_DATA}/deterministic" 1 1 qmc.in.xml )
#    RUN_QMC_APP( deterministic-4-1 "${QMC_DATA}/deterministic" 4 1 qmc.in.xml )
#    RUN_QMC_APP( deterministic-1-4 "${QMC_DATA}/deterministic" 1 4 qmc.in.xml )
#
#    # Run the statistical tests
#    RUN_QMC_APP( statistical_vmc-1-4 "${QMC_DATA}/statistical_vmc" 4 1 qmc.in.xml )
#    RUN_QMC_APP( statistical_vmc-4-4 "${QMC_DATA}/statistical_vmc" 1 4 qmc.in.xml )
#
#    # Add the example dmc test:
#    #   mpirun -np 16 qmcapp dmc.xml >& test.out
#    #   check_scalars.py  --ns 2  -p Li  -s '1 3' -e 10  --le '-7.478011 0.000035 -7.478059 0.000035'
#    RUN_QMC_APP( dmc-16 "${QMC_DATA}/dmc" 16 1 dmc.xml )
#    ADD_TEST( NAME dmc-16-check 
#        COMMAND check_scalars.py --ns 2 -p Li -s "1 3" -e 10 --le "-7.478011 0.000035 -7.478059 0.000035"
#        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/dmc-16"
#    )
#    SET_PROPERTY( TEST dmc-16-check APPEND PROPERTY DEPENDS dmc-16 )

#   COPY_DIRECTORY("${CMAKE_SOURCE_DIR}/examples/molecules/tests/LiH_dimer_ae" "${CMAKE_CURRENT_BINARY_DIR}/LiH_dimer_ae") 
#    SET (USE_EXT_MPI_FOR_SERIAL_TESTS)

    RUN_QMC_APP(LiH_dimer_ae-vmc_hf_noj-16-1 "${CMAKE_SOURCE_DIR}/examples/molecules/tests/LiH_dimer_ae" 16 1 vmc_hf_noj.in.xml) 
    ADD_TEST( NAME LiH_dimer_ae-vmc_hf_noj-16-1-kinetic
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --ke "7.9914014383 0.03"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/LiH_dimer_ae-vmc_hf_noj-16-1"
    )
    SET_PROPERTY( TEST LiH_dimer_ae-vmc_hf_noj-16-1-kinetic APPEND PROPERTY DEPENDS LiH_dimer_ae-vmc_hf_noj-16-1 )
   ADD_TEST( NAME LiH_dimer_ae-vmc_hf_noj-16-1-totenergy
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --le "-7.9873236457 0.003"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/LiH_dimer_ae-vmc_hf_noj-16-1"
    )
    SET_PROPERTY( TEST LiH_dimer_ae-vmc_hf_noj-16-1-totenergy APPEND PROPERTY DEPENDS LiH_dimer_ae-vmc_hf_noj-16-1 )
   ADD_TEST( NAME LiH_dimer_ae-vmc_hf_noj-16-1-eeenergy
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --ee "3.4888166386 0.0024"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/LiH_dimer_ae-vmc_hf_noj-16-1"
    )
    SET_PROPERTY( TEST LiH_dimer_ae-vmc_hf_noj-16-1-totenergy APPEND PROPERTY DEPENDS LiH_dimer_ae-vmc_hf_noj-16-1 )
    
    RUN_QMC_APP(LiH_dimer_ae-vmc_hf_noj-4-4 "${CMAKE_SOURCE_DIR}/examples/molecules/tests/LiH_dimer_ae" 4 4 vmc_hf_noj.in.xml) 
    ADD_TEST( NAME LiH_dimer_ae-vmc_hf_noj-4-4-totenergy
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --le "-7.9873236457 0.003"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/LiH_dimer_ae-vmc_hf_noj-4-4"
    )
    SET_PROPERTY( TEST LiH_dimer_ae-vmc_hf_noj-4-4-totenergy APPEND PROPERTY DEPENDS LiH_dimer_ae-vmc_hf_noj-4-4 )
    
    RUN_QMC_APP(LiH_dimer_ae-vmc_hf_noj-1-16 "${CMAKE_SOURCE_DIR}/examples/molecules/tests/LiH_dimer_ae" 1 16 vmc_hf_noj.in.xml) 
    ADD_TEST( NAME LiH_dimer_ae-vmc_hf_noj-1-16-totenergy
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --le "-7.9873236457 0.003"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/LiH_dimer_ae-vmc_hf_noj-1-16"
    )
    SET_PROPERTY( TEST LiH_dimer_ae-vmc_hf_noj-1-16-totenergy APPEND PROPERTY DEPENDS LiH_dimer_ae-vmc_hf_noj-1-16 )

   RUN_QMC_APP(LiH_dimer_ae-vmc_hf_noj_nodrift-16-1 "${CMAKE_SOURCE_DIR}/examples/molecules/tests/LiH_dimer_ae" 16 1 vmc_hf_noj_nodrift.in.xml) 
    ADD_TEST( NAME LiH_dimer_ae-vmc_hf_noj_nodrift-16-1-totenergy
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj_nodrift -e 2 --le "-7.9873236457 0.003"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/LiH_dimer_ae-vmc_hf_noj_nodrift-16-1"
    )
    SET_PROPERTY( TEST LiH_dimer_ae-vmc_hf_noj_nodrift-16-1-totenergy APPEND PROPERTY DEPENDS LiH_dimer_ae-vmc_hf_noj_nodrift-16-1 )

   RUN_QMC_APP(LiH_dimer_pp-vmc_hf_noj-1-16 "${CMAKE_SOURCE_DIR}/examples/molecules/tests/LiH_dimer_pp" 1 16 vmc_hf_noj.in.xml) 
    ADD_TEST( NAME LiH_dimer_pp-vmc_hf_noj-1-16-totenergy
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --le "-0.7507264125 0.00043"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/LiH_dimer_pp-vmc_hf_noj-1-16"
    )
    SET_PROPERTY( TEST LiH_dimer_pp-vmc_hf_noj-1-16-totenergy APPEND PROPERTY DEPENDS LiH_dimer_pp-vmc_hf_noj-1-16 )
    ADD_TEST( NAME LiH_dimer_pp-vmc_hf_noj-1-16-kinetic
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --ke "0.6354491487 0.0013"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/LiH_dimer_pp-vmc_hf_noj-1-16"
    )
    SET_PROPERTY( TEST LiH_dimer_pp-vmc_hf_noj-1-16-kinetic APPEND PROPERTY DEPENDS LiH_dimer_pp-vmc_hf_noj-1-16 )
    ADD_TEST( NAME LiH_dimer_pp-vmc_hf_noj-1-16-eeenergy
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --ee "0.4814819907 0.00061"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/LiH_dimer_pp-vmc_hf_noj-1-16"
    )
    SET_PROPERTY( TEST LiH_dimer_pp-vmc_hf_noj-1-16-eeenergy APPEND PROPERTY DEPENDS LiH_dimer_pp-vmc_hf_noj-1-16 )
    ADD_TEST( NAME LiH_dimer_pp-vmc_hf_noj-1-16-potential
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --lp "-1.3861755612 0.0013"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/LiH_dimer_pp-vmc_hf_noj-1-16"
    )
    SET_PROPERTY( TEST LiH_dimer_pp-vmc_hf_noj-1-16-potential APPEND PROPERTY DEPENDS LiH_dimer_pp-vmc_hf_noj-1-16 )
    
endif()


