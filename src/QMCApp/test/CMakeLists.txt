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

#   COPY_DIRECTORY("${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae" "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_dimer_ae") 
#    SET (USE_EXT_MPI_FOR_SERIAL_TESTS)

IF(NOT QMC_CUDA)

IF (NOT QMC_COMPLEX)
#
# LiH molecular dimer gaussian tests, all electron
# Also check results for different number of mpi tasks and threads keeping total constant
#
    RUN_QMC_APP(short-LiH_dimer_ae-vmc_hf_noj-16-1 "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae" 16 1 vmc_hf_noj.in.xml) 
    ADD_TEST( NAME short-LiH_dimer_ae-vmc_hf_noj-16-1-kinetic
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --ke "7.9914014383 0.03"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_dimer_ae-vmc_hf_noj-16-1"
    )
    SET_PROPERTY( TEST short-LiH_dimer_ae-vmc_hf_noj-16-1-kinetic APPEND PROPERTY DEPENDS short-LiH_dimer_ae-vmc_hf_noj-16-1 )
   ADD_TEST( NAME short-LiH_dimer_ae-vmc_hf_noj-16-1-totenergy
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --le "-7.9873236457 0.003"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_dimer_ae-vmc_hf_noj-16-1"
    )
    SET_PROPERTY( TEST short-LiH_dimer_ae-vmc_hf_noj-16-1-totenergy APPEND PROPERTY DEPENDS short-LiH_dimer_ae-vmc_hf_noj-16-1 )
   ADD_TEST( NAME short-LiH_dimer_ae-vmc_hf_noj-16-1-eeenergy
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --ee "3.4888166386 0.0024"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_dimer_ae-vmc_hf_noj-16-1"
    )
    SET_PROPERTY( TEST short-LiH_dimer_ae-vmc_hf_noj-16-1-totenergy APPEND PROPERTY DEPENDS short-LiH_dimer_ae-vmc_hf_noj-16-1 )
  ADD_TEST( NAME short-LiH_dimer_ae-vmc_hf_noj-16-1-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p vmc_hf_noj --ts "1600000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_dimer_ae-vmc_hf_noj-16-1"
     )
  SET_PROPERTY( TEST short-LiH_dimer_ae-vmc_hf_noj-16-1-samples APPEND PROPERTY DEPENDS short-LiH_dimer_ae-vmc_hf_noj-16-1 )


    RUN_QMC_APP(short-LiH_dimer_ae-vmc_hf_noj-4-4 "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae" 4 4 vmc_hf_noj.in.xml) 
    ADD_TEST( NAME short-LiH_dimer_ae-vmc_hf_noj-4-4-totenergy
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --le "-7.9873236457 0.003"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_dimer_ae-vmc_hf_noj-4-4"
    )
    SET_PROPERTY( TEST short-LiH_dimer_ae-vmc_hf_noj-4-4-totenergy APPEND PROPERTY DEPENDS short-LiH_dimer_ae-vmc_hf_noj-4-4 )
  ADD_TEST( NAME short-LiH_dimer_ae-vmc_hf_noj-4-4-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p vmc_hf_noj --ts "1600000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_dimer_ae-vmc_hf_noj-4-4"
     )
  SET_PROPERTY( TEST short-LiH_dimer_ae-vmc_hf_noj-4-4-samples APPEND PROPERTY DEPENDS short-LiH_dimer_ae-vmc_hf_noj-4-4 )
    
    RUN_QMC_APP(short-LiH_dimer_ae-vmc_hf_noj-1-16 "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae" 1 16 vmc_hf_noj.in.xml) 
    ADD_TEST( NAME short-LiH_dimer_ae-vmc_hf_noj-1-16-totenergy
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --le "-7.9873236457 0.003"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_dimer_ae-vmc_hf_noj-1-16"
    )
    SET_PROPERTY( TEST short-LiH_dimer_ae-vmc_hf_noj-1-16-totenergy APPEND PROPERTY DEPENDS short-LiH_dimer_ae-vmc_hf_noj-1-16 )
  ADD_TEST( NAME short-LiH_dimer_ae-vmc_hf_noj-1-16-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p vmc_hf_noj --ts "1600000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_dimer_ae-vmc_hf_noj-1-16"
     )
  SET_PROPERTY( TEST short-LiH_dimer_ae-vmc_hf_noj-1-16-samples APPEND PROPERTY DEPENDS short-LiH_dimer_ae-vmc_hf_noj-1-16 )

   RUN_QMC_APP(short-LiH_dimer_ae-vmc_hf_noj_nodrift-16-1 "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae" 16 1 vmc_hf_noj_nodrift.in.xml) 
    ADD_TEST( NAME short-LiH_dimer_ae-vmc_hf_noj_nodrift-16-1-totenergy
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj_nodrift -e 2 --le "-7.9873236457 0.003"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_dimer_ae-vmc_hf_noj_nodrift-16-1"
    )
    SET_PROPERTY( TEST short-LiH_dimer_ae-vmc_hf_noj_nodrift-16-1-totenergy APPEND PROPERTY DEPENDS short-LiH_dimer_ae-vmc_hf_noj_nodrift-16-1 )
  ADD_TEST( NAME short-LiH_dimer_ae-vmc_hf_noj_nodrift-16-1-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p vmc_hf_noj_nodrift --ts "1600000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_dimer_ae-vmc_hf_noj_nodrift-16-1"
     )
  SET_PROPERTY( TEST short-LiH_dimer_ae-vmc_hf_noj_nodrift-16-1-samples APPEND PROPERTY DEPENDS short-LiH_dimer_ae-vmc_hf_noj_nodrift-16-1 )

#
# LiH molecular dimer pseudopotential tests
#
   RUN_QMC_APP(short-LiH_dimer_pp-vmc_hf_noj-1-16 "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp" 1 16 vmc_hf_noj.in.xml) 
    ADD_TEST( NAME short-LiH_dimer_pp-vmc_hf_noj-1-16-totenergy
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --le "-0.7507264125 0.00043"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_dimer_pp-vmc_hf_noj-1-16"
    )
    SET_PROPERTY( TEST short-LiH_dimer_pp-vmc_hf_noj-1-16-totenergy APPEND PROPERTY DEPENDS short-LiH_dimer_pp-vmc_hf_noj-1-16 )
    ADD_TEST( NAME short-LiH_dimer_pp-vmc_hf_noj-1-16-kinetic
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --ke "0.6354491487 0.0013"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_dimer_pp-vmc_hf_noj-1-16"
    )
    SET_PROPERTY( TEST short-LiH_dimer_pp-vmc_hf_noj-1-16-kinetic APPEND PROPERTY DEPENDS short-LiH_dimer_pp-vmc_hf_noj-1-16 )
    ADD_TEST( NAME short-LiH_dimer_pp-vmc_hf_noj-1-16-eeenergy
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --ee "0.4814819907 0.00061"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_dimer_pp-vmc_hf_noj-1-16"
    )
    SET_PROPERTY( TEST short-LiH_dimer_pp-vmc_hf_noj-1-16-eeenergy APPEND PROPERTY DEPENDS short-LiH_dimer_pp-vmc_hf_noj-1-16 )
    ADD_TEST( NAME short-LiH_dimer_pp-vmc_hf_noj-1-16-potential
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p vmc_hf_noj -e 2 --lp "-1.3861755612 0.0013"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_dimer_pp-vmc_hf_noj-1-16"
    )
    SET_PROPERTY( TEST short-LiH_dimer_pp-vmc_hf_noj-1-16-potential APPEND PROPERTY DEPENDS short-LiH_dimer_pp-vmc_hf_noj-1-16 )
  ADD_TEST( NAME short-LiH_dimer_pp-vmc_hf_noj-1-16-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p vmc_hf_noj --ts "1600000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_dimer_pp-vmc_hf_noj-1-16"
     )
  SET_PROPERTY( TEST short-LiH_dimer_pp-vmc_hf_noj-1-16-samples APPEND PROPERTY DEPENDS short-LiH_dimer_pp-vmc_hf_noj-1-16 )
 ELSE()
    MESSAGE("Skipping all gaussian basis tests because they are not supported by complex build (QMC_COMPLEX=1)")  
ENDIF()
ELSE()
    MESSAGE("Skipping all gaussian basis tests because they are not supported by CUDA build (QMC_CUDA=1)")  
ENDIF()

#
# LiH solid einspline wavefunction tests for different k-points
# Wavefunctions from quantum espresso. Reference energies obtained via pw2casino option
#
    RUN_QMC_APP(short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp" 1 16 hf_vmc_LiH-gamma.xml) 
    ADD_TEST( NAME short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16-kinetic
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p hf_vmc_LiH-gamma -e 2 --ke "7.4145221704067854 0.024"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16"
    )
    SET_PROPERTY( TEST short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16-kinetic APPEND PROPERTY DEPENDS short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 )
    ADD_TEST( NAME short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16-localpotential
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p hf_vmc_LiH-gamma -e 2 --lpp "-11.718808431876706 0.025"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16"
    )
    SET_PROPERTY( TEST short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16-localpotential APPEND PROPERTY DEPENDS short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 )
  ADD_TEST( NAME short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p hf_vmc_LiH-gamma --ts "1600000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16"
     )
  SET_PROPERTY( TEST short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16-samples APPEND PROPERTY DEPENDS short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj-1-16 )

   RUN_QMC_APP(short-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp" 1 16 hf_vmc_LiH-x.xml) 
    ADD_TEST( NAME short-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16-kinetic
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p hf_vmc_LiH-x -e 2 --ke "7.3772282047321536 0.024"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16"
    )
    SET_PROPERTY( TEST short-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16-kinetic APPEND PROPERTY DEPENDS short-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 )
    ADD_TEST( NAME short-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16-localpotential
        COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p hf_vmc_LiH-x -e 2 --lpp "-11.414410381100344 0.025"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16"
    )
    SET_PROPERTY( TEST short-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16-localpotential APPEND PROPERTY DEPENDS short-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 )
  ADD_TEST( NAME short-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p hf_vmc_LiH-x --ts "1600000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16"
     )
  SET_PROPERTY( TEST short-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16-samples APPEND PROPERTY DEPENDS short-LiH_solid_1x1x1_pp-x-vmc_hf_noj-1-16 )

# Arbitrary k-point requires complex build to run. Ensure fails with real/non-complex build
    RUN_QMC_APP(short-LiH_solid_1x1x1_pp-arb-vmc_hf_noj-1-16 "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp" 1 16 hf_vmc_LiH-arb.xml) 
    IF (NOT QMC_COMPLEX)
       #message ("Not a complex build")
       SET_PROPERTY(TEST short-LiH_solid_1x1x1_pp-arb-vmc_hf_noj-1-16 APPEND PROPERTY WILL_FAIL TRUE)
    ELSE()
       #message ("Is a complex build")
       ADD_TEST( NAME short-LiH_solid_1x1x1_pp-arb-vmc_hf_noj-1-16-kinetic
           COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p hf_vmc_LiH-arb -e 2 --ke "7.4201207734399359 0.024"
           WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_solid_1x1x1_pp-arb-vmc_hf_noj-1-16"
       )
       SET_PROPERTY( TEST short-LiH_solid_1x1x1_pp-arb-vmc_hf_noj-1-16-kinetic APPEND PROPERTY DEPENDS short-LiH_solid_1x1x1_pp-arb-vmc_hf_noj-1-16 )
       ADD_TEST( NAME short-LiH_solid_1x1x1_pp-arb-vmc_hf_noj-1-16-localpotential
           COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p hf_vmc_LiH-arb -e 2 --lpp "-11.615276430529056 0.025"
           WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_solid_1x1x1_pp-arb-vmc_hf_noj-1-16"
       )
       SET_PROPERTY( TEST short-LiH_solid_1x1x1_pp-arb-vmc_hf_noj-1-16-localpotential APPEND PROPERTY DEPENDS short-LiH_solid_1x1x1_pp-arb-vmc_hf_noj-1-16 )
     ADD_TEST( NAME short-LiH_solid_1x1x1_pp-arb-vmc_hf_noj-1-16-samples
            COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p hf_vmc_LiH-arb --ts "1600000 0.0"
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-LiH_solid_1x1x1_pp-arb-vmc_hf_noj-1-16"
        )
     SET_PROPERTY( TEST short-LiH_solid_1x1x1_pp-arb-vmc_hf_noj-1-16-samples APPEND PROPERTY DEPENDS short-LiH_solid_1x1x1_pp-arb-vmc_hf_noj-1-16 )
    ENDIF() 

#
# Short tests, about 2-3 seconds on 16 cores
#
  RUN_QMC_APP(short-bccH_1x1x1_ae-vmc_sdj-1-16 "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae" 1 16 qmc_short.in.xml) 
  ADD_TEST( NAME short-bccH_1x1x1_ae-vmc_sdj-1-16-totenergy 
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p qmc_short -e 2 --le "-1.834032 0.0003"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-bccH_1x1x1_ae-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST short-bccH_1x1x1_ae-vmc_sdj-1-16-totenergy APPEND PROPERTY DEPENDS short-bccH_1x1x1_ae-vmc_sdj-1-16 )
  ADD_TEST( NAME short-bccH_1x1x1_ae-vmc_sdj-1-16-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p qmc_short --ts "960000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-bccH_1x1x1_ae-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST short-bccH_1x1x1_ae-vmc_sdj-1-16-samples APPEND PROPERTY DEPENDS short-bccH_1x1x1_ae-vmc_sdj-1-16 )


  RUN_QMC_APP(short-diamondC_1x1x1_pp-vmc_sdj-1-16 "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp" 1 16 qmc_short.in.xml) 
  ADD_TEST( NAME short-diamondC_1x1x1_pp-vmc_sdj-1-16-totenergy 
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p qmc_short -e 2 --le "-10.491445 0.002"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-diamondC_1x1x1_pp-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST short-diamondC_1x1x1_pp-vmc_sdj-1-16-totenergy APPEND PROPERTY DEPENDS short-diamondC_1x1x1_pp-vmc_sdj-1-16 )
  ADD_TEST( NAME short-diamondC_1x1x1_pp-vmc_sdj-1-16-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p qmc_short --ts "128000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-diamondC_1x1x1_pp-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST short-diamondC_1x1x1_pp-vmc_sdj-1-16-samples APPEND PROPERTY DEPENDS short-diamondC_1x1x1_pp-vmc_sdj-1-16 )
  
  RUN_QMC_APP(short-diamondC_2x1x1_pp-vmc_sdj-1-16 "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp" 1 16 qmc_short.in.xml) 
  ADD_TEST( NAME short-diamondC_2x1x1_pp-vmc_sdj-1-16-totenergy 
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p qmc_short -e 2 --le "-21.668293 0.007"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-diamondC_2x1x1_pp-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST short-diamondC_2x1x1_pp-vmc_sdj-1-16-totenergy APPEND PROPERTY DEPENDS short-diamondC_2x1x1_pp-vmc_sdj-1-16 )
  ADD_TEST( NAME short-diamondC_2x1x1_pp-vmc_sdj-1-16-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p qmc_short --ts "32000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-diamondC_2x1x1_pp-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST short-diamondC_2x1x1_pp-vmc_sdj-1-16-samples APPEND PROPERTY DEPENDS short-diamondC_2x1x1_pp-vmc_sdj-1-16 )
  
  RUN_QMC_APP(short-hcpBe_1x1x1_pp-vmc_sdj-1-16 "${CMAKE_SOURCE_DIR}/tests/solids/hcpBe_1x1x1_pp" 1 16 qmc_short.in.xml) 
  ADD_TEST( NAME short-hcpBe_1x1x1_pp-vmc_sdj-1-16-totenergy 
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p qmc_short -e 2 --le "-1.481656 0.0005"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-hcpBe_1x1x1_pp-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST short-hcpBe_1x1x1_pp-vmc_sdj-1-16-totenergy APPEND PROPERTY DEPENDS short-hcpBe_1x1x1_pp-vmc_sdj-1-16 )
  ADD_TEST( NAME short-hcpBe_1x1x1_pp-vmc_sdj-1-16-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p qmc_short --ts "416000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-hcpBe_1x1x1_pp-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST short-hcpBe_1x1x1_pp-vmc_sdj-1-16-samples APPEND PROPERTY DEPENDS short-hcpBe_1x1x1_pp-vmc_sdj-1-16 )
  
  
  RUN_QMC_APP(short-monoO_1x1x1_pp-vmc_sdj-1-16 "${CMAKE_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp" 1 16 qmc_short.in.xml)
  ADD_TEST( NAME short-monoO_1x1x1_pp-vmc_sdj-1-16-totenergy 
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p qmc_short -e 2 --le "-31.776596 0.006"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-monoO_1x1x1_pp-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST short-monoO_1x1x1_pp-vmc_sdj-1-16-totenergy APPEND PROPERTY DEPENDS short-monoO_1x1x1_pp-vmc_sdj-1-16 )
  ADD_TEST( NAME short-monoO_1x1x1_pp-vmc_sdj-1-16-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p qmc_short --ts "64000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/short-monoO_1x1x1_pp-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST short-monoO_1x1x1_pp-vmc_sdj-1-16-samples APPEND PROPERTY DEPENDS short-monoO_1x1x1_pp-vmc_sdj-1-16 )
  


#
# Long variants of above, about 1/2 hour with 16 cores
#
  SET( CTEST_TEST_TIMEOUT 36000 )  # 10h Should make a macro to avoid changing global settings
  RUN_QMC_APP(long-bccH_1x1x1_ae-vmc_sdj-1-16 "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae" 1 16 qmc_long.in.xml) 
  ADD_TEST( NAME long-bccH_1x1x1_ae-vmc_sdj-1-16-totenergy 
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p qmc_long -e 2 --le "-1.834032 0.00001"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/long-bccH_1x1x1_ae-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST long-bccH_1x1x1_ae-vmc_sdj-1-16-totenergy APPEND PROPERTY DEPENDS long-bccH_1x1x1_ae-vmc_sdj-1-16 )
  ADD_TEST( NAME long-bccH_1x1x1_ae-vmc_sdj-1-16-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p qmc_long --ts "921600000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/long-bccH_1x1x1_ae-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST long-bccH_1x1x1_ae-vmc_sdj-1-16-samples APPEND PROPERTY DEPENDS long-bccH_1x1x1_ae-vmc_sdj-1-16 )
  
  RUN_QMC_APP(long-diamondC_1x1x1_pp-vmc_sdj-1-16 "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp" 1 16 qmc_long.in.xml) 
  ADD_TEST( NAME long-diamondC_1x1x1_pp-vmc_sdj-1-16-totenergy 
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p qmc_long -e 2 --le "-10.491445 0.000065"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/long-diamondC_1x1x1_pp-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST long-diamondC_1x1x1_pp-vmc_sdj-1-16-totenergy APPEND PROPERTY DEPENDS long-diamondC_1x1x1_pp-vmc_sdj-1-16 )
  ADD_TEST( NAME long-diamondC_1x1x1_pp-vmc_sdj-1-16-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p qmc_long --ts "122880000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/long-diamondC_1x1x1_pp-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST long-diamondC_1x1x1_pp-vmc_sdj-1-16-samples APPEND PROPERTY DEPENDS long-diamondC_1x1x1_pp-vmc_sdj-1-16 )
  
  RUN_QMC_APP(long-diamondC_2x1x1_pp-vmc_sdj-1-16 "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp" 1 16 qmc_long.in.xml) 
  ADD_TEST( NAME long-diamondC_2x1x1_pp-vmc_sdj-1-16-totenergy 
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p qmc_long -e 2 --le "-21.668293 0.00023"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/long-diamondC_2x1x1_pp-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST long-diamondC_2x1x1_pp-vmc_sdj-1-16-totenergy APPEND PROPERTY DEPENDS long-diamondC_2x1x1_pp-vmc_sdj-1-16 )
  ADD_TEST( NAME long-diamondC_2x1x1_pp-vmc_sdj-1-16-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p qmc_long --ts "30720000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/long-diamondC_2x1x1_pp-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST long-diamondC_2x1x1_pp-vmc_sdj-1-16-samples APPEND PROPERTY DEPENDS long-diamondC_2x1x1_pp-vmc_sdj-1-16 )
  
  RUN_QMC_APP(long-hcpBe_1x1x1_pp-vmc_sdj-1-16 "${CMAKE_SOURCE_DIR}/tests/solids/hcpBe_1x1x1_pp" 1 16 qmc_long.in.xml) 
  ADD_TEST( NAME long-hcpBe_1x1x1_pp-vmc_sdj-1-16-totenergy 
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p qmc_long -e 2 --le "-1.481656 0.000017"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/long-hcpBe_1x1x1_pp-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST long-hcpBe_1x1x1_pp-vmc_sdj-1-16-totenergy APPEND PROPERTY DEPENDS long-hcpBe_1x1x1_pp-vmc_sdj-1-16 )
  ADD_TEST( NAME long-hcpBe_1x1x1_pp-vmc_sdj-1-16-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p qmc_long --ts "399360000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/long-hcpBe_1x1x1_pp-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST long-hcpBe_1x1x1_pp-vmc_sdj-1-16-samples APPEND PROPERTY DEPENDS long-hcpBe_1x1x1_pp-vmc_sdj-1-16 )
  
  RUN_QMC_APP(long-monoO_1x1x1_pp-vmc_sdj-1-16 "${CMAKE_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp" 1 16 qmc_long.in.xml)
  ADD_TEST( NAME long-monoO_1x1x1_pp-vmc_sdj-1-16-totenergy 
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p qmc_long -e 2 --le "-31.776596 0.0002"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/long-monoO_1x1x1_pp-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST long-monoO_1x1x1_pp-vmc_sdj-1-16-totenergy APPEND PROPERTY DEPENDS long-monoO_1x1x1_pp-vmc_sdj-1-16 )
  ADD_TEST( NAME long-monoO_1x1x1_pp-vmc_sdj-1-16-samples
         COMMAND ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 0 -p qmc_long --ts "61440000 0.0"
         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/long-monoO_1x1x1_pp-vmc_sdj-1-16"
     )
  SET_PROPERTY( TEST long-monoO_1x1x1_pp-vmc_sdj-1-16-samples APPEND PROPERTY DEPENDS long-monoO_1x1x1_pp-vmc_sdj-1-16 )

endif()
