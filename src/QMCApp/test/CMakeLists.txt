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


# Add a test run and associated scalar checks
# BASE_NAME - name of test (number of MPI processes, number of threads, and value to check (if applicable)
#             will be appended to get the full test name)
# BASE_DIR - source location of test input files
# PREFIX - prefix for output files
# INPUT_FILE - XML input file to QMCPACK
# PROCS - number of MPI processes
# THREADS - number of OpenMP threads
# SCALAR_VALUES - list of output values to check with check_scalars.py
#                 The list entries alternate between the value name and the value (usually a string with the
#                 both the average and error).
# SHOULD_SUCCEED - whether the test is expected to pass or fail.  Expected failing tests will not have
#                  the scalar tests added.

FUNCTION(QMC_RUN_AND_CHECK BASE_NAME BASE_DIR PREFIX INPUT_FILE PROCS THREADS SCALAR_VALUES SHOULD_SUCCEED)
    # Map from name of check to appropriate flag for check_scalars.py
    LIST(APPEND SCALAR_CHECK_TYPE "kinetic" "totenergy" "eeenergy" "samples" "potential" "localpotential" "flux")
    LIST(APPEND CHECK_SCALAR_FLAG "--ke"    "--le"      "--ee"     "--ts"    "--lp"      "--lpp"          "--fl")

    SET( TEST_ADDED FALSE )
    SET( FULL_NAME "${BASE_NAME}-${PROCS}-${THREADS}" )
    MESSAGE("Adding test ${FULL_NAME}")
    RUN_QMC_APP(${FULL_NAME} ${BASE_DIR} ${PROCS} ${THREADS} TEST_ADDED ${INPUT_FILE})


    IF ( TEST_ADDED AND NOT SHOULD_SUCCEED)
        SET_PROPERTY(TEST ${FULL_NAME} APPEND PROPERTY WILL_FAIL TRUE)
        #MESSAGE("Test ${FULL_NAME} should fail")
    ENDIF()

    IF ( TEST_ADDED AND SHOULD_SUCCEED)
        FOREACH(SCALAR_CHECK IN LISTS SCALAR_CHECK_TYPE)
            LIST(FIND ${SCALAR_VALUES} ${SCALAR_CHECK} IDX1)
            IF (IDX1 GREATER -1)
                LIST(FIND SCALAR_CHECK_TYPE ${SCALAR_CHECK} IDX)
                LIST(GET CHECK_SCALAR_FLAG ${IDX} FLAG)

                MATH( EXPR IDX2 "${IDX1} + 1")
                LIST(GET ${SCALAR_VALUES} ${IDX2} VALUE)

                SET( TEST_NAME "${FULL_NAME}-${SCALAR_CHECK}" )
                #MESSAGE("Adding scalar check ${TEST_NAME}")
                SET(CHECK_CMD ${CMAKE_SOURCE_DIR}/utils/check_scalars.py --ns 3 -p ${PREFIX} -e 2 ${FLAG} ${VALUE})
                #MESSAGE("check command = ${CHECK_CMD}")
                ADD_TEST( NAME ${TEST_NAME}
                    COMMAND ${CHECK_CMD}
                    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${FULL_NAME}"
                )
                SET_PROPERTY( TEST ${TEST_NAME} APPEND PROPERTY DEPENDS ${FULL_NAME} )
            ENDIF()
        ENDFOREACH(SCALAR_CHECK)
    ENDIF()
ENDFUNCTION()

IF(NOT QMC_CUDA)

IF (NOT QMC_COMPLEX)
#
# LiH molecular dimer gaussian tests, all electron
# Also check results for different number of mpi tasks and threads keeping total constant
#
  LIST(APPEND LIH_SCALARS "kinetic" "7.9914014383 0.03") # kinetic energy
  LIST(APPEND LIH_SCALARS "totenergy" "-7.9873236457 0.003") # total energy
  LIST(APPEND LIH_SCALARS "eeenergy" "3.4888166386 0.0024") # e-e energy
  LIST(APPEND LIH_SCALARS "samples" "1600000 0.0") # samples
  LIST(APPEND LIH_SCALARS "flux" "0.0 0.1") # sampling check, should be zero

  QMC_RUN_AND_CHECK(short-LiH_dimer_ae-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml
                    16 1
                    LIH_SCALARS
                    TRUE)

  QMC_RUN_AND_CHECK(short-LiH_dimer_ae-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml
                    4 4
                    LIH_SCALARS
                    TRUE)

  QMC_RUN_AND_CHECK(short-LiH_dimer_ae-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml
                    1 16
                    LIH_SCALARS
                    TRUE)

#
# LiH molecular dimer pseudopotential tests
#
  LIST(APPEND LIH_PP_SCALARS "kinetic" "0.6354491487 0.0013")
  LIST(APPEND LIH_PP_SCALARS "totenergy" "-0.7507264125 0.00043")
  LIST(APPEND LIH_PP_SCALARS "eeenergy" "0.4814819907 0.00061")
  LIST(APPEND LIH_PP_SCALARS "potential" "-1.3861755612 0.0013")
  LIST(APPEND LIH_PP_SCALARS "samples" "1600000 0.0")
  LIST(APPEND LIH_PP_SCALARS "flux" "0.0 0.0025")

  QMC_RUN_AND_CHECK(short-LiH_dimer_pp-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml
                    1 16
                    LIH_PP_SCALARS
                    TRUE)

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
  LIST(APPEND LIH_GAMMA_SCALARS "kinetic" "7.4145221704067854 0.024")
  LIST(APPEND LIH_GAMMA_SCALARS "localpotential" "-11.718808431876706 0.025")
  LIST(APPEND LIH_GAMMA_SCALARS "samples" "1600000 0.0")
  LIST(APPEND LIH_GAMMA_SCALARS "flux" "0.0 0.05")

  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-gamma
                    hf_vmc_LiH-gamma.xml
                    1 16
                    LIH_GAMMA_SCALARS
                    TRUE)


  LIST(APPEND LIH_X_SCALARS "kinetic" "7.3772282047321536 0.024")
  LIST(APPEND LIH_X_SCALARS "localpotential" "-11.414410381100344 0.025")
  LIST(APPEND LIH_X_SCALARS "samples" "1600000 0.0")
  LIST(APPEND LIH_X_SCALARS "flux" "0.0 0.05")

  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-x-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-x
                    hf_vmc_LiH-x.xml
                    1 16
                    LIH_X_SCALARS
                    TRUE)


# Arbitrary k-point requires complex build to run. Ensure fails with real/non-complex build
  LIST(APPEND LIH_ARB_SCALARS "kinetic" "7.4201207734399359 0.024")
  LIST(APPEND LIH_ARB_SCALARS "localpotential" "-11.615276430529056 0.025")
  LIST(APPEND LIH_ARB_SCALARS "samples" "1600000 0.0") # samples
  LIST(APPEND LIH_ARB_SCALARS "flux" "0.0 0.05")

  IF (QMC_COMPLEX)
    SET(LIH_ARB_SUCCEED TRUE)
  ELSE()
    SET(LIH_ARB_SUCCEED FALSE)
  ENDIF()


  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-arb-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-arb
                    hf_vmc_LiH-arb.xml
                    1 16
                    LIH_ARB_SCALARS
                    ${LIH_ARB_SUCCEED})

#
# Short tests, about 2-3 seconds on 16 cores
#
  LIST(APPEND BCC_H_SCALARS "totenergy" "-1.834032 0.0003")
  LIST(APPEND BCC_H_SCALARS "samples" "960000 0.0")
  LIST(APPEND BCC_H_SCALARS "flux" "0.0 0.03")

  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    BCC_H_SCALARS
                    TRUE)


  LIST(APPEND DIAMOND_SCALARS "totenergy" "-10.491445 0.002")
  LIST(APPEND DIAMOND_SCALARS "samples" "128000 0.0")
  LIST(APPEND DIAMOND_SCALARS "flux" "0.0 0.4")

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    DIAMOND_SCALARS
                    TRUE)


  LIST(APPEND DIAMOND2_SCALARS "totenergy" "-21.668293 0.007")
  LIST(APPEND DIAMOND2_SCALARS "samples" "32000 0.0")
  LIST(APPEND DIAMOND2_SCALARS "flux" "0.0 0.7")

  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    DIAMOND2_SCALARS
                    TRUE)


  LIST(APPEND HCP_BE_SCALARS "totenergy" "-1.481656 0.0005")
  LIST(APPEND HCP_BE_SCALARS "samples" "416000 0.0")
  LIST(APPEND HCP_BE_SCALARS "flux" "0.0 0.003")

  QMC_RUN_AND_CHECK(short-hcpBe_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/hcpBe_1x1x1_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    HCP_BE_SCALARS
                    TRUE)


  LIST(APPEND MONO_O_SCALARS "totenergy" "-31.776596 0.006")
  LIST(APPEND MONO_O_SCALARS "samples" "64000 0.0")
  LIST(APPEND HCP_BE_SCALARS "flux" "0.0 0.5")

  QMC_RUN_AND_CHECK(short-monoO_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    MONO_O_SCALARS
                    TRUE)



#
# Long variants of above, about 1/2 hour with 16 cores
#
  LIST(APPEND LONG_BCC_H_SCALARS "totenergy" "-1.834032 0.00001")
  LIST(APPEND LONG_BCC_H_SCALARS "samples" "921600000 0.0")

  QMC_RUN_AND_CHECK(long-bccH_1x1x1_ae-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    LONG_BCC_H_SCALARS
                    TRUE)


  LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-10.491445 0.000065")
  LIST(APPEND LONG_DIAMOND_SCALARS "samples" "122880000 0.0")
  LIST(APPEND LONG_DIAMOND_SCALARS "flux" "0.0 0.03")

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    LONG_DIAMOND_SCALARS
                    TRUE)
  

  LIST(APPEND LONG_DIAMOND2_SCALARS "totenergy" "-21.668293 0.00023")
  LIST(APPEND LONG_DIAMOND2_SCALARS "samples" "30720000 0.0")
  LIST(APPEND LONG_DIAMOND2_SCALARS "flux" "0.0 0.05")

  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    LONG_DIAMOND2_SCALARS
                    TRUE)


  LIST(APPEND LONG_HCP_BE_SCALARS "totenergy" "-1.481656 0.000017")
  LIST(APPEND LONG_HCP_BE_SCALARS "samples" "399360000 0.0")
  LIST(APPEND LONG_HCP_BE_SCALARS "flux" "0.0 0.0001")

  QMC_RUN_AND_CHECK(long-hcpBe_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/hcpBe_1x1x1_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    LONG_HCP_BE_SCALARS
                    TRUE)


  LIST(APPEND LONG_MONO_O_SCALARS "totenergy" "-31.776596 0.0002")
  LIST(APPEND LONG_MONO_O_SCALARS "samples" "61440000 0.0")
  LIST(APPEND LONG_MONO_O_SCALARS "flux" "0.0 0.05")

  QMC_RUN_AND_CHECK(long-monoO_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    LONG_MONO_O_SCALARS
                    TRUE)

endif()
