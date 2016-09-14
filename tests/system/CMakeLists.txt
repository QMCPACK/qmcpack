#############################################################
# Add tests to ctest
#############################################################

message("Adding tests for QMCPACK")

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
                    0 # VMC
                    TRUE)

  QMC_RUN_AND_CHECK(short-LiH_dimer_ae-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml
                    4 4
                    LIH_SCALARS
                    0 # VMC
                    TRUE)

  QMC_RUN_AND_CHECK(short-LiH_dimer_ae-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml
                    1 16
                    LIH_SCALARS
                    0 # VMC
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
                    0 # VMC
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
  LIST(APPEND LIH_GAMMA_SCALARS "localecp" "-11.718808431876706 0.025")
  LIST(APPEND LIH_GAMMA_SCALARS "samples" "1600000 0.0")
  LIST(APPEND LIH_GAMMA_SCALARS "flux" "0.0 0.05")

  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-gamma
                    hf_vmc_LiH-gamma.xml
                    1 16
                    LIH_GAMMA_SCALARS
                    0 # VMC
                    TRUE)


  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-gamma-drift-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-gamma-drift
                    hf_vmc_LiH-gamma-drift.xml
                    1 16
                    LIH_GAMMA_SCALARS
                    0 # VMC
                    TRUE)


  LIST(APPEND LIH_X_SCALARS "kinetic" "7.3772282047321536 0.024")
  LIST(APPEND LIH_X_SCALARS "localecp" "-11.414410381100344 0.025")
  LIST(APPEND LIH_X_SCALARS "samples" "1600000 0.0")
  LIST(APPEND LIH_X_SCALARS "flux" "0.0 0.05")

  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-x-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-x
                    hf_vmc_LiH-x.xml
                    1 16
                    LIH_X_SCALARS
                    0 # VMC
                    TRUE)


# Arbitrary k-point requires complex build to run. Ensure fails with real/non-complex build
  LIST(APPEND LIH_ARB_SCALARS "kinetic" "7.4201207734399359 0.024")
  LIST(APPEND LIH_ARB_SCALARS "localecp" "-11.615276430529056 0.025")
  LIST(APPEND LIH_ARB_SCALARS "samples" "1600000 0.0") # samples

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
                    0 # VMC
                    ${LIH_ARB_SUCCEED})


  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-arb-drift-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-arb-drift
                    hf_vmc_LiH-arb-drift.xml
                    1 16
                    LIH_ARB_SCALARS
                    0 # VMC
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
                    0 # VMC
                    TRUE)

# Plane wave test of bccH (no pp, so faster than typical solids)
# Currently disabled since take 90s to run
#  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-vmc_sdj-pw
#                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
#                    qmc_short_pw
#                    qmc_short_pw.in.xml
#                    1 16
#                    BCC_H_SCALARS
#                    0 # VMC
#                    TRUE)


  LIST(APPEND DIAMOND_SCALARS "totenergy" "-10.491445 0.002")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "11.434392 0.018")
  LIST(APPEND DIAMOND_SCALARS "potential" "-21.925822 0.018")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-2.700534 0.0043")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-12.77566752 0.001")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-7.046740 0.020")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "0.597119 0.0056")
  LIST(APPEND DIAMOND_SCALARS "samples" "128000 0.0")
  LIST(APPEND DIAMOND_SCALARS "flux" "0.0 0.4")

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    DIAMOND_SCALARS
                    0 # VMC
                    TRUE)


  LIST(APPEND DIAMOND2_SCALARS "totenergy" "-21.668293 0.007")
  LIST(APPEND DIAMOND2_SCALARS "kinetic" "20.065151 0.044")
  LIST(APPEND DIAMOND2_SCALARS "potential" "-41.733377 0.045")
  LIST(APPEND DIAMOND2_SCALARS "eeenergy" "-5.428561 0.012")
  LIST(APPEND DIAMOND2_SCALARS "ionion" "-25.5513269 0.001")
  LIST(APPEND DIAMOND2_SCALARS "localecp" "-12.29402 0.051")
  LIST(APPEND DIAMOND2_SCALARS "nonlocalecp" "1.540530 0.018")
  LIST(APPEND DIAMOND2_SCALARS "samples" "32000 0.0")
  LIST(APPEND DIAMOND2_SCALARS "flux" "0.0 0.7")

  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    DIAMOND2_SCALARS
                    0 # VMC
                    TRUE)

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_sdj-meshf
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short-meshf
                    qmc_short-meshf.in.xml
                    1 16
                    DIAMOND_SCALARS
                    0 # VMC
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
                    0 # VMC
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
                    0 # VMC
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
                    0 # VMC
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
                    0 # VMC
                    TRUE)

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc_sdj-meshf
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_long-meshf
                    qmc_long-meshf.in.xml
                    1 16
                    LONG_DIAMOND_SCALARS
                    0 # VMC
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
                    0 # VMC
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
                    0 # VMC
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
                    0 # VMC
                    TRUE)

endif()
