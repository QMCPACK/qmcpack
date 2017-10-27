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



#
# Homogeneous electron gas (HEG) tests
#
IF (NOT QMC_CUDA)

IF (NOT QMC_COMPLEX)
#
# Tests for HEG with 14 electrons at the gamma point and various wavefunctions
# See tests/heg/heg_14_gamma/README for full details
#

#
# Short HEG 14 gamma tests
#

#
#  HEG14G - non-interacting Slater only VMC (exact)
#
  LIST(APPEND HEG14GSNI_SCALARS "totenergy" " 0.627711  0.000001") # total energy
  LIST(APPEND HEG14GSNI_SCALARS "variance"  " 0.000000  0.000001") # energy variance

  QMC_RUN_AND_CHECK(short-heg_14_gamma-ni
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_NI
                    heg-short-NI.xml
                    1 16
                    TRUE
                    0 HEG14GSNI_SCALARS
                    )

#
#  HEG14G - non-interacting Slater only DMC (exact)
#
  LIST(APPEND HEG14GSNI_DMC_SCALARS "totenergy" " 0.627711  0.000001") # total energy
  LIST(APPEND HEG14GSNI_DMC_SCALARS "variance"  " 0.000000  0.000001") # energy variance

  QMC_RUN_AND_CHECK(short-heg_14_gamma-ni_dmc
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_NI_dmc
                    heg-short-NI-dmc.xml
                    1 16
                    TRUE
                    1 HEG14GSNI_DMC_SCALARS # DMC
                    )

#
#  HEG14G - Slater only VMC (Hartree-Fock)
#
  LIST(APPEND HEG14GSHF_SCALARS "totenergy" "-0.812484  0.000271") # total energy
  LIST(APPEND HEG14GSHF_SCALARS "variance"  " 0.193169  0.018612") # energy variance

  QMC_RUN_AND_CHECK(short-heg_14_gamma-hf
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_HF
                    heg-short-HF.xml
                    1 16
                    TRUE
                    0 HEG14GSHF_SCALARS # VMC
                    )

#
#  HEG14G - Slater-Jastrow VMC
#
  LIST(APPEND HEG14GSSJ_SCALARS "totenergy" "-1.073323  0.000271") # total energy
  LIST(APPEND HEG14GSSJ_SCALARS "variance"  " 0.024574  0.000181") # energy variance

  QMC_RUN_AND_CHECK(short-heg_14_gamma-sj
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_SJ
                    heg-short-SJ.xml
                    1 16
                    TRUE
                    0 HEG14GSSJ_SCALARS # VMC
                    )

#
#  HEG14G - Slater-Jastrow VMC with newer input format
#
  LIST(APPEND HEG14GSSJ_NEW_SCALARS "totenergy" "-1.073286  0.000271") # total energy
  LIST(APPEND HEG14GSSJ_NEW_SCALARS "variance"  " 0.024564  0.000342") # energy variance

  QMC_RUN_AND_CHECK(short-heg_14_gamma-sj_new
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_SJ_new
                    heg-short-SJ-new.xml
                    1 16
                    TRUE
                    0 HEG14GSSJ_NEW_SCALARS # VMC
                    )

#
#  HEG14G - Slater-Jastrow-Backflow VMC
#           Run 16x1, 4x4, 1x16 to test SJB cloning
#
  LIST(APPEND HEG14GSSJB_SCALARS "totenergy" "-1.084963  0.000251") # total energy
  LIST(APPEND HEG14GSSJB_SCALARS "variance"  " 0.022667  0.000764") # energy variance

  QMC_RUN_AND_CHECK(short-heg_14_gamma-sjb
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_SJB
                    heg-short-SJB.xml
                    16 1
                    TRUE
                    0 HEG14GSSJB_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-heg_14_gamma-sjb
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_SJB
                    heg-short-SJB.xml
                    4 4
                    TRUE
                    0 HEG14GSSJB_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-heg_14_gamma-sjb
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_SJB
                    heg-short-SJB.xml
                    1 16
                    TRUE
                    0 HEG14GSSJB_SCALARS # VMC
                    )

#
#  HEG14G - Slater-Jastrow DMC
#
  LIST(APPEND HEG14GSSJ_DMC_SCALARS "totenergy" "-1.110199  0.000322") # total energy
  LIST(APPEND HEG14GSSJ_DMC_SCALARS "variance"  " 0.022894  0.000060") # energy variance

  QMC_RUN_AND_CHECK(short-heg_14_gamma-sj_dmc
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_short_SJ_dmc
                    heg-short-SJ-dmc.xml
                    1 16
                    TRUE
                    2 HEG14GSSJ_DMC_SCALARS # DMC
                    )

#
# Long HEG 14 gamma tests
#

#
#  HEG14G - non-interacting Slater only VMC (exact)
#
  LIST(APPEND HEG14GLNI_SCALARS "totenergy" " 0.627711  0.000001") # total energy
  LIST(APPEND HEG14GLNI_SCALARS "variance"  " 0.000000  0.000001") # energy variance

  QMC_RUN_AND_CHECK(long-heg_14_gamma-ni
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_NI
                    heg-long-NI.xml
                    1 16
                    TRUE
                    0 HEG14GLNI_SCALARS # VMC
                    )

#
#  HEG14G - non-interacting Slater only DMC (exact)
#
  LIST(APPEND HEG14GLNI_DMC_SCALARS "totenergy" " 0.627711  0.000001") # total energy
  LIST(APPEND HEG14GLNI_DMC_SCALARS "variance"  " 0.000000  0.000001") # energy variance

  QMC_RUN_AND_CHECK(long-heg_14_gamma-ni_dmc
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_NI_dmc
                    heg-long-NI-dmc.xml
                    1 16
                    TRUE
                    1 HEG14GLNI_DMC_SCALARS # DMC
                    )

#
#  HEG14G - Slater only VMC (Hartree-Fock)
#
  LIST(APPEND HEG14GLHF_SCALARS "totenergy" "-0.812484  0.000090") # total energy
  LIST(APPEND HEG14GLHF_SCALARS "variance"  " 0.193169  0.006142") # energy variance

  QMC_RUN_AND_CHECK(long-heg_14_gamma-hf
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_HF
                    heg-long-HF.xml
                    1 16
                    TRUE
                    0 HEG14GLHF_SCALARS # VMC
                    )

#
#  HEG14G - Slater-Jastrow VMC
#
  LIST(APPEND HEG14GLSJ_SCALARS "totenergy" "-1.073323  0.000090") # total energy
  LIST(APPEND HEG14GLSJ_SCALARS "variance"  " 0.024574  0.000060") # energy variance

  QMC_RUN_AND_CHECK(long-heg_14_gamma-sj
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_SJ
                    heg-long-SJ.xml
                    1 16
                    TRUE
                    0 HEG14GLSJ_SCALARS # VMC
                    )

#
#  HEG14G - Slater-Jastrow VMC with newer input format
#
  LIST(APPEND HEG14GLSJ_NEW_SCALARS "totenergy" "-1.073286  0.000090") # total energy
  LIST(APPEND HEG14GLSJ_NEW_SCALARS "variance"  " 0.024564  0.000113") # energy variance

  QMC_RUN_AND_CHECK(long-heg_14_gamma-sj_new
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_SJ_new
                    heg-long-SJ-new.xml
                    1 16
                    TRUE
                    0 HEG14GLSJ_NEW_SCALARS # VMC
                    )

#
#  HEG14G - Slater-Jastrow-Backflow VMC
#           Run 16x1, 4x4, 1x16 to test SJB cloning
#
  LIST(APPEND HEG14GLSJB_SCALARS "totenergy" "-1.084963  0.000083") # total energy
  LIST(APPEND HEG14GLSJB_SCALARS "variance"  " 0.022667  0.000252") # energy variance

  QMC_RUN_AND_CHECK(long-heg_14_gamma-sjb
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_SJB
                    heg-long-SJB.xml
                    16 1
                    TRUE
                    0 HEG14GLSJB_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(long-heg_14_gamma-sjb
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_SJB
                    heg-long-SJB.xml
                    4 4
                    TRUE
                    0 HEG14GLSJB_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(long-heg_14_gamma-sjb
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_SJB
                    heg-long-SJB.xml
                    1 16
                    TRUE
                    0 HEG14GLSJB_SCALARS # VMC
                    )

#
#  HEG14G - Slater-Jastrow DMC
#
  LIST(APPEND HEG14GLSJ_DMC_SCALARS "totenergy" "-1.110199  0.000106") # total energy
  LIST(APPEND HEG14GLSJ_DMC_SCALARS "variance"  " 0.022894  0.000020") # energy variance

  QMC_RUN_AND_CHECK(long-heg_14_gamma-sj_dmc
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_14_gamma"
                    heg_long_SJ_dmc
                    heg-long-SJ-dmc.xml
                    1 16
                    TRUE
                    2 HEG14GLSJ_DMC_SCALARS # DMC
                    )

 ELSE()
    MESSAGE("Skipping all HEG tests because they are not supported by complex build (QMC_COMPLEX=1)")  
 ENDIF()
ELSE()
    MESSAGE("Skipping all HEG tests because they are not supported by CUDA build (QMC_CUDA=1)")  
ENDIF()



IF (NOT QMC_CUDA)

IF (NOT QMC_COMPLEX)
#
# Li2 molecular dimer STO tests, all electron
#   Reproduces UNR 1993
#

# Short tests
  LIST(APPEND LI2STOSJS_VMC_SCALARS "totenergy" "-14.947214  0.000502") # total energy
  LIST(APPEND LI2STOSJS_VMC_SCALARS "variance"  "  0.170618  0.001035") # energy variance
  LIST(APPEND LI2STOSJS_DMC_SCALARS "totenergy" "-14.989438  0.000643") # total energy
  LIST(APPEND LI2STOSJS_DMC_SCALARS "variance"  "  0.151118  0.000312") # energy variance

  QMC_RUN_AND_CHECK(short-li2_sto-sj_dmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/Li2_STO_ae"
                    Li2.STO.short
                    Li2.STO.short.in.xml
                    1 16
                    TRUE
                    0 LI2STOSJS_VMC_SCALARS # VMC
                    2 LI2STOSJS_DMC_SCALARS # DMC
                    )

  QMC_RUN_AND_CHECK(short-li2_sto-sj_dmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/Li2_STO_ae"
                    Li2.STO.short
                    Li2.STO.short.in.xml
                    4 4
                    TRUE
                    0 LI2STOSJS_VMC_SCALARS # VMC
                    2 LI2STOSJS_DMC_SCALARS # DMC
                    )

# Long tests
  LIST(APPEND LI2STOSJL_VMC_SCALARS "totenergy" "-14.947214  0.000166") # total energy
  LIST(APPEND LI2STOSJL_VMC_SCALARS "variance"  "  0.170618  0.000342") # energy variance
  LIST(APPEND LI2STOSJL_DMC_SCALARS "totenergy" "-14.989438  0.000212") # total energy
  LIST(APPEND LI2STOSJL_DMC_SCALARS "variance"  "  0.151118  0.000103") # energy variance

  QMC_RUN_AND_CHECK(long-li2_sto-sj_dmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/Li2_STO_ae"
                    Li2.STO.long
                    Li2.STO.long.in.xml
                    1 16
                    TRUE
                    0 LI2STOSJL_VMC_SCALARS # VMC
                    2 LI2STOSJL_DMC_SCALARS # DMC
                    )

# C2 TZ PP multi determinants
# Short tests
  LIST(APPEND C2MSDJS_VMC_SCALARS "totenergy" "-11.064964  0.000865") # total energy
  LIST(APPEND C2MSDJS_VMC_SCALARS "variance"  "  0.149682  0.003928") # energy variance

  QMC_RUN_AND_CHECK(short-C2_pp-msdj_vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C2_pp"
                    C2_TZ_PP_MSDJ
                    qmc_short_pp_msdj.xml
                    1 16
                    TRUE
                    1 C2MSDJS_VMC_SCALARS # VMC
                    )

  LIST(APPEND C2MSDJL_VMC_SCALARS "totenergy" "-11.064964  0.000274") # total energy
  LIST(APPEND C2MSDJL_VMC_SCALARS "variance"  "  0.149682  0.001242") # energy variance

  QMC_RUN_AND_CHECK(long-C2_pp-msdj_vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C2_pp"
                    C2_TZ_PP_MSDJ
                    qmc_long_pp_msdj.xml
                    1 16
                    TRUE
                    1 C2MSDJL_VMC_SCALARS # VMC
                    )

#
# LiH molecular dimer gaussian tests, all electron
# Also check results for different number of mpi tasks and threads keeping total constant
#
  LIST(APPEND LIH_SCALARS "kinetic" "7.9914014383 0.03") # kinetic energy
  LIST(APPEND LIH_SCALARS "totenergy" "-7.9873236457 0.003") # total energy
  LIST(APPEND LIH_SCALARS "eeenergy" "3.4888166386 0.0024") # e-e energy
  LIST(APPEND LIH_SCALARS "samples" "1600000 0.0") # samples
#  LIST(APPEND LIH_SCALARS "flux" "0.0 0.1") # sampling check, should be zero

  QMC_RUN_AND_CHECK(short-LiH_dimer_ae-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml
                    16 1
                    TRUE
                    0 LIH_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-LiH_dimer_ae-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml
                    4 4
                    TRUE
                    0 LIH_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-LiH_dimer_ae-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml
                    1 16
                    TRUE
                    0 LIH_SCALARS # VMC
                    )

#
# LiH molecular dimer pseudopotential tests
#
  LIST(APPEND LIH_PP_SCALARS "kinetic" "0.6354491487 0.0013")
  LIST(APPEND LIH_PP_SCALARS "totenergy" "-0.7507264125 0.00043")
  LIST(APPEND LIH_PP_SCALARS "eeenergy" "0.4814819907 0.00061")
  LIST(APPEND LIH_PP_SCALARS "potential" "-1.3861755612 0.0013")
  LIST(APPEND LIH_PP_SCALARS "samples" "1600000 0.0")
#  LIST(APPEND LIH_PP_SCALARS "flux" "0.0 0.0025")

  QMC_RUN_AND_CHECK(short-LiH_dimer_pp-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml
                    1 16
                    TRUE
                    0 LIH_PP_SCALARS # VMC
                    )

#
# 4 H atom cluster gaussian tests, all electron
# Also check results for different number of mpi tasks and threads keeping total constant
#
  LIST(APPEND H4_OPT_SCALARS "totenergy" "-2.140786 0.0009") # total energy

  QMC_RUN_AND_CHECK(short-H4-opt-linear
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-linear
                    optm-linear.xml
                    16 1
                    TRUE
                    5 H4_OPT_SCALARS # OPT step 5
                    )

IF(BUILD_LMYENGINE_INTERFACE)
  QMC_RUN_AND_CHECK(short-H4-opt-adaptive
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-adaptive
                    optm-adaptive.xml
                    16 1
                    TRUE
                    5 H4_OPT_SCALARS # OPT step 5
                    )

  QMC_RUN_AND_CHECK(short-H4-opt-adaptive
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-adaptive
                    optm-adaptive.xml
                    1 16
                    TRUE
                    5 H4_OPT_SCALARS # OPT step 5
                    )
ENDIF(BUILD_LMYENGINE_INTERFACE)

  QMC_RUN_AND_CHECK(short-H4-opt-OneShiftOnly
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-OneShiftOnly
                    optm-OneShiftOnly.xml
                    16 1
                    TRUE
                    5 H4_OPT_SCALARS # OPT step 5
                    )

  QMC_RUN_AND_CHECK(short-H4-opt-OneShiftOnly
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-OneShiftOnly
                    optm-OneShiftOnly.xml
                    1 16
                    TRUE
                    5 H4_OPT_SCALARS # OPT step 5
                    )

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
#  LIST(APPEND LIH_GAMMA_SCALARS "flux" "0.0 0.05")

  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-gamma
                    hf_vmc_LiH-gamma.xml
                    1 16
                    TRUE
                    0 LIH_GAMMA_SCALARS # VMC
                    )


  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-gamma-drift-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-gamma-drift
                    hf_vmc_LiH-gamma-drift.xml
                    1 16
                    TRUE
                    0 LIH_GAMMA_SCALARS # VMC
                    )


  LIST(APPEND LIH_X_SCALARS "kinetic" "7.3772282047321536 0.024")
  LIST(APPEND LIH_X_SCALARS "localecp" "-11.414410381100344 0.025")
  LIST(APPEND LIH_X_SCALARS "samples" "1600000 0.0")
#  LIST(APPEND LIH_X_SCALARS "flux" "0.0 0.05")

  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-x-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-x
                    hf_vmc_LiH-x.xml
                    1 16
                    TRUE
                    0 LIH_X_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-x-drift-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-x-drift
                    hf_vmc_LiH-x-drift.xml
                    1 16
                    TRUE
                    0 LIH_X_SCALARS # VMC
                    )

IF(ENABLE_SOA)
  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_hybridrep_pp-x-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-x
                    hf_vmc_LiH-x_hybridrep.xml
                    1 16
                    TRUE
                    0 LIH_X_SCALARS # VMC
                    )
ENDIF(ENABLE_SOA)

# Arbitrary k-point requires complex build to run. Ensure fails with real/non-complex build
  LIST(APPEND LIH_ARB_SCALARS "kinetic" "7.4201207734399359 0.024")
  LIST(APPEND LIH_ARB_SCALARS "localecp" "-11.615276430529056 0.025")
  LIST(APPEND LIH_ARB_SCALARS "samples" "1600000 0.0") # samples
#  LIST(APPEND LIH_ARB_SCALARS "flux" "0.0 0.05")

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
                    ${LIH_ARB_SUCCEED}
                    0 LIH_ARB_SCALARS # VMC
                    )


  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-arb-drift-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-arb-drift
                    hf_vmc_LiH-arb-drift.xml
                    1 16
                    ${LIH_ARB_SUCCEED}
                    0 LIH_ARB_SCALARS # VMC
                    )

# DMC tests for LiH, locality and non-local (t-moves) variants

  LIST(APPEND LIH_GAMMA_DMC_SCALARS "totenergy" "-8.54823 0.004")
  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-gamma-dmc-hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_dmc_short_LiH-gamma
                    hf_vmc_dmc_short_LiH-gamma.xml
                    1 16
                    TRUE
                    1 LIH_GAMMA_DMC_SCALARS # DMC
                    )

  LIST(APPEND LIH_X_DMC_SCALARS "totenergy" "-8.31413 0.004")
  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-x-dmc-hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_dmc_short_LiH-x
                    hf_vmc_dmc_short_LiH-x.xml
                    1 16
                    TRUE
                    1 LIH_X_DMC_SCALARS # DMC
                    )

  LIST(APPEND LIH_ARB_DMC_SCALARS "totenergy" "-8.45401 0.004")
  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-arb-dmc-hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_dmc_short_LiH-arb
                    hf_vmc_dmc_short_LiH-arb.xml
                    1 16
                    ${LIH_ARB_SUCCEED}
                    1 LIH_ARB_DMC_SCALARS # DMC
                    )

  LIST(APPEND LIH_GAMMA_DMCNL_SCALARS "totenergy" "-8.54979 0.004")
  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-gamma-dmcnl-hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_dmcnl_short_LiH-gamma
                    hf_vmc_dmcnl_short_LiH-gamma.xml
                    1 16
                    TRUE
                    1 LIH_GAMMA_DMCNL_SCALARS # DMC
                    )

  LIST(APPEND LIH_X_DMCNL_SCALARS "totenergy" "-8.31369 0.004")
  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-x-dmcnl-hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_dmcnl_short_LiH-x
                    hf_vmc_dmcnl_short_LiH-x.xml
                    1 16
                    TRUE
                    1 LIH_X_DMCNL_SCALARS # DMC
                    )

  LIST(APPEND LIH_ARB_DMCNL_SCALARS "totenergy" "-8.45174 0.004")
  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-arb-dmcnl-hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_dmcnl_short_LiH-arb
                    hf_vmc_dmcnl_short_LiH-arb.xml
                    1 16
                    ${LIH_ARB_SUCCEED}
                    1 LIH_ARB_DMCNL_SCALARS # DMC
                    )


#
# Short tests, about 2-3 seconds on 16 cores
#
  LIST(APPEND BCC_H_SCALARS "totenergy" "-1.834032 0.0003")
  LIST(APPEND BCC_H_SCALARS "kinetic" "0.186673807952 0.003")
  LIST(APPEND BCC_H_SCALARS "potential" "-2.02070672998 0.003")
  LIST(APPEND BCC_H_SCALARS "eeenergy" "-0.775870180954 0.003")
  LIST(APPEND BCC_H_SCALARS "ionion" "-0.96289961022 0.001")
  LIST(APPEND BCC_H_SCALARS "samples" "960000 0.0")
#  LIST(APPEND BCC_H_SCALARS "flux" "0.0 0.03")

  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 BCC_H_SCALARS # VMC
                    )

IF(NOT MIXED_PRECISION)
  SIMPLE_RUN_AND_CHECK(short-bccH_3x3x3_ae-deriv
    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_3x3x3_ae"
    deriv.xml
    1 1
    check_deriv.py
  )

  SIMPLE_RUN_AND_CHECK(short-bccH_3x3x3_ae-grad_lap
    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_3x3x3_ae"
    grad_lap.xml
    1 1
    check_grad_lap.py
  )

  SIMPLE_RUN_AND_CHECK(short-bccH_3x3x3_ae-gamma-deriv
    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_3x3x3_ae"
    gamma_deriv.xml
    1 1
    check_deriv.py
  )

  SIMPLE_RUN_AND_CHECK(short-bccH_3x3x3_ae-not_orth-deriv
    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_3x3x3_ae"
    not-orth_deriv.xml
    1 1
    check_deriv.py
  )

  SIMPLE_RUN_AND_CHECK(short-bccH_2x2x2_ae-grad-lap
    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_2x2x2_ae"
    grad_lap.xml
    1 1
    ../short-bccH_3x3x3_ae-deriv-1-1/check_grad_lap.py
  )

  SIMPLE_RUN_AND_CHECK(short-bccH_2x2x2_ae-deriv
    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_2x2x2_ae"
    deriv.xml
    1 1
    ../short-bccH_3x3x3_ae-deriv-1-1/check_deriv.py
  )

  SIMPLE_RUN_AND_CHECK(short-bccH_2x2x2_ae-gamma-deriv
    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_2x2x2_ae"
    gamma_deriv.xml
    1 1
    ../short-bccH_3x3x3_ae-deriv-1-1/check_deriv.py
  )
ENDIF(NOT MIXED_PRECISION)

# qmc_ref_vmc_dmc  series 1  -1.841352 +/- 0.000078
# qmc_short_vmc_dmc  series 1  -1.84196 +/- 0.00073
  LIST(APPEND BCC_H_DMC_SCALARS "totenergy" "-1.84196 0.00073")
  LIST(APPEND BCC_H_DMC_SCALARS "potential" "-2.0187  0.0033")
  LIST(APPEND BCC_H_DMC_SCALARS "kinetic"    "0.1767  0.0035")
  LIST(APPEND BCC_H_DMC_SCALARS "eeenergy"  "-0.77985 0.00069")
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_vmc_dmc
                    qmc_short_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 BCC_H_DMC_SCALARS # DMC
                    )

  IF(NOT QMC_CUDA)

# Check VMCUpdateAll driver using all-particle updates
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-vmc-all-nodrift_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_vmc_all_nodrift
                    qmc_short_vmc_all_nodrift.in.xml
                    1 16
                    TRUE
                    0 BCC_H_SCALARS # VMC
                    )

# Check VMCUpdateAllWithDrift driver using all-particle updates
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-vmc-all_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_vmc_all
                    qmc_short_vmc_all.in.xml
                    1 16
                    TRUE
                    0 BCC_H_SCALARS # VMC
                    )

# Check DMCUpdateAllWithRejection driver using all-particle updates
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-dmc-all_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_vmc_dmc_all
                    qmc_short_vmc_dmc_all.in.xml
                    1 16
                    TRUE
                    1 BCC_H_DMC_SCALARS # DMC
                    )

#Now for the correlated sampling tests
  LIST(APPEND BCC_H_CSVMC_SCALARS "totenergy_A" "-1.83428 0.00013")
  LIST(APPEND BCC_H_CSVMC_SCALARS "totenergy_B" "-1.83390 0.00012")
  LIST(APPEND BCC_H_CSVMC_SCALARS "dtotenergy_AB" "-0.000384 0.000068")
  LIST(APPEND BCC_H_CSVMC_SCALARS "ionion_A" "-0.96289962 0.0000001")
  LIST(APPEND BCC_H_CSVMC_SCALARS "ionion_B" "-0.95816559 0.0000001")
  LIST(APPEND BCC_H_CSVMC_SCALARS "dionion_AB" "-0.004734021 0.0000001")
  LIST(APPEND BCC_H_CSVMC_SCALARS "eeenergy_A" "-0.77607 0.00013")
  LIST(APPEND BCC_H_CSVMC_SCALARS "eeenergy_B" "-0.77546 0.00013")
  LIST(APPEND BCC_H_CSVMC_SCALARS "deeenergy_AB" "-0.000603 0.000035")

# Check CSVMCUpdatePbyP driver using particle-by-particle updates
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-csvmc-pbyp-nodrift_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_csvmc_pbyp_nodrift
                    qmc_short_csvmc_pbyp_nodrift.in.xml
                    16 1
                    TRUE
                    0 BCC_H_CSVMC_SCALARS # VMC
                    )

# Check CSVMCUpdatePbyPWithDrift driver using all-particle updates
#  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-csvmc-pbyp_sdj
#                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
#                    qmc_short_csvmc_pbyp
#                    qmc_short_csvmc_pbyp.in.xml
#                    16 1
#                    TRUE
#                    0 BCC_H_CSVMC_SCALARS # VMC
#                    )

# Check CSVMCUpdateAll driver using all-particle updates
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-csvmc-all-nodrift_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_csvmc_all_nodrift
                    qmc_short_csvmc_all_nodrift.in.xml
                    16 1
                    TRUE
                    0 BCC_H_CSVMC_SCALARS # VMC
                    )

# Check CSVMCUpdateAllWithDrift driver using all-particle updates
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-csvmc-all_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_csvmc_all
                    qmc_short_csvmc_all.in.xml
                    16 1
                    TRUE
                    0 BCC_H_CSVMC_SCALARS # VMC
                    )


  LIST(APPEND BCC_H_RMC_SCALARS "totenergy" "-1.8420 0.0011")
  #check RMC on 16 cores.  1 MPI
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-rmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_vmc_rmc
                    qmc_short_vmc_rmc.in.xml
                    1 16
                    TRUE
                    1 BCC_H_RMC_SCALARS #series# 
                    )
  #check RMC on 4 cores, 4 mpi
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-rmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_vmc_rmc
                    qmc_short_vmc_rmc.in.xml
                    4 4
                    TRUE
                    1 BCC_H_RMC_SCALARS #series#
                    )
  ENDIF()

                     
# Plane wave test of bccH (no pp, so faster than typical solids)
# Reduced number of blocks and samples compared to other short bccH tests
  LIST(APPEND BCC_H_PW_SCALARS "totenergy" "-1.834032 0.001")
  LIST(APPEND BCC_H_PW_SCALARS "samples" "96000 0.0")
#  LIST(APPEND BCC_H_PW_SCALARS "flux" "0.0 0.1")
  IF(NOT QMC_CUDA)
# No GPU implementation. Would be fast if implemented
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-vmc_sdj-pw
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_pw
                    qmc_short_pw.in.xml
                    1 16
                    TRUE
                    0 BCC_H_PW_SCALARS # VMC
                    )
  ENDIF()

  LIST(APPEND DIAMOND_SCALARS "totenergy" "-10.491445 0.002")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "11.434392 0.018")
  LIST(APPEND DIAMOND_SCALARS "potential" "-21.925822 0.018")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-2.700534 0.0043")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-12.77566752 0.001")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-7.046740 0.020")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "0.597119 0.0056")
  LIST(APPEND DIAMOND_SCALARS "samples" "128000 0.0")
#  LIST(APPEND DIAMOND_SCALARS "flux" "0.0 0.4")

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )

IF(ENABLE_SOA)
  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_hybridrep_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short
                    qmc_short_hybridrep.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )
ENDIF(ENABLE_SOA)

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_sdj-meshf
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short-meshf
                    qmc_short-meshf.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )

# Reference OPT run in qmc-ref
  LIST(APPEND DIAMOND_OPT_SCALARS "totenergy"   "-10.49370 0.0024")
  LIST(APPEND DIAMOND_OPT_SCALARS "kinetic"     "11.5504 0.021")
  LIST(APPEND DIAMOND_OPT_SCALARS "localecp"    "-7.1979 0.023")
  LIST(APPEND DIAMOND_OPT_SCALARS "nonlocalecp" "0.61540 0.0068")
  LIST(APPEND DIAMOND_OPT_SCALARS "samples"     "128000 0.0")

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-opt_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short
                    qmc_short_opt.in.xml
                    1 16
                    TRUE
                    3 DIAMOND_OPT_SCALARS # VMC
                    )

# Reference DMC run in qmc-ref "-10.531583 0.000265"
  LIST(APPEND DIAMOND_DMC_SCALARS "totenergy" "-10.531583 0.0024")
  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short_vmc_dmc
                    qmc_short_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 DIAMOND_DMC_SCALARS # DMC
                    )

# Coverage test - shorter version of dmc test
  SET(FULL_NAME "coverage-diamondC_1x1x1_pp-vmc_sdj-1-1")
  SET(TEST_ADDED FALSE)
  RUN_QMC_APP(${FULL_NAME}
              "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
              1 1
              TEST_ADDED
              "qmc_short_vmc_dmc_coverage.in.xml")

  IF (TEST_ADDED)
    SET_PROPERTY(TEST ${FULL_NAME} APPEND PROPERTY LABELS "coverage")
  ENDIF()


  LIST(APPEND DIAMOND2_SCALARS "totenergy" "-21.668293 0.007")
  LIST(APPEND DIAMOND2_SCALARS "kinetic" "20.065151 0.044")
  LIST(APPEND DIAMOND2_SCALARS "potential" "-41.733377 0.045")
  LIST(APPEND DIAMOND2_SCALARS "eeenergy" "-5.428561 0.012")
  LIST(APPEND DIAMOND2_SCALARS "ionion" "-25.5513269 0.001")
  LIST(APPEND DIAMOND2_SCALARS "localecp" "-12.29402 0.051")
  LIST(APPEND DIAMOND2_SCALARS "nonlocalecp" "1.540530 0.018")
  LIST(APPEND DIAMOND2_SCALARS "samples" "32000 0.0")
#  LIST(APPEND DIAMOND2_SCALARS "flux" "0.0 0.7")

  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 DIAMOND2_SCALARS # VMC
                    )

IF(ENABLE_SOA)
  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_hybridrep_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_short
                    qmc_short_hybridrep.in.xml
                    1 16
                    TRUE
                    0 DIAMOND2_SCALARS # VMC
                    )
ENDIF(ENABLE_SOA)

# Reference DMC run in qmc-ref "-21.844975 0.001113"
  LIST(APPEND DIAMOND2_DMC_SCALARS "totenergy" "-21.844975 0.00779")
  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_short_vmc_dmc
                    qmc_short_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 DIAMOND2_DMC_SCALARS # DMC
                    )

  LIST(APPEND DIAMOND2_DMC_SCALARS "totenergy" "-21.844975 0.00779")
  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp-dmc-reconf_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_short_vmc_dmc
                    qmc_short_vmc_dmc_reconf.in.xml
                    1 16
                    TRUE
                    1 DIAMOND2_DMC_SCALARS # DMC
                    )


  LIST(APPEND HCP_BE_SCALARS "totenergy" "-1.481656 0.0005")
  LIST(APPEND HCP_BE_SCALARS "samples" "416000 0.0")
#  LIST(APPEND HCP_BE_SCALARS "flux" "0.0 0.003")

  QMC_RUN_AND_CHECK(short-hcpBe_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/hcpBe_1x1x1_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 HCP_BE_SCALARS # VMC
                    )


  LIST(APPEND MONO_O_SCALARS "totenergy" "-31.776596 0.006")
  LIST(APPEND MONO_O_SCALARS "samples" "64000 0.0")
#  LIST(APPEND MONO_O_SCALARS "flux" "0.0 0.5")

  QMC_RUN_AND_CHECK(short-monoO_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 MONO_O_SCALARS # VMC
                    )



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
                    TRUE
                    0 LONG_BCC_H_SCALARS # VMC
                    )


# qmc_ref_vmc_dmc  series 1  -1.841352 +/- 0.000078
# qmc_long_vmc_dmc  series 1   -1.841204 +/- 0.000285
  LIST(APPEND LONG_BCC_H_DMC_SCALARS "totenergy" "-1.841204 0.000285")
  QMC_RUN_AND_CHECK(long-bccH_1x1x1_ae-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 LONG_BCC_H_DMC_SCALARS # DMC
                    )

  IF(NOT QMC_CUDA)
  LIST(APPEND LONG_BCC_H_RMC_SCALARS "totenergy" "-1.84150 0.00011")
  LIST(APPEND LONG_BCC_H_RMC_SCALARS "potential" "-2.01964 0.00072")
  LIST(APPEND LONG_BCC_H_RMC_SCALARS "potential_pure" "-2.01907 0.00085")
  LIST(APPEND LONG_BCC_H_RMC_SCALARS "kinetic_mixed" "0.17814 0.00073")
  LIST(APPEND LONG_BCC_H_RMC_SCALARS "kinetic_pure" "0.17063 0.00086")
  LIST(APPEND LONG_BCC_H_RMC_SCALARS "eeenergy_mixed" "-0.77991 0.00016")
  LIST(APPEND LONG_BCC_H_RMC_SCALARS "eeenergy_pure" "-0.78335 0.00019")
  QMC_RUN_AND_CHECK(long-bccH_1x1x1_ae-rmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_long_vmc_rmc
                    qmc_long_vmc_rmc.in.xml
                    1 16
                    TRUE
                    1 LONG_BCC_H_RMC_SCALARS # DMC
                    )
  ENDIF()

 
  LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-10.491445 0.000065")
  LIST(APPEND LONG_DIAMOND_SCALARS "samples" "122880000 0.0")
  LIST(APPEND LONG_DIAMOND_SCALARS "flux" "0.0 0.03")

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc_sdj-meshf
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_long-meshf
                    qmc_long-meshf.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )

# Reference DMC run in qmc-ref "-10.531583 0.000265"
  LIST(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-10.531583 0.000815")
  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 LONG_DIAMOND_DMC_SCALARS # DMC
                    )


  LIST(APPEND LONG_DIAMOND2_SCALARS "totenergy" "-21.668293 0.00023")
  LIST(APPEND LONG_DIAMOND2_SCALARS "samples" "30720000 0.0")
  LIST(APPEND LONG_DIAMOND2_SCALARS "flux" "0.0 0.05")

  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND2_SCALARS # VMC
                    )

# Reference DMC run in qmc-ref "-21.844975 0.001113"
  LIST(APPEND LONG_DIAMOND2_DMC_SCALARS "totenergy" "-21.844975 0.002879")
  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 LONG_DIAMOND2_DMC_SCALARS # DMC
                    )

  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp-dmc-reconf_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc_reconf.in.xml
                    1 16
                    TRUE
                    1 LONG_DIAMOND2_DMC_SCALARS # DMC
                    )


  LIST(APPEND LONG_HCP_BE_SCALARS "totenergy" "-1.481656 0.000017")
  LIST(APPEND LONG_HCP_BE_SCALARS "samples" "399360000 0.0")
  LIST(APPEND LONG_HCP_BE_SCALARS "flux" "0.0 0.0001")

  QMC_RUN_AND_CHECK(long-hcpBe_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/hcpBe_1x1x1_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_HCP_BE_SCALARS # VMC
                    )


  LIST(APPEND LONG_MONO_O_SCALARS "totenergy" "-31.776596 0.0002")
  LIST(APPEND LONG_MONO_O_SCALARS "samples" "61440000 0.0")
  LIST(APPEND LONG_MONO_O_SCALARS "flux" "0.0 0.05")

  QMC_RUN_AND_CHECK(long-monoO_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_MONO_O_SCALARS # VMC
                    )

endif()
