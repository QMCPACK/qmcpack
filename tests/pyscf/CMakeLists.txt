#############################################################
# Add tests to ctest
#############################################################

include("${qmcpack_SOURCE_DIR}/CMake/run_pyscf.cmake")
set(LAST_TEST_NAME "NONE")
set(PYSCF_TEST_NAME "NONE")

################################################################################

run_pyscf_test(pyscf-diamond_1x1x1_pp_LCAO ${qmcpack_SOURCE_DIR}/examples/solids/pyscf-inputs/diamondC_1x1x1_pp diamondC_1x1x1_pp_LCAO
               PYSCF_TEST_NAME)

set_tests_properties(${PYSCF_TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "successfully saved to QMCPACK HDF5")

softlink_h5(${PYSCF_TEST_NAME} diamondC_1x1x1_pp-vmc_gaussian_sdj-1-16 C_Diamond C_Diamond.h5 LAST_TEST_NAME)

# LCAO and Bspline test should use the same references.
list(APPEND DIAMOND_SCALARS "totenergy" "-10.495941 0.0065")
list(APPEND DIAMOND_SCALARS "kinetic" "11.51198 0.078")
list(APPEND DIAMOND_SCALARS "potential" "-22.00792 0.079")
list(APPEND DIAMOND_SCALARS "eeenergy" "-2.68548 0.018")
list(APPEND DIAMOND_SCALARS "ionion" "-12.77566 0.0001")
list(APPEND DIAMOND_SCALARS "localecp" "-7.2011 0.090")
list(APPEND DIAMOND_SCALARS "nonlocalecp" "0.65437 0.024")
list(APPEND DIAMOND_SCALARS "samples" "16000 0.0")

set(THIS_TEST_NAME "${PYSCF_TEST_NAME}-diamondC_1x1x1_pp-vmc_gaussian_sdj")
qmc_run_and_check(
  ${THIS_TEST_NAME}
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
  qmc_short
  qmc_short.in.xml
  1
  16
  TRUE
  0
  DIAMOND_SCALARS # VMC
)

if(${THIS_TEST_NAME})
  set_tests_properties(${THIS_TEST_NAME}-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})
  set_property(
    TEST ${THIS_TEST_NAME}-1-16
    APPEND
    PROPERTY LABELS "converter")
endif()

################################################################################

run_pyscf_test(pyscf-diamond_1x1x1_pp_Bspline ${qmcpack_SOURCE_DIR}/examples/solids/pyscf-inputs/diamondC_1x1x1_pp
               diamondC_1x1x1_pp_Bspline PYSCF_TEST_NAME)

set_tests_properties(${PYSCF_TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "successfully saved to QMCPACK HDF5")

softlink_h5(${PYSCF_TEST_NAME} diamondC_1x1x1_pp-vmc_gaussian-bspline_sdj-1-16 C_Diamond C_Diamond-Bspline.h5
            LAST_TEST_NAME)

set(THIS_TEST_NAME "${PYSCF_TEST_NAME}-diamondC_1x1x1_pp-vmc_gaussian-bspline_sdj")
qmc_run_and_check(
  ${THIS_TEST_NAME}
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
  qmc_short-bspline
  qmc_short-bspline.in.xml
  1
  16
  TRUE
  0
  DIAMOND_SCALARS # VMC
)

if(${THIS_TEST_NAME})
  set_tests_properties(${THIS_TEST_NAME}-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})
  set_property(TEST ${THIS_TEST_NAME}-1-16 APPEND PROPERTY LABELS "converter")
endif()

################################################################################

# PYSCF_TEST_NAME set to pyscf-he2_1x1x1_sph2cart
# copy src dir to workdir
# add pyscf test: pyscf-he2_1x1x1_sph2cart
#    python he2_1x1x1_sph2cart.py
run_pyscf_test(pyscf-he2_1x1x1_sph ${qmcpack_SOURCE_DIR}/examples/solids/pyscf-inputs/he2_1x1x1_ae
	he2_1x1x1_sph2cart PYSCF_TEST_NAME)

set_property(TEST ${PYSCF_TEST_NAME} APPEND PROPERTY LABELS "s2c")

run_pyscf_s2c(${PYSCF_TEST_NAME} he_sp ${PYSCF_TEST_NAME}_sp_s2c)

softlink_h5(${PYSCF_TEST_NAME} he2_1x1x1_ae-vmc_sp_s2c_noj-1-16 he_sp_cart he_sp_cart.h5
	    LAST_TEST_NAME)

  list(APPEND HE2_SP_SCALARS "totenergy" "-2.32676861 0.0863")
  list(APPEND HE2_SP_SCALARS "kinetic" "16.23198428 0.1117")
  list(APPEND HE2_SP_SCALARS "potential" "-18.55875289 0.1136")
  list(APPEND HE2_SP_SCALARS "eeenergy" "4.59673345 0.0174")
  list(APPEND HE2_SP_SCALARS "ionion" "0.86420131 0.0001")
  list(APPEND HE2_SP_SCALARS "samples" "16000 0.0")

set(THIS_TEST_NAME "${PYSCF_TEST_NAME}-he2_1x1x1_ae-vmc_sp_s2c_noj")

qmc_run_and_check(
  ${THIS_TEST_NAME}
  "${qmcpack_SOURCE_DIR}/tests/solids/he2_1x1x1_ae"
  he2_sp_cart
  he2_sp_cart.qmc.in-wfnoj.xml
  1
  16
  TRUE
  0
  HE2_SP_SCALARS # VMC
  )

if(${THIS_TEST_NAME})
  set_tests_properties(${THIS_TEST_NAME}-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})
  set_property(TEST ${THIS_TEST_NAME}-1-16 APPEND PROPERTY LABELS "converter")
endif()

run_pyscf_s2c(${PYSCF_TEST_NAME} he_sd ${PYSCF_TEST_NAME}_sd_s2c)

softlink_h5(${PYSCF_TEST_NAME} he2_1x1x1_ae-vmc_sd_s2c_noj-1-16 he_sd_cart he_sd_cart.h5
            LAST_TEST_NAME)

  list(APPEND HE2_SD_SCALARS "totenergy" "-2.32942007 0.0835")
  list(APPEND HE2_SD_SCALARS "kinetic" "16.21547981 0.1067")
  list(APPEND HE2_SD_SCALARS "potential" "-18.54489989 0.1189")
  list(APPEND HE2_SD_SCALARS "eeenergy" "4.59377532 0.0184")
  list(APPEND HE2_SD_SCALARS "ionion" "0.86420131 0.0001")
  list(APPEND HE2_SD_SCALARS "samples" "16000 0.0")

set(THIS_TEST_NAME "${PYSCF_TEST_NAME}-he2_1x1x1_ae-vmc_sd_s2c_noj")

qmc_run_and_check(
  ${THIS_TEST_NAME}
  "${qmcpack_SOURCE_DIR}/tests/solids/he2_1x1x1_ae"
  he2_sd_cart
  he2_sd_cart.qmc.in-wfnoj.xml
  1
  16
  TRUE
  0
  HE2_SD_SCALARS # VMC
  )

if(${THIS_TEST_NAME})
  set_tests_properties(${THIS_TEST_NAME}-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})
  set_property(TEST ${THIS_TEST_NAME}-1-16 APPEND PROPERTY LABELS "converter")
endif()

run_pyscf_s2c(${PYSCF_TEST_NAME} he_sf ${PYSCF_TEST_NAME}_sf_s2c)

softlink_h5(${PYSCF_TEST_NAME} he2_1x1x1_ae-vmc_sf_s2c_noj-1-16 he_sf_cart he_sf_cart.h5
            LAST_TEST_NAME)

  list(APPEND HE2_SF_SCALARS "totenergy" "-2.35068897 0.1249")
  list(APPEND HE2_SF_SCALARS "kinetic" "16.20083403 0.1568")
  list(APPEND HE2_SF_SCALARS "potential" "-18.55152300 0.1705")
  list(APPEND HE2_SF_SCALARS "eeenergy" "4.59192604 0.0261")
  list(APPEND HE2_SF_SCALARS "ionion" "0.86420131 0.0001")
  list(APPEND HE2_SF_SCALARS "samples" "16000 0.0")

set(THIS_TEST_NAME "${PYSCF_TEST_NAME}-he2_1x1x1_ae-vmc_sf_s2c_noj")

qmc_run_and_check(
  ${THIS_TEST_NAME}
  "${qmcpack_SOURCE_DIR}/tests/solids/he2_1x1x1_ae"
  he2_sf_cart
  he2_sf_cart.qmc.in-wfnoj.xml
  1
  16
  TRUE
  0
  HE2_SF_SCALARS # VMC
  )

if(${THIS_TEST_NAME})
  set_tests_properties(${THIS_TEST_NAME}-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})
  set_property(TEST ${THIS_TEST_NAME}-1-16 APPEND PROPERTY LABELS "converter")
endif()

run_pyscf_s2c(${PYSCF_TEST_NAME} he_sg ${PYSCF_TEST_NAME}_sg_s2c)

softlink_h5(${PYSCF_TEST_NAME} he2_1x1x1_ae-vmc_sg_s2c_noj-1-16 he_sg_cart he_sg_cart.h5
            LAST_TEST_NAME)

  list(APPEND HE2_SG_SCALARS "totenergy" "-2.30636188 0.1233")
  list(APPEND HE2_SG_SCALARS "kinetic" "16.21415978 0.1493")
  list(APPEND HE2_SG_SCALARS "potential" "-18.52052166 0.1655")
  list(APPEND HE2_SG_SCALARS "eeenergy" "4.59711595 0.0251")
  list(APPEND HE2_SG_SCALARS "ionion" "0.86420131 0.0001")
  list(APPEND HE2_SG_SCALARS "samples" "16000 0.0")

  set(THIS_TEST_NAME "${PYSCF_TEST_NAME}-he2_1x1x1_ae-vmc_sg_s2c_noj")

qmc_run_and_check(
  ${THIS_TEST_NAME}
  "${qmcpack_SOURCE_DIR}/tests/solids/he2_1x1x1_ae"
  he2_sg_cart
  he2_sg_cart.qmc.in-wfnoj.xml
  1
  16
  TRUE
  0
  HE2_SG_SCALARS # VMC
  )

if(${THIS_TEST_NAME})
  set_tests_properties(${THIS_TEST_NAME}-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})
  set_property(TEST ${THIS_TEST_NAME}-1-16 APPEND PROPERTY LABELS "converter")
endif()

run_pyscf_s2c(${PYSCF_TEST_NAME} he_sh ${PYSCF_TEST_NAME}_sh_s2c)

softlink_h5(${PYSCF_TEST_NAME} he2_1x1x1_ae-vmc_sh_s2c_noj-1-16 he_sh_cart he_sh_cart.h5
            LAST_TEST_NAME)

  list(APPEND HE2_SH_SCALARS "totenergy" "-2.31529851 0.0837")
  list(APPEND HE2_SH_SCALARS "kinetic" "16.21254352 0.1095")
  list(APPEND HE2_SH_SCALARS "potential" "-18.52784203 0.1163")
  list(APPEND HE2_SH_SCALARS "eeenergy" "4.59749106 0.0174")
  list(APPEND HE2_SH_SCALARS "ionion" " 0.86420131 0.0001")
  list(APPEND HE2_SH_SCALARS "samples" "16000 0.0")

  set(THIS_TEST_NAME "${PYSCF_TEST_NAME}-he2_1x1x1_ae-vmc_sh_s2c_noj")

qmc_run_and_check(
  ${THIS_TEST_NAME}
  "${qmcpack_SOURCE_DIR}/tests/solids/he2_1x1x1_ae"
  he2_sh_cart
  he2_sh_cart.qmc.in-wfnoj.xml
  1 
  16
  TRUE
  0
  HE2_SH_SCALARS # VMC
  )

if(${THIS_TEST_NAME})
  set_tests_properties(${THIS_TEST_NAME}-1-16 PROPERTIES DEPENDS ${LAST_TEST_NAME})
  set_property(TEST ${THIS_TEST_NAME}-1-16 APPEND PROPERTY LABELS "converter")
endif()
