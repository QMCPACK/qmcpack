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

  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_hybridrep_pp-x-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-x
                    hf_vmc_LiH-x_hybridrep.xml
                    1 16
                    TRUE
                    0 LIH_X_SCALARS # VMC
                    )

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

# multi-slater determinant requires complex build to run. Ensure fails with real/non-complex build
  LIST(APPEND LIH_MD_GAMMA_SCALARS "totenergy" "-8.489858484145001 0.004")
  IF (QMC_COMPLEX)
    QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-gamma-vmc_msdj
                      "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                      md_cc_vmc_LiH-gamma
                      md_cc_vmc_LiH-gamma.xml
                      1 16
                      TRUE
                      1 LIH_MD_GAMMA_SCALARS # VMC
                      )

  ELSE()
    QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_pp-gamma-vmc_msdj
                      "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                      md_rc_vmc_LiH-gamma
                      md_rc_vmc_LiH-gamma.xml
                      1 16
                      TRUE                  
                      1 LIH_MD_GAMMA_SCALARS # VMC
                      )
  ENDIF()


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

# Deterministic test
IF(NOT QMC_CUDA)
 IF(QMC_MIXED_PRECISION)
  LIST(APPEND DET_LIH_GAMMA_SCALARS "totenergy" "-8.99301725 0.00001597")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "kinetic" "7.42550431 0.00001588")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "potential" "-16.41852156 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "eeenergy" "-1.39002015 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "ionion" "-3.68928007 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "localecp" "-11.33922134 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "samples" "9 0.0")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "flux" "0.04160754 0.00002925")
 ELSE()
  LIST(APPEND DET_LIH_GAMMA_SCALARS "totenergy" "-8.99302939 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "kinetic" "7.42549769 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "potential" "-16.41852708 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "eeenergy" "-1.39001949 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "ionion" "-3.68927959 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "localecp" "-11.33922800 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "samples" "9 0.0")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "flux" "0.04162024 0.000001")
 ENDIF()
ELSE()
 IF(QMC_MIXED_PRECISION)
  LIST(APPEND DET_LIH_GAMMA_SCALARS "totenergy" "-8.40171773 0.00001393")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "kinetic" "0.28091468 0.00001482")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "potential" "-8.68263241 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "eeenergy" "-1.40882665 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "ionion" "-3.68927951 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "localecp" "-3.58452626 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "samples" "9 0.0")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "flux" "13.05705762 0.00002873")
 ELSE()
  LIST(APPEND DET_LIH_GAMMA_SCALARS "totenergy" "-8.40171568 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "kinetic" "0.28091457 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "potential" "-8.68263024 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "eeenergy" "-1.40882677 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "ionion" "-3.68927951 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "localecp" "-3.58452396 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "samples" "9 0.0")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "flux" "13.05706070 0.000001")
 ENDIF()
ENDIF()

 QMC_RUN_AND_CHECK(deterministic-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj
                   "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                   det_hf_vmc_LiH-gamma
                   det_hf_vmc_LiH-gamma.xml
                   1 1
                   TRUE
                   0 DET_LIH_GAMMA_SCALARS # VMC
                   )
