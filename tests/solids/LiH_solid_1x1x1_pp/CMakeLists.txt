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

IF(ENABLE_SOA OR QMC_CUDA)
  QMC_RUN_AND_CHECK(short-LiH_solid_1x1x1_hybridrep_pp-x-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    hf_vmc_LiH-x
                    hf_vmc_LiH-x_hybridrep.xml
                    1 16
                    TRUE
                    0 LIH_X_SCALARS # VMC
                    )
ELSE()
  MESSAGE("Skipping LiH_solid_1x1x1_hybridrep_pp-x-vmc_hf_noj test because neither SOA nor GPU is enabled")
ENDIF()

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
  LIST(APPEND DET_LIH_GAMMA_SCALARS "totenergy" "-8.30870557 0.00002841")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "kinetic" "0.65348406 0.00002835")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "potential" "-8.96218962 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "eeenergy" "-0.92293001 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "ionion" "-3.68922767 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "localecp" "-4.35003194 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "samples" "9 0.0")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "flux" "12.49117099 0.00005584")
 ELSE()
  LIST(APPEND DET_LIH_GAMMA_SCALARS "totenergy" "-7.97064375 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "kinetic" "3.17905650 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "potential" "-11.14970025 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "eeenergy" "-0.07789848 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "ionion" "-3.68922683 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "localecp" "-7.38257493 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "samples" "9 0.0")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "flux" "9.26788839 0.000001")
 ENDIF()
ELSE()
 IF(QMC_MIXED_PRECISION)
  LIST(APPEND DET_LIH_GAMMA_SCALARS "totenergy" "-8.62473375 0.0000118")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "kinetic" "18.17890605 0.00001662")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "potential" "-26.80363981 0.00000483")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "eeenergy" "0.27617574 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "ionion" "-3.68922683 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "localecp" "-23.39058871 0.00000489")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "samples" "9 0.0")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "flux" "-26.14287483 0.00004006")
 ELSE()
  LIST(APPEND DET_LIH_GAMMA_SCALARS "totenergy" "-8.62476955 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "kinetic" "18.17889660 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "potential" "-26.80366615 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "eeenergy" "0.27617595 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "ionion" "-3.68922682 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "localecp" "-23.39061528 0.000001")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "samples" "9 0.0")
  LIST(APPEND DET_LIH_GAMMA_SCALARS "flux" "-26.14286403 0.000001")
 ENDIF()
ENDIF()

  QMC_RUN_AND_CHECK(deterministic-short-LiH_solid_1x1x1_pp-gamma-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/solids/LiH_solid_1x1x1_pp"
                    det_hf_vmc_LiH-gamma
                    det_hf_vmc_LiH-gamma.xml
                    1 1
                    TRUE
                    0 DET_LIH_GAMMA_SCALARS # VMC
                    )
