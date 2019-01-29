
# Reference VMC run in qmc-ref "-21.71006 +/- 0.00057"
IF(ENABLE_SOA)
  LIST(APPEND DIAMOND_SCALARS "totenergy" "-21.71006 0.018")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "20.4859  0.089")
  LIST(APPEND DIAMOND_SCALARS "potential" "-42.1959  0.088")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-5.3556 0.020")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-25.55132676 0.0001")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-13.1306 0.13")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "1.8307 0.032")
  LIST(APPEND DIAMOND_SCALARS "samples" "8000 0.0")
#  LIST(APPEND DIAMOND_SCALARS "flux" "0.0 0.4")

  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp_kpts-vmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_kpts"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )

#
# Long tests
#

# Reference VMC run in qmc-ref "-21.70152 +/- 0.00042"
  LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-21.70152   0.008351")
  LIST(APPEND LONG_DIAMOND_SCALARS "samples" "32000 0.0")

  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp_kpts-vmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_kpts"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )


# Reference DMC run in qmc-ref "-21.84425 +/-  0.00083"
  LIST(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-21.84425 0.009024")
  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp_kpts-dmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_kpts"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 LONG_DIAMOND_DMC_SCALARS # DMC
                    )

ELSE()
  MESSAGE("Skipping all tests since PBC gaussians only enabled with SOA")
ENDIF()
