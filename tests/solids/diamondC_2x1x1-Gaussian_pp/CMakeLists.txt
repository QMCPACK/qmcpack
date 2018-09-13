

# Reference VMC run in qmc-ref "-21.70152 +/- 0.00042"
IF(ENABLE_SOA)
  LIST(APPEND DIAMOND_SCALARS "totenergy" "-21.70152 0.015")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "20.6385 0.16 ")
  LIST(APPEND DIAMOND_SCALARS "potential" "-42.3400 0.19")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-5.33913 0.034")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-25.551326 0.0001")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-13.1852 0.20")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "1.7356 0.05")
  LIST(APPEND DIAMOND_SCALARS "samples" "8000 0.0")
#  LIST(APPEND DIAMOND_SCALARS "flux" "0.0 0.4")

  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp-vmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
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

  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp-vmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )


# Reference DMC run in qmc-ref "-21.84712 +/-  0.00068"
  LIST(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-21.84712 0.009024")
  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp-dmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 LONG_DIAMOND_DMC_SCALARS # DMC
                    )

ELSE()
  MESSAGE("Skipping all tests since PBC gaussians only enabled with SOA")
ENDIF()
