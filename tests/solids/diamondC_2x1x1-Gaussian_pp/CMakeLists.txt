

# Reference VMC run in qmc-ref "-21.70152 +/- 0.00042"
IF(ENABLE_SOA)
  LIST(APPEND DIAMOND_SCALARS "totenergy" "-21.701059413048 0.004507698620")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "20.639674255315 0.030277362808 ")
  LIST(APPEND DIAMOND_SCALARS "potential" "-42.340733668372 0.030384241216")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-5.337307250298 0.007597956495")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-25.551326937000 0.0001")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-13.192050121857 0.034990790950")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "1.739950640646 0.011286094181")
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
  LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-21.701059413048   0.002239292047")
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
