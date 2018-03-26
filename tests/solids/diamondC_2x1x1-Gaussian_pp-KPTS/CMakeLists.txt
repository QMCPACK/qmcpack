

IF(ENABLE_SOA)
  LIST(APPEND DIAMOND_SCALARS "totenergy" "-21.702829 0.017")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "20.706 0.067")
  LIST(APPEND DIAMOND_SCALARS "potential" "-42.409 0.067")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-5.361 0.021")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-25.55 0.00")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-13.244 0.079")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "1.748 0.045")
  LIST(APPEND DIAMOND_SCALARS "samples" "6880 0.0")
#  LIST(APPEND DIAMOND_SCALARS "flux" "0.0 0.4")

  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp-vmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )


# Reference DMC run in qmc-ref "-21.807657 0.095185"
  LIST(APPEND DIAMOND_DMC_SCALARS "totenergy" "-21.807657 0.095185")
  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp-dmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
                    qmc_short_vmc_dmc
                    qmc_short_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 DIAMOND_DMC_SCALARS # DMC
                    )


#
# Long tests
#

  LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-21.709193  0.006340")
  LIST(APPEND LONG_DIAMOND_SCALARS "samples" "49600 0.0")

  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp-vmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )


# Reference DMC run in qmc-ref "-21.840087 +/- 0.026449"
  LIST(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-21.840087 0.026449")
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
