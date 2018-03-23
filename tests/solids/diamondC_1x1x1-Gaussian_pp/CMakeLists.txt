


IF(ENABLE_SOA)

  LIST(APPEND DIAMOND_SCALARS "totenergy" "-10.4957 0.0061")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "11.517 0.055")
  LIST(APPEND DIAMOND_SCALARS "potential" "-22.013 0.054")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-2.679 0.013")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-12.78 0.00")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-7.196 0.060")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "0.638 0.025")
  LIST(APPEND DIAMOND_SCALARS "samples" "128000 0.0")
#  LIST(APPEND DIAMOND_SCALARS "flux" "0.0 0.4")

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_gausian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )


# Reference DMC run in qmc-ref "-10.541 0.031"
  LIST(APPEND DIAMOND_DMC_SCALARS "totenergy" "-10.541 0.031")
  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-dmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
                    qmc_short_vmc_dmc
                    qmc_short_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 DIAMOND_DMC_SCALARS # DMC
                    )


#
# Long tests
#

  LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-10.4951 0.0017")
  LIST(APPEND LONG_DIAMOND_SCALARS "samples" "192640 0.0")

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )


# Reference DMC run in qmc-ref "-10.544774 0.014660"
  LIST(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-10.544774 0.014408")
  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-dmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 LONG_DIAMOND_DMC_SCALARS # DMC
                    )

ELSE()
  MESSAGE("Skipping all tests since PBC gaussians only enabled with SOA")
ENDIF()
