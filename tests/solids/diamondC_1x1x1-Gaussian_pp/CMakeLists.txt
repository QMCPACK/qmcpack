
IF(ENABLE_SOA)

# Reference VMC run in qmc-ref "-10.495941  0.000073"
  LIST(APPEND DIAMOND_SCALARS "totenergy" "-10.495853296208 0.001704169726")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "11.510223929044 0.015701736226")
  LIST(APPEND DIAMOND_SCALARS "potential" "-22.006077225239 0.015711745694")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-2.685859088397 0.003828126585")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-12.775667436000 0.0001")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-7.198842031405 0.017987159139")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "0.654291330991 0.005204188163")
  LIST(APPEND DIAMOND_SCALARS "samples" "16000 0.0")



  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )



#
# Long tests
#

  LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-10.495853296208 0.000544919799")
  LIST(APPEND LONG_DIAMOND_SCALARS "samples" "160000 0.0")
  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc_gaussian_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )

# Reference DMC run in qmc-ref "-10.532497 0.000470"
  LIST(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-10.532497 0.006127")

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

