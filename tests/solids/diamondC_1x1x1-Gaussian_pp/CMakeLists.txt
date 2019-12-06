
IF(ENABLE_SOA)

# Reference VMC run in qmc-ref "-10.495941  0.000073"
  LIST(APPEND DIAMOND_SCALARS "totenergy" "-10.495941 0.0065")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "11.51198 0.078")
  LIST(APPEND DIAMOND_SCALARS "potential" "-22.00792 0.079")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-2.68548 0.018")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-12.77566 0.0001")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-7.2011 0.090")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "0.65437 0.024")
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

  LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-10.495941 0.0017")
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

