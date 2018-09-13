
  LIST(APPEND DIAMOND2_SCALARS "totenergy" "-21.668293 0.007")
  LIST(APPEND DIAMOND2_SCALARS "kinetic" "20.065151 0.044")
  LIST(APPEND DIAMOND2_SCALARS "potential" "-41.733377 0.045")
  LIST(APPEND DIAMOND2_SCALARS "eeenergy" "-5.428561 0.012")
  LIST(APPEND DIAMOND2_SCALARS "ionion" "-25.5513269 0.001")
  LIST(APPEND DIAMOND2_SCALARS "localecp" "-12.29402 0.051")
  LIST(APPEND DIAMOND2_SCALARS "nonlocalecp" "1.540530 0.018")
  LIST(APPEND DIAMOND2_SCALARS "samples" "32000 0.0")
#  LIST(APPEND DIAMOND2_SCALARS "flux" "0.0 0.7")

  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 DIAMOND2_SCALARS # VMC
                    )

IF(ENABLE_SOA OR QMC_CUDA)
  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_hybridrep_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_short
                    qmc_short_hybridrep.in.xml
                    1 16
                    TRUE
                    0 DIAMOND2_SCALARS # VMC
                    )
ELSE()
  MESSAGE("Skipping diamondC_2x1x1_hybridrep_pp-vmc_sdj test because neither SOA nor CUDA is enabled")
ENDIF()

# Reference DMC run in qmc-ref "-21.844975 0.001113"
  LIST(APPEND DIAMOND2_DMC_SCALARS "totenergy" "-21.844975 0.00779")
  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_short_vmc_dmc
                    qmc_short_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 DIAMOND2_DMC_SCALARS # DMC
                    )

  LIST(APPEND DIAMOND2_DMC_SCALARS "totenergy" "-21.844975 0.00779")
  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp-dmc-reconf_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_short_vmc_dmc
                    qmc_short_vmc_dmc_reconf.in.xml
                    1 16
                    TRUE
                    1 DIAMOND2_DMC_SCALARS # DMC
                    )


#
# Long tests
#

  LIST(APPEND LONG_DIAMOND2_SCALARS "totenergy" "-21.668293 0.00023")
  LIST(APPEND LONG_DIAMOND2_SCALARS "samples" "30720000 0.0")
  LIST(APPEND LONG_DIAMOND2_SCALARS "flux" "0.0 0.05")

  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND2_SCALARS # VMC
                    )

# Reference DMC run in qmc-ref "-21.844975 0.001113"
  LIST(APPEND LONG_DIAMOND2_DMC_SCALARS "totenergy" "-21.844975 0.002879")
  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 LONG_DIAMOND2_DMC_SCALARS # DMC
                    )

  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp-dmc-reconf_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc_reconf.in.xml
                    1 16
                    TRUE
                    1 LONG_DIAMOND2_DMC_SCALARS # DMC
                    )



