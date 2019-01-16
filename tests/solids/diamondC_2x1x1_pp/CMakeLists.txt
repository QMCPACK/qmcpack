
  LIST(APPEND DIAMOND2_SCALARS "totenergy" "-21.668173343662 0.007001625578")
  LIST(APPEND DIAMOND2_SCALARS "kinetic" "20.065101661866 0.043673620402")
  LIST(APPEND DIAMOND2_SCALARS "potential" "-41.733275005528 0.044730744420")
  LIST(APPEND DIAMOND2_SCALARS "eeenergy" "-5.428882678829 0.011970082813")
  LIST(APPEND DIAMOND2_SCALARS "ionion" "-25.551327183000 0.001")
  LIST(APPEND DIAMOND2_SCALARS "localecp" "-12.293639902627 0.051239542940")
  LIST(APPEND DIAMOND2_SCALARS "nonlocalecp" "1.540574758843 0.018460717355")
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

  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp-delayed_update-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_short
                    qmc_short_delayedupdate.in.xml
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
  LIST(APPEND DIAMOND2_DMC_SCALARS "totenergy" "-21.846906711853 0.011365387412")
  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1_pp"
                    qmc_short_vmc_dmc
                    qmc_short_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 DIAMOND2_DMC_SCALARS # DMC
                    )

  LIST(APPEND DIAMOND2_DMC_SCALARS "totenergy" "-21.846906711853 0.011365387412")
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

  LIST(APPEND LONG_DIAMOND2_SCALARS "totenergy" "-21.668190934054 0.000239004082")
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
  LIST(APPEND LONG_DIAMOND2_DMC_SCALARS "totenergy" "-21.846906711853 0.003503129238")
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



