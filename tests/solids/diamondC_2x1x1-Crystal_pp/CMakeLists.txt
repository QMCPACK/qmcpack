

# Reference VMC run in qmc-ref ""
IF(ENABLE_SOA)
  LIST(APPEND DIAMOND_SCALARS "totenergy" "-21.6834 0.0038")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "20.658 0.019")
  LIST(APPEND DIAMOND_SCALARS "potential" "-42.341 0.018")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-5.3405 0.0042")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-25.57 0.0")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-13.249  0.018")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "1.8217 0.0068")
  LIST(APPEND DIAMOND_SCALARS "samples" "245760 0")

  QMC_RUN_AND_CHECK(short-diamondC_2x1x1_pp-vmc_crystal_sdj
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Crystal_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )

#
# Long tests
#

# Reference VMC run in qmc-ref "
  LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-21.6773 0.0014")
  LIST(APPEND LONG_DIAMOND_SCALARS "samples" "1029120 0")

  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp-vmc_crystal_sdj
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Crystal_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )


# Reference DMC run in qmc-ref "
  LIST(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-21.815 0.019")
  QMC_RUN_AND_CHECK(long-diamondC_2x1x1_pp-dmc_crystal_sdj
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Crystal_pp"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 LONG_DIAMOND_DMC_SCALARS # DMC
                    )

ELSE()
  MESSAGE("Skipping all tests since PBC gaussians only enabled with SOA")
ENDIF()
