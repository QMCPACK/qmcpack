
IF(ENABLE_SOA)

# Reference VMC run in qmc-ref "-10.47755 0.00073"
  LIST(APPEND DIAMOND_SCALARS "totenergy" "-10.47858 0.00094")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "11.5576 0.0065")
  LIST(APPEND DIAMOND_SCALARS "potential" "-22.0362 0.0063")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-2.6766 0.0015")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-12.79 0.00")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-7.2636 0.0080")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "0.6910 0.0021")
  LIST(APPEND DIAMOND_SCALARS "samples" "1536000 0")



  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_crystal_sdj
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Crystal_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )



#
# Long tests
#

  LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-10.47834 0.00025")
  LIST(APPEND LONG_DIAMOND_SCALARS "samples" "15360000 0")
  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc_crystal_sdj
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Crystal_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )

# Reference DMC run in qmc-ref "-10.52158 0.00036"
  LIST(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-10.5176 0.0071")

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-dmc_Crystal_sdj
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Crystal_pp"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 LONG_DIAMOND_DMC_SCALARS # DMC
                    )
ELSE()
  MESSAGE("Skipping all tests since PBC gaussians only enabled with SOA")
ENDIF()

