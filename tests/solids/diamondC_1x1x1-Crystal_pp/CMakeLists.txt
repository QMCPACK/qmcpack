
IF(ENABLE_SOA)

# Reference VMC run in qmc-ref "-10.47816 0.00016"
  LIST(APPEND DIAMOND_SCALARS "totenergy" "-10.47792 0.00048")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "11.5675 0.0047")
  LIST(APPEND DIAMOND_SCALARS "potential" "-22.0455 0.0046")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-2.6724 0.0012")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-12.79 0.00")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-7.2774 0.0048")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "0.6914 0.0012")
  LIST(APPEND DIAMOND_SCALARS "samples" "3264000 0")



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

  LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-10.47816 0.00016")
  LIST(APPEND LONG_DIAMOND_SCALARS "samples" "32640000 0")
  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc_crystal_sdj
      "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Crystal_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )

# Reference DMC run in qmc-ref "-10.532497 0.000470"
  LIST(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-10.5290 0.0039")

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-dmc_crystal_sdj
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

