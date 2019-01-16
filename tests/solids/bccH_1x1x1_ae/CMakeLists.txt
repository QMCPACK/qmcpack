

#
# Short tests, about 2-3 seconds on 16 cores
#
  LIST(APPEND BCC_H_SCALARS "totenergy" "-1.834030571913 0.000127343381")
  LIST(APPEND BCC_H_SCALARS "kinetic" "0.186549257009 0.000898028235")
  LIST(APPEND BCC_H_SCALARS "potential" "-2.020579828921 0.000862358211")
  LIST(APPEND BCC_H_SCALARS "eeenergy" "-0.775851098291 0.000122455703")
  LIST(APPEND BCC_H_SCALARS "ionion" "-0.962899625270 0.001")
  LIST(APPEND BCC_H_SCALARS "samples" "960000 0.0")
#  LIST(APPEND BCC_H_SCALARS "flux" "0.0 0.03")

  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 BCC_H_SCALARS # VMC
                    )


# qmc_ref_vmc_dmc  series 1  -1.841352 +/- 0.000078
# qmc_short_vmc_dmc  series 1  -1.84196 +/- 0.00073
  LIST(APPEND BCC_H_DMC_SCALARS "totenergy" "-1.841388542839 0.000525334459")
  LIST(APPEND BCC_H_DMC_SCALARS "potential" "-2.020107240639  0.003261277449")
  LIST(APPEND BCC_H_DMC_SCALARS "kinetic"    "0.178718697799  0.003328529697")
  LIST(APPEND BCC_H_DMC_SCALARS "eeenergy"  "-0.779885842929 0.000787672689")
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_vmc_dmc
                    qmc_short_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 BCC_H_DMC_SCALARS # DMC
                    )

IF (NOT QMC_CUDA)

# Check VMCUpdateAll driver using all-particle updates
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-vmc-all-nodrift_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_vmc_all_nodrift
                    qmc_short_vmc_all_nodrift.in.xml
                    1 16
                    TRUE
                    0 BCC_H_SCALARS # VMC
                    )

# Check VMCUpdateAllWithDrift driver using all-particle updates
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-vmc-all_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_vmc_all
                    qmc_short_vmc_all.in.xml
                    1 16
                    TRUE
                    0 BCC_H_SCALARS # VMC
                    )

# Check DMCUpdateAllWithRejection driver using all-particle updates
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-dmc-all_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_vmc_dmc_all
                    qmc_short_vmc_dmc_all.in.xml
                    1 16
                    TRUE
                    1 BCC_H_DMC_SCALARS # DMC
                    )

#Now for the correlated sampling tests
  LIST(APPEND BCC_H_CSVMC_SCALARS "totenergy_A" "-1.834030973656 0.000128283895")
  LIST(APPEND BCC_H_CSVMC_SCALARS "totenergy_B" "-1.833774860468 0.000128408851")
  LIST(APPEND BCC_H_CSVMC_SCALARS "dtotenergy_AB" "-0.000256113188 0.000071573900")
  LIST(APPEND BCC_H_CSVMC_SCALARS "ionion_A" "-0.962899621840 0.0000001")
  LIST(APPEND BCC_H_CSVMC_SCALARS "ionion_B" "-0.958165599910 0.0000001")
  LIST(APPEND BCC_H_CSVMC_SCALARS "dionion_AB" "-0.004734021929 0.0000001")
  LIST(APPEND BCC_H_CSVMC_SCALARS "eeenergy_A" "-0.775840944343 0.000122329486")
  LIST(APPEND BCC_H_CSVMC_SCALARS "eeenergy_B" "-0.775180796755 0.000123008833")
  LIST(APPEND BCC_H_CSVMC_SCALARS "deeenergy_AB" "-0.000660147588 0.000036878212")

# Check CSVMCUpdatePbyP driver using particle-by-particle updates
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-csvmc-pbyp-nodrift_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_csvmc_pbyp_nodrift
                    qmc_short_csvmc_pbyp_nodrift.in.xml
                    16 1
                    TRUE
                    0 BCC_H_CSVMC_SCALARS # VMC
                    )

# Check CSVMCUpdatePbyPWithDrift driver using all-particle updates
#  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-csvmc-pbyp_sdj
#                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
#                    qmc_short_csvmc_pbyp
#                    qmc_short_csvmc_pbyp.in.xml
#                    16 1
#                    TRUE
#                    0 BCC_H_CSVMC_SCALARS # VMC
#                    )

# Check CSVMCUpdateAll driver using all-particle updates
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-csvmc-all-nodrift_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_csvmc_all_nodrift
                    qmc_short_csvmc_all_nodrift.in.xml
                    16 1
                    TRUE
                    0 BCC_H_CSVMC_SCALARS # VMC
                    )

# Check CSVMCUpdateAllWithDrift driver using all-particle updates
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-csvmc-all_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_csvmc_all
                    qmc_short_csvmc_all.in.xml
                    16 1
                    TRUE
                    0 BCC_H_CSVMC_SCALARS # VMC
                    )


  LIST(APPEND BCC_H_RMC_SCALARS "totenergy" "-1.8420 0.0011")
  #check RMC on 16 cores.  1 MPI
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-rmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_vmc_rmc
                    qmc_short_vmc_rmc.in.xml
                    1 16
                    TRUE
                    1 BCC_H_RMC_SCALARS #series# 
                    )
  #check RMC on 4 cores, 4 mpi
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-rmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_vmc_rmc
                    qmc_short_vmc_rmc.in.xml
                    4 4
                    TRUE
                    1 BCC_H_RMC_SCALARS #series#
                    )

  # Check CPU time limit
  CPU_LIMIT_RUN("bccH_1x1x1_ae-rmc-1-1"
              "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
              1 1
              120
              qmc_cpu_limit_rmc.xml)

ELSE()
  MESSAGE("Skipping bccH_1x1x1_ae tests for VMC all-particle moves, CSVMC, and RMC because they are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

# Plane wave test of bccH (no pp, so faster than typical solids)
# Reduced number of blocks and samples compared to other short bccH tests
  LIST(APPEND BCC_H_PW_SCALARS "totenergy" "-1.834006336519 0.000225118131")
  LIST(APPEND BCC_H_PW_SCALARS "samples" "96000 0.0")
#  LIST(APPEND BCC_H_PW_SCALARS "flux" "0.0 0.1")
# No GPU implementation. Would be fast if implemented
IF(NOT QMC_CUDA)
  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-vmc_sdj-pw
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_short_pw
                    qmc_short_pw.in.xml
                    1 16
                    TRUE
                    0 BCC_H_PW_SCALARS # VMC
                    )
ELSE()
  MESSAGE("Skipping bccH_1x1x1_ae plane wave tests because they are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()


#
# Long variants of above, about 1/2 hour with 16 cores
#
  LIST(APPEND LONG_BCC_H_SCALARS "totenergy" "-1.834033399740 0.000004159080")
  LIST(APPEND LONG_BCC_H_SCALARS "samples" "921600000 0.0")

  QMC_RUN_AND_CHECK(long-bccH_1x1x1_ae-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_BCC_H_SCALARS # VMC
                    )


# qmc_ref_vmc_dmc  series 1  -1.841352 +/- 0.000078
# qmc_long_vmc_dmc  series 1   -1.841204 +/- 0.000285
  LIST(APPEND LONG_BCC_H_DMC_SCALARS "totenergy" "-1.841388542839 0.000233993188")
  QMC_RUN_AND_CHECK(long-bccH_1x1x1_ae-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 LONG_BCC_H_DMC_SCALARS # DMC
                    )

IF(NOT QMC_CUDA)
  LIST(APPEND LONG_BCC_H_RMC_SCALARS "totenergy" "-1.84150 0.00011")
  LIST(APPEND LONG_BCC_H_RMC_SCALARS "potential" "-2.01964 0.00072")
  LIST(APPEND LONG_BCC_H_RMC_SCALARS "potential_pure" "-2.01907 0.00085")
  LIST(APPEND LONG_BCC_H_RMC_SCALARS "kinetic_mixed" "0.17814 0.00073")
  LIST(APPEND LONG_BCC_H_RMC_SCALARS "kinetic_pure" "0.17063 0.00086")
  LIST(APPEND LONG_BCC_H_RMC_SCALARS "eeenergy_mixed" "-0.77991 0.00016")
  LIST(APPEND LONG_BCC_H_RMC_SCALARS "eeenergy_pure" "-0.78335 0.00019")
  QMC_RUN_AND_CHECK(long-bccH_1x1x1_ae-rmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
                    qmc_long_vmc_rmc
                    qmc_long_vmc_rmc.in.xml
                    1 16
                    TRUE
                    1 LONG_BCC_H_RMC_SCALARS # DMC
                    )
ENDIF()

