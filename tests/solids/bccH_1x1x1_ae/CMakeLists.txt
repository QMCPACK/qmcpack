

#
# Short tests, about 2-3 seconds on 16 cores
#
  LIST(APPEND BCC_H_SCALARS "totenergy" "-1.834032 0.0003")
  LIST(APPEND BCC_H_SCALARS "kinetic" "0.186673807952 0.003")
  LIST(APPEND BCC_H_SCALARS "potential" "-2.02070672998 0.003")
  LIST(APPEND BCC_H_SCALARS "eeenergy" "-0.775870180954 0.003")
  LIST(APPEND BCC_H_SCALARS "ionion" "-0.96289961022 0.001")
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
  LIST(APPEND BCC_H_DMC_SCALARS "totenergy" "-1.84196 0.00073")
  LIST(APPEND BCC_H_DMC_SCALARS "potential" "-2.0187  0.0033")
  LIST(APPEND BCC_H_DMC_SCALARS "kinetic"    "0.1767  0.0035")
  LIST(APPEND BCC_H_DMC_SCALARS "eeenergy"  "-0.77985 0.00069")
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
  LIST(APPEND BCC_H_CSVMC_SCALARS "totenergy_A" "-1.83428 0.00013")
  LIST(APPEND BCC_H_CSVMC_SCALARS "totenergy_B" "-1.83390 0.00012")
  LIST(APPEND BCC_H_CSVMC_SCALARS "dtotenergy_AB" "-0.000384 0.000068")
  LIST(APPEND BCC_H_CSVMC_SCALARS "ionion_A" "-0.96289962 0.0000001")
  LIST(APPEND BCC_H_CSVMC_SCALARS "ionion_B" "-0.95816559 0.0000001")
  LIST(APPEND BCC_H_CSVMC_SCALARS "dionion_AB" "-0.004734021 0.0000001")
  LIST(APPEND BCC_H_CSVMC_SCALARS "eeenergy_A" "-0.77607 0.00013")
  LIST(APPEND BCC_H_CSVMC_SCALARS "eeenergy_B" "-0.77546 0.00013")
  LIST(APPEND BCC_H_CSVMC_SCALARS "deeenergy_AB" "-0.000603 0.000035")

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
  SET(FULL_NAME "cpu_limit-bccH_1x1x1_ae-rmc-1-1")
  SET(TEST_ADDED FALSE)
  RUN_QMC_APP(${FULL_NAME}
              "${CMAKE_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
              1 1
              TEST_ADDED
              qmc_cpu_limit_rmc.xml)

  IF (TEST_ADDED)
    SET_PROPERTY(TEST ${FULL_NAME} APPEND PROPERTY TIMEOUT 120)
    SET_PROPERTY(TEST ${FULL_NAME} APPEND PROPERTY PASS_REGULAR_EXPRESSION "Time limit reached for")
  ENDIF()

ELSE()
  MESSAGE("Skipping bccH_1x1x1_ae tests for VMC all-particle moves, CSVMC, and RMC because they are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

# Plane wave test of bccH (no pp, so faster than typical solids)
# Reduced number of blocks and samples compared to other short bccH tests
  LIST(APPEND BCC_H_PW_SCALARS "totenergy" "-1.834032 0.001")
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
  LIST(APPEND LONG_BCC_H_SCALARS "totenergy" "-1.834032 0.00001")
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
  LIST(APPEND LONG_BCC_H_DMC_SCALARS "totenergy" "-1.841204 0.000285")
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

