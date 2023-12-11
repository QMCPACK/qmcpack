#
# Short tests, about 2-3 seconds on 16 cores
#
list(APPEND BCC_H_SCALARS "totenergy" "-1.834032 0.0003")
list(APPEND BCC_H_SCALARS "kinetic" "0.186673807952 0.003")
list(APPEND BCC_H_SCALARS "potential" "-2.02070672998 0.003")
list(APPEND BCC_H_SCALARS "eeenergy" "-0.775870180954 0.003")
list(APPEND BCC_H_SCALARS "ionion" "-0.96289961022 0.001")
list(APPEND BCC_H_SCALARS "samples" "960000 0.0")
#  LIST(APPEND BCC_H_SCALARS "flux" "0.0 0.03")

qmc_run_and_check(
  short-bccH_1x1x1_ae-vmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
  qmc_short
  qmc_short.in.xml
  1
  16
  TRUE
  0
  BCC_H_SCALARS # VMC
)

# qmc_ref_vmc_dmc  series 1  -1.841352 +/- 0.000078
# qmc_short_vmc_dmc  series 1  -1.84196 +/- 0.00073
list(APPEND BCC_H_DMC_SCALARS "totenergy" "-1.84196 0.00073")
list(APPEND BCC_H_DMC_SCALARS "potential" "-2.0187  0.004")
list(APPEND BCC_H_DMC_SCALARS "kinetic" "0.1767  0.004")
list(APPEND BCC_H_DMC_SCALARS "eeenergy" "-0.77985 0.00069")
qmc_run_and_check(
  short-bccH_1x1x1_ae-dmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
  qmc_short_vmc_dmc
  qmc_short_vmc_dmc.in.xml
  1
  16
  TRUE
  1
  BCC_H_DMC_SCALARS # DMC
)

# Check VMCUpdateAll driver using all-particle updates
qmc_run_and_check(
  short-bccH_1x1x1_ae-vmc-all-nodrift_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
  qmc_short_vmc_all_nodrift
  qmc_short_vmc_all_nodrift.in.xml
  1
  16
  TRUE
  0
  BCC_H_SCALARS # VMC
)

# Check VMCUpdateAllWithDrift driver using all-particle updates
qmc_run_and_check(
  short-bccH_1x1x1_ae-vmc-all_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
  qmc_short_vmc_all
  qmc_short_vmc_all.in.xml
  1
  16
  TRUE
  0
  BCC_H_SCALARS # VMC
)

# Check DMCUpdateAllWithRejection driver using all-particle updates
qmc_run_and_check(
  short-bccH_1x1x1_ae-dmc-all_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
  qmc_short_vmc_dmc_all
  qmc_short_vmc_dmc_all.in.xml
  1
  16
  TRUE
  1
  BCC_H_DMC_SCALARS # DMC
)

# Deterministic correlated sampling test
list(APPEND DET_BCC_H_CSVMC_SCALARS "totenergy_A" "-1.5879406376 0.00001")
list(APPEND DET_BCC_H_CSVMC_SCALARS "totenergy_B" "-1.5845867400 0.00001")
list(APPEND DET_BCC_H_CSVMC_SCALARS "dtotenergy_AB" "-0.0033713724 0.00001")
list(APPEND DET_BCC_H_CSVMC_SCALARS "ionion_A" "-0.9628996202 0.0000001")
list(APPEND DET_BCC_H_CSVMC_SCALARS "ionion_B" "-0.9581655983 0.0000001")
list(APPEND DET_BCC_H_CSVMC_SCALARS "dionion_AB" "-0.0047340219 0.0000001")

qmc_run_and_check(
  deterministic-bccH_1x1x1_ae-csvmc-all_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
  det_qmc_short_csvmc_all
  det_qmc_short_csvmc_all.in.xml
  1
  1
  TRUE
  0
  DET_BCC_H_CSVMC_SCALARS # VMC
)

#Now for the correlated sampling tests
list(APPEND BCC_H_CSVMC_SCALARS "totenergy_A" "-1.83428 0.00013")
list(APPEND BCC_H_CSVMC_SCALARS "totenergy_B" "-1.83390 0.00012")
list(APPEND BCC_H_CSVMC_SCALARS "dtotenergy_AB" "-0.000384 0.000068")
list(APPEND BCC_H_CSVMC_SCALARS "ionion_A" "-0.96289962 0.0000001")
list(APPEND BCC_H_CSVMC_SCALARS "ionion_B" "-0.95816559 0.0000001")
list(APPEND BCC_H_CSVMC_SCALARS "dionion_AB" "-0.004734021 0.0000001")
list(APPEND BCC_H_CSVMC_SCALARS "eeenergy_A" "-0.77607 0.00013")
list(APPEND BCC_H_CSVMC_SCALARS "eeenergy_B" "-0.77546 0.00013")
list(APPEND BCC_H_CSVMC_SCALARS "deeenergy_AB" "-0.000603 0.000035")

# Check CSVMCUpdatePbyP driver using particle-by-particle updates
qmc_run_and_check(
  short-bccH_1x1x1_ae-csvmc-pbyp-nodrift_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
  qmc_short_csvmc_pbyp_nodrift
  qmc_short_csvmc_pbyp_nodrift.in.xml
  16
  1
  TRUE
  0
  BCC_H_CSVMC_SCALARS # VMC
)

# Check CSVMCUpdatePbyPWithDrift driver using all-particle updates
#  QMC_RUN_AND_CHECK(short-bccH_1x1x1_ae-csvmc-pbyp_sdj
#                    "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
#                    qmc_short_csvmc_pbyp
#                    qmc_short_csvmc_pbyp.in.xml
#                    16 1
#                    TRUE
#                    0 BCC_H_CSVMC_SCALARS # VMC
#                    )

# Check CSVMCUpdateAll driver using all-particle updates
qmc_run_and_check(
  short-bccH_1x1x1_ae-csvmc-all-nodrift_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
  qmc_short_csvmc_all_nodrift
  qmc_short_csvmc_all_nodrift.in.xml
  16
  1
  TRUE
  0
  BCC_H_CSVMC_SCALARS # VMC
)

# Check CSVMCUpdateAllWithDrift driver using all-particle updates
qmc_run_and_check(
  short-bccH_1x1x1_ae-csvmc-all_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
  qmc_short_csvmc_all
  qmc_short_csvmc_all.in.xml
  16
  1
  TRUE
  0
  BCC_H_CSVMC_SCALARS # VMC
)

list(APPEND BCC_H_RMC_SCALARS "totenergy" "-1.8420 0.0011")
#check RMC on 16 cores.  1 MPI
qmc_run_and_check(
  short-bccH_1x1x1_ae-rmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
  qmc_short_vmc_rmc
  qmc_short_vmc_rmc.in.xml
  1
  16
  TRUE
  1
  BCC_H_RMC_SCALARS #series#
)
#check RMC on 4 cores, 4 mpi
qmc_run_and_check(
  short-bccH_1x1x1_ae-rmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
  qmc_short_vmc_rmc
  qmc_short_vmc_rmc.in.xml
  4
  4
  TRUE
  1
  BCC_H_RMC_SCALARS #series#
)

# Check CPU time limit
cpu_limit_run("bccH_1x1x1_ae-rmc-1-1" "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae" 1 1 120 qmc_cpu_limit_rmc.xml)

# Plane wave test of bccH (no pp, so faster than typical solids)
# Reduced number of blocks and samples compared to other short bccH tests
list(APPEND BCC_H_PW_SCALARS "totenergy" "-1.834032 0.001")
list(APPEND BCC_H_PW_SCALARS "samples" "96000 0.0")
#  LIST(APPEND BCC_H_PW_SCALARS "flux" "0.0 0.1")
# No GPU implementation. Would be fast if implemented
qmc_run_and_check(
  short-bccH_1x1x1_ae-vmc_sdj-pw
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
  qmc_short_pw
  qmc_short_pw.in.xml
  1
  16
  TRUE
  0
  BCC_H_PW_SCALARS # VMC
)

#
# Long variants of above, about 1/2 hour with 16 cores
#
list(APPEND LONG_BCC_H_SCALARS "totenergy" "-1.834032 0.00001")
list(APPEND LONG_BCC_H_SCALARS "samples" "921600000 0.0")

qmc_run_and_check(
  long-bccH_1x1x1_ae-vmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
  qmc_long
  qmc_long.in.xml
  1
  16
  TRUE
  0
  LONG_BCC_H_SCALARS # VMC
)

# qmc_ref_vmc_dmc  series 1  -1.841352 +/- 0.000078
# qmc_long_vmc_dmc  series 1   -1.841204 +/- 0.000285
list(APPEND LONG_BCC_H_DMC_SCALARS "totenergy" "-1.841204 0.000285")
qmc_run_and_check(
  long-bccH_1x1x1_ae-dmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
  qmc_long_vmc_dmc
  qmc_long_vmc_dmc.in.xml
  1
  16
  TRUE
  1
  LONG_BCC_H_DMC_SCALARS # DMC
)

list(APPEND LONG_BCC_H_RMC_SCALARS "totenergy" "-1.84150 0.00011")
list(APPEND LONG_BCC_H_RMC_SCALARS "potential" "-2.01964 0.00072")
list(APPEND LONG_BCC_H_RMC_SCALARS "potential_pure" "-2.01907 0.00085")
list(APPEND LONG_BCC_H_RMC_SCALARS "kinetic_mixed" "0.17814 0.00073")
list(APPEND LONG_BCC_H_RMC_SCALARS "kinetic_pure" "0.17063 0.00086")
list(APPEND LONG_BCC_H_RMC_SCALARS "eeenergy_mixed" "-0.77991 0.00016")
list(APPEND LONG_BCC_H_RMC_SCALARS "eeenergy_pure" "-0.78335 0.00019")
qmc_run_and_check(
  long-bccH_1x1x1_ae-rmc_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/bccH_1x1x1_ae"
  qmc_long_vmc_rmc
  qmc_long_vmc_rmc.in.xml
  1
  16
  TRUE
  1
  LONG_BCC_H_RMC_SCALARS # DMC
)
