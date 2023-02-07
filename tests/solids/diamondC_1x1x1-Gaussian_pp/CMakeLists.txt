# Reference VMC run in qmc-ref "-10.495941  0.000073"
list(APPEND DIAMOND_SCALARS "totenergy" "-10.495941 0.0065")
list(APPEND DIAMOND_SCALARS "kinetic" "11.51198 0.078")
list(APPEND DIAMOND_SCALARS "potential" "-22.00792 0.079")
list(APPEND DIAMOND_SCALARS "eeenergy" "-2.68548 0.018")
list(APPEND DIAMOND_SCALARS "ionion" "-12.77566 0.0001")
list(APPEND DIAMOND_SCALARS "localecp" "-7.2011 0.090")
list(APPEND DIAMOND_SCALARS "nonlocalecp" "0.65437 0.024")
list(APPEND DIAMOND_SCALARS "samples" "16000 0.0")

qmc_run_and_check(
  short-diamondC_1x1x1_pp-vmc_gaussian_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
  qmc_short
  qmc_short.in.xml
  1
  16
  TRUE
  0
  DIAMOND_SCALARS # VMC
)

#
# Long tests
#

list(APPEND LONG_DIAMOND_SCALARS "totenergy" "-10.495941 0.0017")
list(APPEND LONG_DIAMOND_SCALARS "samples" "160000 0.0")
qmc_run_and_check(
  long-diamondC_1x1x1_pp-vmc_gaussian_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
  qmc_long
  qmc_long.in.xml
  1
  16
  TRUE
  0
  LONG_DIAMOND_SCALARS # VMC
)

# Reference DMC run in qmc-ref "-10.532497 0.000470"
list(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-10.532497 0.006127")

qmc_run_and_check(
  long-diamondC_1x1x1_pp-dmc_gaussian_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
  qmc_long_vmc_dmc
  qmc_long_vmc_dmc.in.xml
  1
  16
  TRUE
  1
  LONG_DIAMOND_DMC_SCALARS # DMC
)

# Deterministic test
if(QMC_MIXED_PRECISION)
  list(APPEND DET_DIAMOND_SCALARS "totenergy" "-10.26910069 0.000222")
  list(APPEND DET_DIAMOND_SCALARS "kinetic" "10.15426733 0.000221")
  list(APPEND DET_DIAMOND_SCALARS "potential" "-20.42328009 0.000005")
  list(APPEND DET_DIAMOND_SCALARS "eeenergy" "-1.329163913 0.000003")
  list(APPEND DET_DIAMOND_SCALARS "ionion" "-12.77567012 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "localecp" "-6.762353003 0.000007")
  list(APPEND DET_DIAMOND_SCALARS "nonlocalecp" "0.44390473 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "samples" "9 0.0")
else()
  list(APPEND DET_DIAMOND_SCALARS "totenergy" "-10.26894006 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "kinetic" "10.15433651 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "potential" "-20.42327657 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "eeenergy" "-1.32916126 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "ionion" "-12.77566743 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "localecp" "-6.76235283 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "nonlocalecp" "0.44390498 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "samples" "9 0.0")
endif()

qmc_run_and_check(
  deterministic-diamondC_1x1x1_pp-vmc_gaussian_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
  det_qmc_short
  det_qmc_short.in.xml
  1
  1
  TRUE
  0
  DET_DIAMOND_SCALARS # VMC
)
