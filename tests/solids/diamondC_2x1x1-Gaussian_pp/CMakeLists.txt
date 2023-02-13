# Reference VMC run in qmc-ref "-21.70152 +/- 0.00042"
list(APPEND DIAMOND_SCALARS "totenergy" "-21.70152 0.015")
list(APPEND DIAMOND_SCALARS "kinetic" "20.6385 0.16 ")
list(APPEND DIAMOND_SCALARS "potential" "-42.3400 0.19")
list(APPEND DIAMOND_SCALARS "eeenergy" "-5.33913 0.034")
list(APPEND DIAMOND_SCALARS "ionion" "-25.551326 0.0001")
list(APPEND DIAMOND_SCALARS "localecp" "-13.1852 0.20")
list(APPEND DIAMOND_SCALARS "nonlocalecp" "1.7356 0.05")
list(APPEND DIAMOND_SCALARS "samples" "8000 0.0")
#  LIST(APPEND DIAMOND_SCALARS "flux" "0.0 0.4")

qmc_run_and_check(
  short-diamondC_2x1x1_pp-vmc_gaussian_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
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

# Reference VMC run in qmc-ref "-21.70152 +/- 0.00042"
list(APPEND LONG_DIAMOND_SCALARS "totenergy" "-21.70152   0.008351")
list(APPEND LONG_DIAMOND_SCALARS "samples" "32000 0.0")

qmc_run_and_check(
  long-diamondC_2x1x1_pp-vmc_gaussian_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
  qmc_long
  qmc_long.in.xml
  1
  16
  TRUE
  0
  LONG_DIAMOND_SCALARS # VMC
)

# Reference DMC run in qmc-ref "-21.84712 +/-  0.00068"
list(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-21.84712 0.009024")
qmc_run_and_check(
  long-diamondC_2x1x1_pp-dmc_gaussian_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
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
  list(APPEND DET_DIAMOND_SCALARS "totenergy" "-22.7014044 0.000516")
  list(APPEND DET_DIAMOND_SCALARS "kinetic" "15.76862594 0.000516")
  list(APPEND DET_DIAMOND_SCALARS "potential" "-38.46978416 0.000025")
  list(APPEND DET_DIAMOND_SCALARS "eeenergy" "-5.67930456 0.000009")
  list(APPEND DET_DIAMOND_SCALARS "ionion" "-25.55133282 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "localecp" "-12.0143307 0.000051")
  list(APPEND DET_DIAMOND_SCALARS "nonlocalecp" "4.77523301 0.000026")
  list(APPEND DET_DIAMOND_SCALARS "samples" "9 0.0")
else()
  list(APPEND DET_DIAMOND_SCALARS "totenergy" "-22.70038582 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "kinetic" "15.76939694 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "potential" "-38.46978275 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "eeenergy" "-5.67929898 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "ionion" "-25.55132719 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "localecp" "-12.01437633 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "nonlocalecp" "4.77521942 0.000001")
  list(APPEND DET_DIAMOND_SCALARS "samples" "9 0.0")
endif()

qmc_run_and_check(
  deterministic-diamondC_2x1x1_pp-vmc_gaussian_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
  det_qmc_short
  det_qmc_short.in.xml
  1
  1
  TRUE
  0
  DET_DIAMOND_SCALARS # VMC
)
