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

qmc_run_and_check(
  short-diamondC_2x1x1_pp-vmcbatch_gaussian_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
  qmc_short_vmcbatch
  qmc_short_vmcbatch.in.xml
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

qmc_run_and_check(
  long-diamondC_2x1x1_pp-vmcbatch_gaussian_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
  qmc_long_vmcbatch
  qmc_long_vmcbatch.in.xml
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

qmc_run_and_check(
  long-diamondC_2x1x1_pp-dmcbatch_gaussian_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
  qmc_long_vmcbatch_dmcbatch
  qmc_long_vmcbatch_dmcbatch.in.xml
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

qmc_run_and_check(
  deterministic-diamondC_2x1x1_pp-vmcbatch_gaussian_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
  det_qmc_short_vmcbatch
  det_qmc_short_vmcbatch.in.xml
  1
  1
  TRUE
  0
  DET_DIAMOND_SCALARS # VMC
)

if(QMC_MIXED_PRECISION)
  list(APPEND DET_DIAMOND-r1-t4_SCALARS "totenergy" "-22.32437951 0.000516")
  list(APPEND DET_DIAMOND-r1-t4_SCALARS "kinetic" "20.77633955 0.000516")
  list(APPEND DET_DIAMOND-r1-t4_SCALARS "potential" "-43.10071906 0.000025")
  list(APPEND DET_DIAMOND-r1-t4_SCALARS "eeenergy" "-5.57610243 0.000009")
  list(APPEND DET_DIAMOND-r1-t4_SCALARS "ionion" "-25.55133196 0.000001")
  list(APPEND DET_DIAMOND-r1-t4_SCALARS "localecp" "-12.98422021 0.000051")
  list(APPEND DET_DIAMOND-r1-t4_SCALARS "nonlocalecp" "1.01093553 0.000026")
else()
  list(APPEND DET_DIAMOND-r1-t4_SCALARS "totenergy" "-21.89379749 0.000001")
  list(APPEND DET_DIAMOND-r1-t4_SCALARS "kinetic" "21.08911929 0.000001")
  list(APPEND DET_DIAMOND-r1-t4_SCALARS "potential" "-42.98291678 0.000001")
  list(APPEND DET_DIAMOND-r1-t4_SCALARS "eeenergy" "-5.70907337 0.000001")
  list(APPEND DET_DIAMOND-r1-t4_SCALARS "ionion" "-25.55132697 0.000001")
  list(APPEND DET_DIAMOND-r1-t4_SCALARS "localecp" "-13.72628312 0.000001")
  list(APPEND DET_DIAMOND-r1-t4_SCALARS "nonlocalecp" "2.00376669 0.000001")
endif()

qmc_run_and_check(
  deterministic-diamondC_2x1x1_pp-vmcbatch-dmcbatch_gaussian_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
  det_qmc_short_vmcbatch_dmcbatch
  det_qmc_short_vmcbatch_dmcbatch.in.xml
  1
  4
  TRUE
  1
  DET_DIAMOND-r1-t4_SCALARS #DMC
)  

if(QMC_MIXED_PRECISION)
  list(APPEND DET_DIAMOND-r2-t2_SCALARS "totenergy" "-21.94662338 0.000516")
  list(APPEND DET_DIAMOND-r2-t2_SCALARS "kinetic" "19.71552170 0.000516")
  list(APPEND DET_DIAMOND-r2-t2_SCALARS "potential" "-41.66214508 0.000025")
  list(APPEND DET_DIAMOND-r2-t2_SCALARS "eeenergy" "-5.30987675 0.000009")
  list(APPEND DET_DIAMOND-r2-t2_SCALARS "ionion" "-25.55133196 0.000001")
  list(APPEND DET_DIAMOND-r2-t2_SCALARS "localecp" "-11.87056876 0.000051")
  list(APPEND DET_DIAMOND-r2-t2_SCALARS "nonlocalecp" "1.06963239 0.000026")
else()
  list(APPEND DET_DIAMOND-r2-t2_SCALARS "totenergy" "-21.59886973 0.000001")
  list(APPEND DET_DIAMOND-r2-t2_SCALARS "kinetic" "20.69547774 0.000001")
  list(APPEND DET_DIAMOND-r2-t2_SCALARS "potential" "-42.29434746 0.000001")
  list(APPEND DET_DIAMOND-r2-t2_SCALARS "eeenergy" "-6.07824372 0.000001")
  list(APPEND DET_DIAMOND-r2-t2_SCALARS "ionion" "-25.55132697 0.000001")
  list(APPEND DET_DIAMOND-r2-t2_SCALARS "localecp" "-12.41866755 0.000001")
  list(APPEND DET_DIAMOND-r2-t2_SCALARS "nonlocalecp" "1.75389079 0.000001")
endif()

qmc_run_and_check(
  deterministic-diamondC_2x1x1_pp-vmcbatch-dmcbatch_gaussian_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
  det_qmc_short_vmcbatch_dmcbatch
  det_qmc_short_vmcbatch_dmcbatch.in.xml
  2
  2
  TRUE
  1
  DET_DIAMOND-r2-t2_SCALARS #DMC
)

if(QMC_MIXED_PRECISION)
  list(APPEND DET_DIAMOND-r2-t3_SCALARS "totenergy" "-21.86219836 0.000516")
  list(APPEND DET_DIAMOND-r2-t3_SCALARS "kinetic" "19.82298133 0.000516")
  list(APPEND DET_DIAMOND-r2-t3_SCALARS "potential" "-41.68517969 0.000025")
  list(APPEND DET_DIAMOND-r2-t3_SCALARS "eeenergy" "-5.83941772 0.000009")
  list(APPEND DET_DIAMOND-r2-t3_SCALARS "ionion" "-25.55133196 0.000001")
  list(APPEND DET_DIAMOND-r2-t3_SCALARS "localecp" "-13.02546579 0.000051")
  list(APPEND DET_DIAMOND-r2-t3_SCALARS "nonlocalecp" "2.73103577 0.000026")
else()
  list(APPEND DET_DIAMOND-r2-t3_SCALARS "totenergy" "-21.96961866 0.000001")
  list(APPEND DET_DIAMOND-r2-t3_SCALARS "kinetic" "18.90428814 0.000001")
  list(APPEND DET_DIAMOND-r2-t3_SCALARS "potential" "-40.87390680 0.000001")
  list(APPEND DET_DIAMOND-r2-t3_SCALARS "eeenergy" "-4.71429887 0.000001")
  list(APPEND DET_DIAMOND-r2-t3_SCALARS "ionion" "-25.55132697 0.000001")
  list(APPEND DET_DIAMOND-r2-t3_SCALARS "localecp" "-12.72947287 0.000001")
  list(APPEND DET_DIAMOND-r2-t3_SCALARS "nonlocalecp" "2.12119191 0.000001")
endif()

qmc_run_and_check(
  deterministic-diamondC_2x1x1_pp-vmcbatch-dmcbatch_gaussian_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp"
  det_qmc_short_vmcbatch_dmcbatch
  det_qmc_short_vmcbatch_dmcbatch.in.xml
  2
  3
  TRUE
  1
  DET_DIAMOND-r2-t3_SCALARS #DMC
)

