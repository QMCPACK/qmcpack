#
# HEG J2 RPA test
#

#
# Tests for HEG with 54 electrons and an RPA two body Jastrow
# See tests/heg/heg_54_J2rpa/README for full details
#

# real but noisy error bar is 0.012941 (23 blocks only)
# resurrect the test then create better reference data
list(APPEND HEG54J2RPA_SCALARS "totenergy" "-4.209725  0.02") # total energy

qmc_run_and_check(
  developer-heg_54_J2rpa
  "${qmcpack_SOURCE_DIR}/tests/heg/heg_54_J2rpa"
  heg.rpa.rs5.xml
  16
  1
  TRUE
  1
  HEG54J2RPA_SCALARS)

# Fixed-seed references after correcting accepted k-space densities and
# proposed gradients, checked against full recomputation (PR #6386).
if(QMC_COMPLEX)
  if(QMC_MIXED_PRECISION)
    list(APPEND DET_HEG54J2RPA_SCALARS "totenergy" "-3.89404867 0.000007")
    list(APPEND DET_HEG54J2RPA_SCALARS "kinetic" "3.43353584 0.000007")
    list(APPEND DET_HEG54J2RPA_SCALARS "potential" "-7.32758451 0.000005")
    list(APPEND DET_HEG54J2RPA_SCALARS "eeenergy" "-7.32758451 0.000005")
    list(APPEND DET_HEG54J2RPA_BATCHED_SCALARS "totenergy" "-3.89404920 0.000007")
    list(APPEND DET_HEG54J2RPA_BATCHED_SCALARS "kinetic" "3.43353510 0.000007")
    list(APPEND DET_HEG54J2RPA_BATCHED_SCALARS "potential" "-7.32758430 0.000005")
    list(APPEND DET_HEG54J2RPA_BATCHED_SCALARS "eeenergy" "-7.32758430 0.000005")
  else()
    list(APPEND DET_HEG54J2RPA_SCALARS "totenergy" "-3.89404977 0.000001")
    list(APPEND DET_HEG54J2RPA_SCALARS "kinetic" "3.43353327 0.000001")
    list(APPEND DET_HEG54J2RPA_SCALARS "potential" "-7.32758304 0.000001")
    list(APPEND DET_HEG54J2RPA_SCALARS "eeenergy" "-7.32758304 0.000001")
    list(APPEND DET_HEG54J2RPA_BATCHED_SCALARS "totenergy" "-3.89404977 0.000001")
    list(APPEND DET_HEG54J2RPA_BATCHED_SCALARS "kinetic" "3.43353327 0.000001")
    list(APPEND DET_HEG54J2RPA_BATCHED_SCALARS "potential" "-7.32758304 0.000001")
    list(APPEND DET_HEG54J2RPA_BATCHED_SCALARS "eeenergy" "-7.32758304 0.000001")
  endif()
  qmc_run_and_check(
    deterministic-heg_54_J2rpa
    "${qmcpack_SOURCE_DIR}/tests/heg/heg_54_J2rpa"
    det_heg.rpa.rs5.xml
    1
    1
    TRUE
    1
    DET_HEG54J2RPA_SCALARS)

  qmc_run_and_check(
    deterministic-heg_54_J2rpa-batched
    "${qmcpack_SOURCE_DIR}/tests/heg/heg_54_J2rpa"
    det_heg.rpa.rs5_batched.xml
    1
    1
    TRUE
    1
    DET_HEG54J2RPA_BATCHED_SCALARS)
else()
  ## real code should not be ran with a twist
  ##if(QMC_MIXED_PRECISION)
  ## list(APPEND DET_HEG54J2RPA_SCALARS "totenergy" "-3.92025567 0.000005")
  ## list(APPEND DET_HEG54J2RPA_SCALARS "kinetic" "3.63462568 0.000005")
  ## list(APPEND DET_HEG54J2RPA_SCALARS "potential" "-7.55488135 0.000005")
  ## list(APPEND DET_HEG54J2RPA_SCALARS "eeenergy" "-7.55488135 0.000005")
  ##else()
  ## list(APPEND DET_HEG54J2RPA_SCALARS "totenergy" "-3.92025567 0.000001")
  ## list(APPEND DET_HEG54J2RPA_SCALARS "kinetic" "3.63462568 0.000001")
  ## list(APPEND DET_HEG54J2RPA_SCALARS "potential" "-7.55488135 0.000001")
  ## list(APPEND DET_HEG54J2RPA_SCALARS "eeenergy" "-7.55488135 0.000001")
  ##endif()
endif()

get_property(CURRENT_TESTS DIRECTORY PROPERTY TESTS)
add_labels_for_tests("${CURRENT_TESTS}")
