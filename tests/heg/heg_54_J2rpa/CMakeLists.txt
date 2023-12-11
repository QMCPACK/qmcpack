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
  HEG.rs5.rpa
  heg.rpa.rs5.xml
  16
  1
  TRUE
  1
  HEG54J2RPA_SCALARS)

if(QMC_COMPLEX)
  if(QMC_MIXED_PRECISION)
    list(APPEND DET_HEG54J2RPA_SCALARS "totenergy" "-3.91513517 0.000005")
    list(APPEND DET_HEG54J2RPA_SCALARS "kinetic" "3.27073155 0.000005")
    list(APPEND DET_HEG54J2RPA_SCALARS "potential" "-7.18586671 0.000005")
    list(APPEND DET_HEG54J2RPA_SCALARS "eeenergy" "-7.18586671 0.000005")
  else()
    list(APPEND DET_HEG54J2RPA_SCALARS "totenergy" "-3.91514075 0.000001")
    list(APPEND DET_HEG54J2RPA_SCALARS "kinetic" "3.27072556 0.000001")
    list(APPEND DET_HEG54J2RPA_SCALARS "potential" "-7.18586631 0.000001")
    list(APPEND DET_HEG54J2RPA_SCALARS "eeenergy" "-7.18586631 0.000001")
  endif()
  qmc_run_and_check(
    deterministic-heg_54_J2rpa
    "${qmcpack_SOURCE_DIR}/tests/heg/heg_54_J2rpa"
    det_HEG.rs5.rpa
    det_heg.rpa.rs5.xml
    1
    1
    TRUE
    1
    DET_HEG54J2RPA_SCALARS)
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
