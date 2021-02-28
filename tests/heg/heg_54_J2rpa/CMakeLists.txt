

#
# HEG J2 RPA test
#
IF (NOT QMC_CUDA)

#
# Tests for HEG with 54 electrons and an RPA two body Jastrow
# See tests/heg/heg_54_J2rpa/README for full details
#

  # real but noisy error bar is 0.012941 (23 blocks only)
  # resurrect the test then create better reference data
  LIST(APPEND HEG54J2RPA_SCALARS "totenergy" "-4.209725  0.02") # total energy

  QMC_RUN_AND_CHECK(developer-heg_54_J2rpa
                    "${CMAKE_SOURCE_DIR}/tests/heg/heg_54_J2rpa"
                    HEG.rs5.rpa
                    heg.rpa.rs5.xml
                    16 1
                    TRUE
                    1 HEG54J2RPA_SCALARS
                    )

ELSE()
    MESSAGE("Skipping HEG J2 RPA test because they are not supported by CUDA build (QMC_CUDA=1)")  
ENDIF()

