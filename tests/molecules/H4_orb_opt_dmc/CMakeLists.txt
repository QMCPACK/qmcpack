
IF (NOT QMC_CUDA)
IF (NOT QMC_COMPLEX)
#
# H4 starting from converged wave function, doing a quick adaptive linear method run (using the same
# spo set for both spin-up and spin-down determinants) and then doing a DMC calculation
#
  LIST(APPEND H4_ORB_OPT_DMC_SCALARS "totenergy" "-2.191461 0.0048") # total energy

IF(BUILD_LMYENGINE_INTERFACE)
  QMC_RUN_AND_CHECK(short-H4-orb-opt-dmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_orb_opt_dmc"
                    H4
                    h4_orb_opt.xml
                    16 1
                    H4_ORB_OPT_DMC_SCALARS
                    3 # Test the DMC results
                    TRUE)

  QMC_RUN_AND_CHECK(short-H4-orb-opt-dmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_orb_opt_dmc"
                    H4
                    h4_orb_opt.xml
                    1 16
                    H4_ORB_OPT_DMC_SCALARS
                    3 # Test the DMC results
                    TRUE)
ENDIF(BUILD_LMYENGINE_INTERFACE)

ENDIF()
ENDIF()
