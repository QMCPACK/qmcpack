
IF (NOT QMC_CUDA)
IF (NOT QMC_COMPLEX)
#
# H2 test in a DZV basis set, using optimizable determinants, and the adaptive linear method
#
  LIST(APPEND H2_OPT_SCALARS "totenergy" "-1.124563 0.00084") # total energy

IF(BUILD_LMYENGINE_INTERFACE)
  QMC_RUN_AND_CHECK(short-H2-orb-opt
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H2_orb_opt"
                    H2
                    h2_orb_opt.xml
                    16 1
                    H2_OPT_SCALARS
                    10 # Final VMC step
                    TRUE)
ENDIF(BUILD_LMYENGINE_INTERFACE)
ENDIF()
ENDIF()
