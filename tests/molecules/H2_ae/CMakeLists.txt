
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
  IF(BUILD_LMYENGINE_INTERFACE)

  #
  # H2 test in a DZV basis set, using optimizable determinants, and the adaptive linear method
  #
  LIST(APPEND H2_OPT_SCALARS "totenergy" "-1.124563 0.00084") # total energy

  QMC_RUN_AND_CHECK(short-H2-orb-opt
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H2_ae"
                    H2
                    h2_orb_opt.xml
                    16 1
                    ${SUCCESS_STATUS_MP}
                    10 H2_OPT_SCALARS # Final VMC step
                    TRUE)

  ELSE()
    MESSAGE_VERBOSE("Skipping H2 test because lmyengine interface was not built (BUILD_LMYENGINE_INTERFACE=0)")
  ENDIF(BUILD_LMYENGINE_INTERFACE)
 ELSE()
    MESSAGE_VERBOSE("Skipping H2 tests because orbital rotation is not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE_VERBOSE("Skipping H2 tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

