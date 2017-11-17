
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)

#
# H2 starting from unoptimized CIS coefficients and orbital
# Targetting the first excited state, and using the adaptive linear method
# The Jastrow is already optimized (by itself). This test optimizes both orbitals and CIS coefficients
#
  LIST(APPEND H2_FDLR_SCALARS "totenergy" "-0.64755 0.0015") # total energy

  IF(BUILD_LMYENGINE_INTERFACE)
    QMC_RUN_AND_CHECK(short-H2-FDLR
                      "${CMAKE_SOURCE_DIR}/tests/molecules/H2_FDLR_opt"
                      H2
                      h2_fdlr.xml
                      16 1
                      ${MP_SUCCESS}
                      6 H2_FDLR_SCALARS # Final VMC step
                      TRUE)
  ELSE()
   MESSAGE("Skipping H4_FDLR test because lmyengine interface was not built (BUILD_LMYENGINE_INTERFACE=0)")
  ENDIF(BUILD_LMYENGINE_INTERFACE)
 ELSE()
    MESSAGE("Skipping H4_FDLR tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping H4_FDLR tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

