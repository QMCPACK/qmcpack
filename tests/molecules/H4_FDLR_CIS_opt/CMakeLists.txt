
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)

#
# H4 starting from unoptimized CIS coefficients, targetting the first excited state, and using the adaptive linear method
#
  LIST(APPEND H4_FDLR_SCALARS "totenergy" "-1.82507816 0.001") # total energy

  IF(BUILD_LMYENGINE_INTERFACE)
    QMC_RUN_AND_CHECK(short-H4-FDLR
                      "${CMAKE_SOURCE_DIR}/tests/molecules/H4_FDLR_CIS_opt"
                      H4
                      H4_FDLR_CIS_opt.xml
                      16 1
                      ${MP_SUCCESS}
                      8 H4_FDLR_SCALARS # Final VMC step
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

