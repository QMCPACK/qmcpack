
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)

#
# H4 starting from perturbed orbitals, using optimizable determinants, and the adaptive linear method
#
  LIST(APPEND H4_ORB_OPT_SCALARS "totenergy" "-1.93045653 0.00837934") # total energy

  IF(BUILD_LMYENGINE_INTERFACE)
    QMC_RUN_AND_CHECK(short-H4-orb-opt
                      "${CMAKE_SOURCE_DIR}/tests/molecules/H4_orb_opt"
                      H4
                      h4_orb_opt.xml
                      16 1
                      ${MP_SUCCESS}
                      14 H4_ORB_OPT_SCALARS # Final VMC step
                      TRUE)
  ELSE()
   MESSAGE("Skipping H4_orb_opt test because lmyengine interface was not built (BUILD_LMYENGINE_INTERFACE=0)")
  ENDIF(BUILD_LMYENGINE_INTERFACE)
 ELSE()
    MESSAGE("Skipping H4_orb_opt tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping H4_orb_opt tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

