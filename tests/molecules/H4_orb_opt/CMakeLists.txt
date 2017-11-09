
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
                    H4_ORB_OPT_SCALARS
                    14 # Final VMC step
                    TRUE)
ENDIF(BUILD_LMYENGINE_INTERFACE)

ENDIF()
ENDIF()
