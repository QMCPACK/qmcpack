
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
   IF(NOT QMC_MIXED_PRECISION)
#
#  Co atom L2 tests
#

#
# Deterministic tests
#

  LIST(APPEND DET_CO_L2_STD_VMC_1_1_SCALARS "totenergy" "-144.46113333 0.000001")
  LIST(APPEND DET_CO_L2_STD_VMC_1_1_SCALARS "kinetic"   "  52.20699764 0.000001")
  LIST(APPEND DET_CO_L2_STD_VMC_1_1_SCALARS "potential" "-196.66813096 0.000001")
  LIST(APPEND DET_CO_L2_STD_VMC_1_1_SCALARS "localecp"  "-266.85397100 0.000001")
  LIST(APPEND DET_CO_L2_STD_VMC_1_1_SCALARS "samples"   " 768.00000000 0.000001")

  LIST(APPEND DET_CO_L2_STD_DMC_1_1_SCALARS "totenergy" "-145.07726145 0.000001")
  LIST(APPEND DET_CO_L2_STD_DMC_1_1_SCALARS "kinetic"   "  75.69234768 0.000001")
  LIST(APPEND DET_CO_L2_STD_DMC_1_1_SCALARS "potential" "-220.76960913 0.000001")
  LIST(APPEND DET_CO_L2_STD_DMC_1_1_SCALARS "localecp"  "-297.60619331 0.000001")
  LIST(APPEND DET_CO_L2_STD_DMC_1_1_SCALARS "samples"   "1127.00000000 0.000001")

  QMC_RUN_AND_CHECK(deterministic-Co_L2_std-vmc_dmc
                    "${qmcpack_SOURCE_DIR}/tests/molecules/Co_L2"
                    det_dmc_std_diff
                    det_dmc_std_diff.in.xml
                    1 1
                    TRUE
                    0 DET_CO_L2_STD_VMC_1_1_SCALARS # VMC
                    1 DET_CO_L2_STD_DMC_1_1_SCALARS # DMC
                    )

  LIST(APPEND DET_CO_L2_L2_VMC_1_1_SCALARS "totenergy" "-144.46113333 0.000001")
  LIST(APPEND DET_CO_L2_L2_VMC_1_1_SCALARS "kinetic"   "  52.20699764 0.000001")
  LIST(APPEND DET_CO_L2_L2_VMC_1_1_SCALARS "potential" "-196.66813096 0.000001")
  LIST(APPEND DET_CO_L2_L2_VMC_1_1_SCALARS "localecp"  "-266.85397100 0.000001")
  LIST(APPEND DET_CO_L2_L2_VMC_1_1_SCALARS "samples"   " 768.00000000 0.000001")

  LIST(APPEND DET_CO_L2_L2_DMC_1_1_SCALARS "totenergy" "-145.41190911 0.000001")
  LIST(APPEND DET_CO_L2_L2_DMC_1_1_SCALARS "kinetic"   "  69.14536606 0.000001")
  LIST(APPEND DET_CO_L2_L2_DMC_1_1_SCALARS "potential" "-214.55727517 0.000001")
  LIST(APPEND DET_CO_L2_L2_DMC_1_1_SCALARS "localecp"  "-291.18582120 0.000001")
  LIST(APPEND DET_CO_L2_L2_DMC_1_1_SCALARS "samples"   "1144.00000000 0.000001")

  QMC_RUN_AND_CHECK(deterministic-Co_L2_L2-vmc_dmc
                    "${qmcpack_SOURCE_DIR}/tests/molecules/Co_L2"
                    det_dmc_L2_diff
                    det_dmc_L2_diff.in.xml
                    1 1
                    TRUE
                    0 DET_CO_L2_L2_VMC_1_1_SCALARS # VMC
                    1 DET_CO_L2_L2_DMC_1_1_SCALARS # DMC
                    )


  ENDIF()

 ELSE()
    MESSAGE_VERBOSE("Skipping Co_L2 test because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE_VERBOSE("Skipping Co_L2 test because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

