
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
#
# 4 H atom cluster gaussian tests, all electron
# Also check results for different number of mpi tasks and threads keeping total constant
#
  LIST(APPEND H4_OPT_SCALARS "totenergy" "-2.140786 0.0009") # total energy

  QMC_RUN_AND_CHECK(short-H4-opt-linear
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-linear
                    optm-linear.xml
                    16 1
                    ${MP_SUCCESS}
                    5 H4_OPT_SCALARS # OPT step 5
                    )

  IF(BUILD_LMYENGINE_INTERFACE)
    QMC_RUN_AND_CHECK(short-H4-opt-adaptive
                      "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                      H4-adaptive
                      optm-adaptive.xml
                      16 1
                      ${MP_SUCCESS}
                      5 H4_OPT_SCALARS # OPT step 5
                      )

    QMC_RUN_AND_CHECK(short-H4-opt-adaptive
                      "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                      H4-adaptive
                      optm-adaptive.xml
                      1 16
                      ${MP_SUCCESS}
                      5 H4_OPT_SCALARS # OPT step 5
                      )
  ELSE()
   MESSAGE("Skipping H4_ae adaptive optimizer tests because lmyengine interface was not built (BUILD_LMYENGINE_INTERFACE=0)")
  ENDIF(BUILD_LMYENGINE_INTERFACE)

  QMC_RUN_AND_CHECK(short-H4-opt-OneShiftOnly
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-OneShiftOnly
                    optm-OneShiftOnly.xml
                    16 1
                    ${MP_SUCCESS}
                    5 H4_OPT_SCALARS # OPT step 5
                    )

  QMC_RUN_AND_CHECK(short-H4-opt-OneShiftOnly
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-OneShiftOnly
                    optm-OneShiftOnly.xml
                    1 16
                    ${MP_SUCCESS}
                    5 H4_OPT_SCALARS # OPT step 5
                    )
 ELSE()
    MESSAGE("Skipping H4_ae tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping H4_ae tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

