
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
    #
    # H4 starting from perturbed orbitals, using optimizable determinants, and the adaptive linear method
    #
    LIST(APPEND H4_ORB_OPT_SCALARS "totenergy" "-2.0889 0.001") # total energy

    QMC_RUN_AND_CHECK(short-H4-orb-opt
                      "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                      H4
                      h4_orb_opt.xml
                      16 1
                      ${MP_SUCCESS}
                      14 H4_ORB_OPT_SCALARS # Final VMC step
                      TRUE)

    #
    # H4 starting from converged wave function, doing a quick adaptive linear method run (using the same
    # spo set for both spin-up and spin-down determinants) and then doing a DMC calculation
    #
    LIST(APPEND H4_ORB_OPT_DMC_SCALARS "totenergy" "-2.191461 0.0048") # total energy

    QMC_RUN_AND_CHECK(short-H4-orb-opt-dmc
                      "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                      H4
                      h4_orb_opt_dmc.xml
                      16 1
                      ${MP_SUCCESS}
                      3 H4_ORB_OPT_DMC_SCALARS # Test the DMC results
                      TRUE)

    QMC_RUN_AND_CHECK(short-H4-orb-opt-dmc
                      "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                      H4
                      h4_orb_opt_dmc.xml
                      1 16
                      ${MP_SUCCESS}
                      3 H4_ORB_OPT_DMC_SCALARS # Test the DMC results
                      TRUE)
    #
    # H4 starting from unoptimized CIS coefficients, targetting the first excited state, and using the adaptive linear method
    #
    LIST(APPEND H4_FDLR_SCALARS "totenergy" "-1.82507816 0.001") # total energy

    QMC_RUN_AND_CHECK(short-H4-FDLR
                      "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                      H4
                      H4_FDLR_CIS_opt.xml
                      16 1
                      ${MP_SUCCESS}
                      8 H4_FDLR_SCALARS # Final VMC step
                      TRUE)

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

