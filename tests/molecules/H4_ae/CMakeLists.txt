
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

  QMC_RUN_AND_CHECK(short-H4-opt-linear-rescale
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-linear-rescale
                    optm-linear-rescale.xml
                    4 4
                    ${MP_SUCCESS}
                    5 H4_OPT_SCALARS # OPT step 5
                    )

  QMC_RUN_AND_CHECK(short-H4-opt-linear-linemin
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-linear-linemin
                    optm-linear-linemin.xml
                    16 1
                    ${MP_SUCCESS}
                    5 H4_OPT_SCALARS # OPT step 5
                    )

  QMC_RUN_AND_CHECK(short-H4-opt-cslinear-quartic
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-cslinear-quartic
                    optm-cslinear-quartic.xml
                    4 4
                    ${MP_SUCCESS}
                    7 H4_OPT_SCALARS # OPT step 7
                    )

  QMC_RUN_AND_CHECK(short-H4-opt-cslinear-rescale
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-cslinear-rescale
                    optm-cslinear-rescale.xml
                    16 1
                    ${MP_SUCCESS}
                    7 H4_OPT_SCALARS # OPT step 7
                    )

  QMC_RUN_AND_CHECK(short-H4-opt-cslinear-linemin
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-cslinear-linemin
                    optm-cslinear-linemin.xml
                    16 1
                    ${MP_SUCCESS}
                    7 H4_OPT_SCALARS # OPT step 7
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

    QMC_RUN_AND_CHECK(short-H4-RHF-orb-opt
                      "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                      H4
                      H4_RHF_orb_opt.xml
                      16 1
                      ${MP_SUCCESS}
                      14 H4_ORB_OPT_SCALARS # Final VMC step
                      TRUE)

    QMC_RUN_AND_CHECK(short-H4-orb-opt
                      "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                      H4
                      H4_orb_opt.xml
                      16 1
                      ${MP_SUCCESS}
                      14 H4_ORB_OPT_SCALARS # Final VMC step
                      TRUE)

    QMC_RUN_AND_CHECK(short-H4-orb-opt
                      "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                      H4
                      H4_orb_opt.xml
                      1 16
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
                      H4_orb_opt_dmc.xml
                      16 1
                      ${MP_SUCCESS}
                      3 H4_ORB_OPT_DMC_SCALARS # Test the DMC results
                      TRUE)

    QMC_RUN_AND_CHECK(short-H4-orb-opt-dmc
                      "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                      H4
                      H4_orb_opt_dmc.xml
                      1 16
                      ${MP_SUCCESS}
                      3 H4_ORB_OPT_DMC_SCALARS # Test the DMC results
                      TRUE)
    #
    # H4 starting from unoptimized CIS coefficients, targeting the first excited state, and using the adaptive linear method
    #

    LIST(APPEND H4_CJS_SCALARS "totenergy" "-2.0438 0.0045") # total energy

    QMC_RUN_AND_CHECK(short-H4-cjs
                      "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                      H4
                      H4_cjs_vmc.xml
                      16 1
                      ${MP_SUCCESS}
                      0 H4_CJS_SCALARS
                      TRUE)


  ELSE()
   MESSAGE("Skipping H4_ae adaptive optimizer tests because lmyengine interface was not built (BUILD_LMYENGINE_INTERFACE=0)")
  ENDIF(BUILD_LMYENGINE_INTERFACE)

  LIST(APPEND H4_DESCENT_SCALARS "totenergy" "-2.047099 0.005")

  QMC_RUN_AND_CHECK(short-H4-descent
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-descent
                    optm-descent.xml
                    16 1
                    ${MP_SUCCESS}
                    9 H4_DESCENT_SCALARS
                    TRUE)

  QMC_RUN_AND_CHECK(short-H4-descent
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-descent
                    optm-descent.xml
                    1 16
                    ${MP_SUCCESS}
                    9 H4_DESCENT_SCALARS
                    TRUE)

  LIST(APPEND H4_TARGET_DESCENT_SCALARS "totenergy" "-2.04604154 0.003")

  QMC_RUN_AND_CHECK(short-H4-target-descent
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-target-descent
                    optm-target-descent.xml
                    16 1
                    ${MP_SUCCESS}
                    9 H4_TARGET_DESCENT_SCALARS
                    TRUE)

  QMC_RUN_AND_CHECK(short-H4-target-descent
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-target-descent
                    optm-target-descent.xml
                    1 16
                    ${MP_SUCCESS}
                    9 H4_TARGET_DESCENT_SCALARS
                    TRUE)

  LIST(APPEND H4_HYBRID_SCALARS "totenergy" "-2.140736 0.003")

  QMC_RUN_AND_CHECK(short-H4-hybrid
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-hybrid
                    optm-hybrid.xml
                    16 1
                    ${MP_SUCCESS}
		    9 H4_HYBRID_SCALARS
                    TRUE)

  QMC_RUN_AND_CHECK(short-H4-hybrid
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-hybrid
                    optm-hybrid.xml
                    1 16
                    ${MP_SUCCESS}
		    9 H4_HYBRID_SCALARS
                    TRUE)
  
  LIST(APPEND H4_TARGET_HYBRID_SCALARS "totenergy" "-2.13527301 0.003")

  QMC_RUN_AND_CHECK(short-H4-target-hybrid
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-target-hybrid
                    optm-target-hybrid.xml
                    16 1
                    ${MP_SUCCESS}
                    9 H4_TARGET_HYBRID_SCALARS
                    TRUE)

  QMC_RUN_AND_CHECK(short-H4-target-hybrid
                    "${CMAKE_SOURCE_DIR}/tests/molecules/H4_ae"
                    H4-target-hybrid
                    optm-target-hybrid.xml
                    1 16
                    ${MP_SUCCESS}
                    9 H4_TARGET_HYBRID_SCALARS
                    TRUE)

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
    MESSAGE_VERBOSE("Skipping H4_ae tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE_VERBOSE("Skipping H4_ae tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

