if(NOT QMC_COMPLEX)
  #
  # 4 H atom cluster gaussian tests, all electron
  # Also check results for different number of mpi tasks and threads keeping total constant
  #
  list(APPEND H4_OPT_SCALARS "totenergy" "-2.140786 0.0009") # total energy

  qmc_run_and_check(
    short-H4-opt-linear
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-linear
    optm-linear.xml
    16
    1
    ${SUCCESS_STATUS_MP}
    5
    H4_OPT_SCALARS # OPT step 5
  )

  qmc_run_and_check(
    short-H4-opt-linear-rescale
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-linear-rescale
    optm-linear-rescale.xml
    4
    4
    ${SUCCESS_STATUS_MP}
    5
    H4_OPT_SCALARS # OPT step 5
  )

  qmc_run_and_check(
    short-H4-opt-linear-linemin
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-linear-linemin
    optm-linear-linemin.xml
    16
    1
    ${SUCCESS_STATUS_MP}
    5
    H4_OPT_SCALARS # OPT step 5
  )

  if(BUILD_LMYENGINE_INTERFACE)
    qmc_run_and_check(
      short-H4-opt-adaptive
      "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
      H4-adaptive
      optm-adaptive.xml
      16
      1
      ${SUCCESS_STATUS_MP}
      5
      H4_OPT_SCALARS # OPT step 5
    )

    qmc_run_and_check(
      short-H4-opt-adaptive
      "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
      H4-adaptive
      optm-adaptive.xml
      1
      16
      ${SUCCESS_STATUS_MP}
      5
      H4_OPT_SCALARS # OPT step 5
    )

    #
    # H4 starting from perturbed orbitals, using optimizable determinants, and the adaptive linear method
    #
    list(APPEND H4_ORB_OPT_SCALARS "totenergy" "-2.0889 0.001") # total energy

    qmc_run_and_check(
      short-H4-RHF-orb-opt
      "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
      H4
      H4_RHF_orb_opt.xml
      16
      1
      ${SUCCESS_STATUS_MP}
      14
      H4_ORB_OPT_SCALARS # Final VMC step
      TRUE)

    qmc_run_and_check(
      short-H4-orb-opt
      "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
      H4
      H4_orb_opt.xml
      16
      1
      ${SUCCESS_STATUS_MP}
      14
      H4_ORB_OPT_SCALARS # Final VMC step
      TRUE)

    qmc_run_and_check(
      short-H4-orb-opt
      "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
      H4
      H4_orb_opt.xml
      1
      16
      ${SUCCESS_STATUS_MP}
      14
      H4_ORB_OPT_SCALARS # Final VMC step
      TRUE)

    #
    # H4 starting from converged wave function, doing a quick adaptive linear method run (using the same
    # spo set for both spin-up and spin-down determinants) and then doing a DMC calculation
    #
    list(APPEND H4_ORB_OPT_DMC_SCALARS "totenergy" "-2.191461 0.007") # total energy

    qmc_run_and_check(
      short-H4-orb-opt-dmc
      "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
      H4
      H4_orb_opt_dmc.xml
      16
      1
      ${SUCCESS_STATUS_MP}
      3
      H4_ORB_OPT_DMC_SCALARS # Test the DMC results
      TRUE)

    qmc_run_and_check(
      short-H4-orb-opt-dmc
      "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
      H4
      H4_orb_opt_dmc.xml
      1
      16
      ${SUCCESS_STATUS_MP}
      3
      H4_ORB_OPT_DMC_SCALARS # Test the DMC results
      TRUE)
    #
    # H4 starting from unoptimized CIS coefficients, targeting the first excited state, and using the adaptive linear method
    #

    list(APPEND H4_CJS_SCALARS "totenergy" "-2.0438 0.0045") # total energy

    qmc_run_and_check(
      short-H4-cjs
      "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
      H4
      H4_cjs_vmc.xml
      16
      1
      ${SUCCESS_STATUS_MP}
      0
      H4_CJS_SCALARS
      TRUE)

  else()
    message(
      VERBOSE
      "Skipping H4_ae adaptive optimizer tests because lmyengine interface was not built (BUILD_LMYENGINE_INTERFACE=0)")
  endif(BUILD_LMYENGINE_INTERFACE)

  list(APPEND H4_DESCENT_SCALARS "totenergy" "-2.047099 0.005")

  qmc_run_and_check(
    short-H4-descent
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-descent
    optm-descent.xml
    16
    1
    ${SUCCESS_STATUS_MP}
    9
    H4_DESCENT_SCALARS
    TRUE)

  qmc_run_and_check(
    short-H4-descent
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-descent
    optm-descent.xml
    1
    16
    ${SUCCESS_STATUS_MP}
    9
    H4_DESCENT_SCALARS
    TRUE)

  list(APPEND H4_TARGET_DESCENT_SCALARS "totenergy" "-2.04604154 0.003")

  qmc_run_and_check(
    short-H4-target-descent
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-target-descent
    optm-target-descent.xml
    16
    1
    ${SUCCESS_STATUS_MP}
    9
    H4_TARGET_DESCENT_SCALARS
    TRUE)

  qmc_run_and_check(
    short-H4-target-descent
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-target-descent
    optm-target-descent.xml
    1
    16
    ${SUCCESS_STATUS_MP}
    9
    H4_TARGET_DESCENT_SCALARS
    TRUE)

  list(APPEND H4_HYBRID_SCALARS "totenergy" "-2.140736 0.003")

  qmc_run_and_check(
    short-H4-hybrid
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-hybrid
    optm-hybrid.xml
    16
    1
    ${SUCCESS_STATUS_MP}
    9
    H4_HYBRID_SCALARS
    TRUE)

  qmc_run_and_check(
    short-H4-hybrid
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-hybrid
    optm-hybrid.xml
    1
    16
    ${SUCCESS_STATUS_MP}
    9
    H4_HYBRID_SCALARS
    TRUE)

  list(APPEND H4_BATCHED_HYBRID_SCALARS "totenergy" "-2.10456 0.01")

  qmc_run_and_check(
    short-H4-batched-hybrid
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-batched-hybrid
    optm-batched-hybrid.xml
    16
    1
    ${SUCCESS_STATUS_MP}
    9
    H4_BATCHED_HYBRID_SCALARS
    TRUE)

  qmc_run_and_check(
    short-H4-batched-hybrid
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-batched-hybrid
    optm-batched-hybrid.xml
    1
    16
    ${SUCCESS_STATUS_MP}
    9
    H4_BATCHED_HYBRID_SCALARS
    TRUE)

  list(APPEND H4_TARGET_HYBRID_SCALARS "totenergy" "-2.13527301 0.003")

  qmc_run_and_check(
    short-H4-target-hybrid
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-target-hybrid
    optm-target-hybrid.xml
    16
    1
    ${SUCCESS_STATUS_MP}
    9
    H4_TARGET_HYBRID_SCALARS
    TRUE)

  qmc_run_and_check(
    short-H4-target-hybrid
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-target-hybrid
    optm-target-hybrid.xml
    1
    16
    ${SUCCESS_STATUS_MP}
    9
    H4_TARGET_HYBRID_SCALARS
    TRUE)

  qmc_run_and_check(
    short-H4-opt-OneShiftOnly
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-OneShiftOnly
    optm-OneShiftOnly.xml
    16
    1
    ${SUCCESS_STATUS_MP}
    5
    H4_OPT_SCALARS # OPT step 5
  )

  qmc_run_and_check(
    short-H4-opt-OneShiftOnly
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-OneShiftOnly
    optm-OneShiftOnly.xml
    1
    16
    ${SUCCESS_STATUS_MP}
    5
    H4_OPT_SCALARS # OPT step 5
  )

  list(APPEND H4_OPT_SCALARS_ONLY_JAS "totenergy" "-2.066504 0.0009") # total energy
  qmc_run_and_check(
    short-H4-opt-OneShiftOnly-onlyjas
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-OneShiftOnly
    optm-OneShiftOnly-onlyjas.xml
    1
    16
    ${SUCCESS_STATUS_MP}
    5
    H4_OPT_SCALARS_ONLY_JAS # OPT step 5
  )

  list(APPEND H4_OPT_SCALARS_ONLY_MSD "totenergy" "-2.138635 0.0009") # total energy
  qmc_run_and_check(
    short-H4-opt-OneShiftOnly-onlymsd
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-OneShiftOnly
    optm-OneShiftOnly-onlymsd.xml
    1
    16
    ${SUCCESS_STATUS_MP}
    5
    H4_OPT_SCALARS_ONLY_MSD # OPT step 5
  )

  qmc_run_and_check(
    short-H4-optbatch-OneShiftOnly-multistep
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-OneShiftOnly
    optm-OneShiftOnly-multistep.xml
    1
    16
    ${SUCCESS_STATUS_MP}
    11
    H4_OPT_SCALARS_ONLY_MSD # OPT step 011
  )

  qmc_run_and_check(
    deterministic-input_check-invalid_variational_subset
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-OneShiftOnly
    det-input-check-optm-invalid-subset.xml
    2
    1
    FALSE)
else()
  message(VERBOSE "Skipping H4_ae tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
endif()
