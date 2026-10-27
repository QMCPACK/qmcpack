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
    short-H4-optbatch-linear
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-linear
    optm-batch-linear.xml
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
    short-H4-optbatch-linear-rescale
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-linear-rescale
    optm-batch-linear-rescale.xml
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
  qmc_run_and_check(
    short-H4-optbatch-linear-linemin
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-linear-linemin
    optm-batch-linear-linemin.xml
    16
    1
    ${SUCCESS_STATUS_MP}
    5
    H4_OPT_SCALARS # OPT step 5
  )


  if(NOT MIXED_PRECISION)
    check_python_reqs(h5py deterministic-H4-descent-vp ADD_DESCENT_VP_TESTS)

    set(H4_DETERMINISTIC_BATCHED_DESCENT_TEST deterministic-H4-batched-descent-r1-t4)
    set(H4_DETERMINISTIC_BATCHED_DESCENT_TEST_ADDED FALSE)
    run_qmc_app(
      ${H4_DETERMINISTIC_BATCHED_DESCENT_TEST}
      "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
      1
      4
      H4_DETERMINISTIC_BATCHED_DESCENT_TEST_ADDED
      det-optm-batched-descent.xml)

    if(H4_DETERMINISTIC_BATCHED_DESCENT_TEST_ADDED)
      add_test(
        NAME ${H4_DETERMINISTIC_BATCHED_DESCENT_TEST}-energy
        COMMAND
          ${Python3_EXECUTABLE} ${qmcpack_SOURCE_DIR}/tests/scripts/check_scalars.py --ns 3 --series 1 -p
          det-H4-batched-descent -e 0 --le "-1.2857706733 0.000001"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${H4_DETERMINISTIC_BATCHED_DESCENT_TEST}")

      add_test(
        NAME ${H4_DETERMINISTIC_BATCHED_DESCENT_TEST}-coefficients
        COMMAND
          ${Python3_EXECUTABLE} ${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae/check_h4_descent.py
          det-H4-batched-descent.s000.opt.xml reference.det-H4-batched-descent.s000.opt.xml
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${H4_DETERMINISTIC_BATCHED_DESCENT_TEST}")

      set_tests_properties(
        ${H4_DETERMINISTIC_BATCHED_DESCENT_TEST}-energy ${H4_DETERMINISTIC_BATCHED_DESCENT_TEST}-coefficients
        PROPERTIES DEPENDS ${H4_DETERMINISTIC_BATCHED_DESCENT_TEST} LABELS
                   "QMCPACK-checking-results;deterministic")

      if(ADD_DESCENT_VP_TESTS)
        add_test(
          NAME ${H4_DETERMINISTIC_BATCHED_DESCENT_TEST}-vp-h5
          COMMAND
            ${Python3_EXECUTABLE} ${qmcpack_SOURCE_DIR}/tests/scripts/check_vp_h5.py
            det-H4-batched-descent.s000.vp.h5 det-H4-batched-descent.s000.opt.xml
          WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${H4_DETERMINISTIC_BATCHED_DESCENT_TEST}")
        set_tests_properties(
          ${H4_DETERMINISTIC_BATCHED_DESCENT_TEST}-vp-h5
          PROPERTIES DEPENDS ${H4_DETERMINISTIC_BATCHED_DESCENT_TEST} LABELS
                     "QMCPACK-checking-results;deterministic")
      endif()
    endif()
  endif()

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

  qmc_run_and_check(
    short-H4-optbatch-OneShiftOnly
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-OneShiftOnly
    optm-batch-OneShiftOnly.xml
    1
    16
    ${SUCCESS_STATUS_MP}
    5
    H4_OPT_SCALARS # OPT step 5
  )

  qmc_run_and_check(
    short-H4-optbatch-OneShiftOnly
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-OneShiftOnly
    optm-batch-OneShiftOnly.xml
    2
    4
    ${SUCCESS_STATUS_MP}
    5
    H4_OPT_SCALARS # OPT step 5
  )

  qmc_run_and_check(
    short-H4-optbatch-OneShiftOnly
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-OneShiftOnly
    optm-batch-OneShiftOnly.xml
    3
    2
    ${SUCCESS_STATUS_MP}
    5
    H4_OPT_SCALARS # OPT step 5
  )

  qmc_run_and_check(
    short-H4-optbatch-OneShiftOnly-general
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-OneShiftOnly-general
    optm-batch-OneShiftOnly-general.xml
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
  qmc_run_and_check(
    short-H4-optbatch-OneShiftOnly-onlyjas
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-OneShiftOnly
    optm-batch-OneShiftOnly-onlyjas.xml
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
    short-H4-optbatch-OneShiftOnly-onlymsd
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    H4-OneShiftOnly
    optm-batch-OneShiftOnly-onlymsd.xml
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

  run_qmc_app(
    deterministic-check-eigenvalue-choice
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_ae"
    1
    1
    TEST_ADDED_EIGENVALUE_CHOICE
    det-check-eigenvalue-choice.xml
    )
  # Check that the general eigenvalue solver was chosen. Overwrite the usual passing regular expression.
  set_property(TEST deterministic-check-eigenvalue-choice PROPERTY PASS_REGULAR_EXPRESSION "Using generalized eigenvalue solver")
  if (NOT ${SUCCESS_STATUS_MP})
    set_property(TEST deterministic-check-eigenvalue-choice PROPERTY WILL_FAIL TRUE)
  endif()

else()
  message(VERBOSE "Skipping H4_ae tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
endif()

get_property(CURRENT_TESTS DIRECTORY PROPERTY TESTS)
add_labels_for_tests("${CURRENT_TESTS}")
