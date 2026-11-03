if(NOT QMC_COMPLEX)

  # 2 water molecules far away test for Tmoves.
  # short tests
  list(APPEND H2OS2_DMCLA_SCALARS "totenergy" "-34.550288  0.0015") # total energy
  list(APPEND H2OS2_DMCLA_SCALARS "nonlocalecp" "3.7922 0.012") # non local ecp

  qmc_run_and_check(
    short-H2O_dimer_sep_pp-j3_dmc_la
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    short-H2O-LA.xml
    1
    16
    TRUE
    2
    H2OS2_DMCLA_SCALARS)

  qmc_run_and_check(
    short-H2O_dimer_sep_pp-j3_dmcbatch_la
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    short-H2O-LA-batch.xml
    1
    16
    TRUE
    2
    H2OS2_DMCLA_SCALARS)

  list(APPEND H2OS2_DMCTM0_SCALARS "totenergy" "-34.556348  0.0015") # total energy
  list(APPEND H2OS2_DMCTM0_SCALARS "nonlocalecp" "3.339723 0.012") # non local ecp

  qmc_run_and_check(
    short-H2O_dimer_sep_pp-j3_dmc_tm0
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    short-H2O-TM0.xml
    1
    16
    TRUE
    2
    H2OS2_DMCTM0_SCALARS)

  qmc_run_and_check(
    short-H2O_dimer_sep_pp-j3_dmcbatch_tm0
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    short-H2O-TM0-batch.xml
    1
    16
    TRUE
    2
    H2OS2_DMCTM0_SCALARS)

  list(APPEND H2OS2_DMCTM1_SCALARS "totenergy" "-34.560303  0.0015") # total energy
  list(APPEND H2OS2_DMCTM1_SCALARS "nonlocalecp" "3.177908 0.012") # non local ecp

  qmc_run_and_check(
    short-H2O_dimer_sep_pp-j3_dmc_tm1
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    short-H2O-TM1.xml
    1
    16
    TRUE
    2
    H2OS2_DMCTM1_SCALARS)

  qmc_run_and_check(
    short-H2O_dimer_sep_pp-j3_dmcbatch_tm1
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    short-H2O-TM1-batch.xml
    1
    16
    TRUE
    2
    H2OS2_DMCTM1_SCALARS)

  list(APPEND H2OS2_DMCTM3_SCALARS "totenergy" "-34.560980  0.0015") # total energy
  list(APPEND H2OS2_DMCTM3_SCALARS "nonlocalecp" "3.174425 0.018") # non local ecp

  qmc_run_and_check(
    short-H2O_dimer_sep_pp-j3_dmc_tm3
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    short-H2O-TM3.xml
    1
    16
    TRUE
    2
    H2OS2_DMCTM3_SCALARS)

  qmc_run_and_check(
    short-H2O_dimer_sep_pp-j3_dmcbatch_tm3
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    short-H2O-TM3-batch.xml
    1
    16
    TRUE
    2
    H2OS2_DMCTM3_SCALARS)

  # long tests
  list(APPEND LONG_H2OS2_DMCLA_SCALARS "totenergy" "-34.550288  0.0005") # total energy
  list(APPEND LONG_H2OS2_DMCLA_SCALARS "nonlocalecp" "3.7922 0.0045") # non local ecp

  qmc_run_and_check(
    long-H2O_dimer_sep_pp-j3_dmc_la
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    long-H2O-LA.xml
    1
    16
    TRUE
    2
    LONG_H2OS2_DMCLA_SCALARS)

  qmc_run_and_check(
    long-H2O_dimer_sep_pp-j3_dmcbatch_la
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    long-H2O-LA-batch.xml
    1
    16
    TRUE
    2
    LONG_H2OS2_DMCLA_SCALARS)

  list(APPEND LONG_H2OS2_DMCTM0_SCALARS "totenergy" "-34.556348  0.0005") # total energy
  list(APPEND LONG_H2OS2_DMCTM0_SCALARS "nonlocalecp" "3.339723 0.0045") # non local ecp

  qmc_run_and_check(
    long-H2O_dimer_sep_pp-j3_dmc_tm0
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    long-H2O-TM0.xml
    1
    16
    TRUE
    2
    LONG_H2OS2_DMCTM0_SCALARS)

  qmc_run_and_check(
    long-H2O_dimer_sep_pp-j3_dmcbatch_tm0
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    long-H2O-TM0-batch.xml
    1
    16
    TRUE
    2
    LONG_H2OS2_DMCTM0_SCALARS)

  list(APPEND LONG_H2OS2_DMCTM1_SCALARS "totenergy" "-34.560303  0.0005") # total energy
  list(APPEND LONG_H2OS2_DMCTM1_SCALARS "nonlocalecp" "3.177908 0.0045") # non local ecp

  qmc_run_and_check(
    long-H2O_dimer_sep_pp-j3_dmc_tm1
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    long-H2O-TM1.xml
    1
    16
    TRUE
    2
    LONG_H2OS2_DMCTM1_SCALARS)

  qmc_run_and_check(
    long-H2O_dimer_sep_pp-j3_dmcbatch_tm1
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    long-H2O-TM1-batch.xml
    1
    16
    TRUE
    2
    LONG_H2OS2_DMCTM1_SCALARS)

  list(APPEND LONG_H2OS2_DMCTM3_SCALARS "totenergy" "-34.560980  0.0005") # total energy
  list(APPEND LONG_H2OS2_DMCTM3_SCALARS "nonlocalecp" "3.174425 0.0045") # non local ecp

  qmc_run_and_check(
    long-H2O_dimer_sep_pp-j3_dmc_tm3
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    long-H2O-TM3.xml
    1
    16
    TRUE
    2
    LONG_H2OS2_DMCTM3_SCALARS)

  qmc_run_and_check(
    long-H2O_dimer_sep_pp-j3_dmcbatch_tm3
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    long-H2O-TM3-batch.xml
    1
    16
    TRUE
    2
    LONG_H2OS2_DMCTM3_SCALARS)

  foreach(LENGTH short long)
    foreach(MOVE LA TM1)
      string(TOLOWER "${MOVE}" MOVE_SUFFIX)
      set(SCALARS H2OS2_DMC${MOVE}_SCALARS)
      if(LENGTH STREQUAL "long")
        set(SCALARS LONG_${SCALARS})
      endif()
      qmc_run_and_check(
        ${LENGTH}-H2O_dimer_sep_pp-j3_dmcbatch_${MOVE_SUFFIX}_reconfig
        "${CMAKE_CURRENT_SOURCE_DIR}"
        H2O.reconfig
        ${LENGTH}-H2O-${MOVE}-batch-reconfig.xml
        1
        16
        TRUE
        2 ${SCALARS})
      add_test_check_fixed_population(${LENGTH}-H2O_dimer_sep_pp-j3_dmcbatch_${MOVE_SUFFIX}_reconfig-r1-t16
                                     H2O.reconfig.s002.dmc.dat 512)
    endforeach()
  endforeach()

  foreach(MODE variable reconfig)
    set(PROFILE_INPUT short-H2O-LA-batch.xml)
    if(MODE STREQUAL "reconfig")
      set(PROFILE_INPUT short-H2O-LA-batch-reconfig.xml)
    endif()
    set(PROFILE_TEST performance-H2O_dimer_sep_pp-dmcbatch_la_${MODE}-r4-t4)
    set(PROFILE_DIR "${CMAKE_CURRENT_BINARY_DIR}/${PROFILE_TEST}")
    copy_directory_maybe_using_symlink("${CMAKE_CURRENT_SOURCE_DIR}" "${PROFILE_DIR}" "${PROFILE_INPUT}")
    run_qmc_app_no_copy(
      ${PROFILE_TEST}
      "${PROFILE_DIR}"
      4
      4
      PROFILE_ADDED
      --enable-timers=medium ${PROFILE_INPUT})
    if(PROFILE_ADDED)
      set_property(TEST ${PROFILE_TEST} APPEND PROPERTY ENVIRONMENT
                   OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1)
    endif()
    if(MODE STREQUAL "reconfig")
      add_test_check_fixed_population(${PROFILE_TEST} H2O.reconfig.s002.dmc.dat 512)
    endif()
  endforeach()

else()
  message(
    VERBOSE
    "Skipping H2O_dimer_sep_pp tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
endif()

get_property(CURRENT_TESTS DIRECTORY PROPERTY TESTS)
add_labels_for_tests("${CURRENT_TESTS}")
