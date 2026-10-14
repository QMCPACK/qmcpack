if(NOT QMC_COMPLEX)

  # 4 H atom cluster gaussian tests, all electron
  #
  # H4 Multi-Slater Jastrow wfn. Fixed CI and J parameters but the orbitals are allowed to relax
  #
  list(APPEND H4_MSD_OPT_ORB_SCALARS "totenergy" "-2.160310 0.001") # total energy

  qmc_run_and_check(
    short-H4-msd-opt_orb
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_msd_ae"
    H4
    H4_msd_opt_orb.xml
    16
    1
    ${SUCCESS_STATUS_MP}
    4
    H4_MSD_OPT_ORB_SCALARS # OPT step 4
  )
  qmc_run_and_check(
    short-H4-msd-optbatch_orb
    "${qmcpack_SOURCE_DIR}/tests/molecules/H4_msd_ae"
    H4
    H4_msd_optbatch_orb.xml
    16
    1
    ${SUCCESS_STATUS_MP}
    4
    H4_MSD_OPT_ORB_SCALARS # OPT step 4
  )

else()
  message(VERBOSE
          "Skipping H4_msd_ae tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
endif()

get_property(CURRENT_TESTS DIRECTORY PROPERTY TESTS)
add_labels_for_tests("${CURRENT_TESTS}")
