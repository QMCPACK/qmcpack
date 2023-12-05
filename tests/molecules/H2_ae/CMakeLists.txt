if(NOT QMC_COMPLEX)
  if(BUILD_LMYENGINE_INTERFACE)

    #
    # H2 test in a DZV basis set, using optimizable determinants, and the adaptive linear method
    #
    list(APPEND H2_OPT_SCALARS "totenergy" "-1.124563 0.00084") # total energy

    qmc_run_and_check(
      short-H2-orb-opt
      "${qmcpack_SOURCE_DIR}/tests/molecules/H2_ae"
      H2
      h2_orb_opt.xml
      16
      1
      ${SUCCESS_STATUS_MP}
      10
      H2_OPT_SCALARS # Final VMC step
      TRUE)

    qmc_run_and_check(
      short-H2-orb-opt-spindepjas
      "${qmcpack_SOURCE_DIR}/tests/molecules/H2_ae"
      H2
      h2_orb_opt_spindepjas.xml
      16
      1
      ${SUCCESS_STATUS_MP}
      10
      H2_OPT_SCALARS # Final VMC step
      TRUE)

  else()
    message(VERBOSE "Skipping H2 test because lmyengine interface was not built (BUILD_LMYENGINE_INTERFACE=0)")
  endif(BUILD_LMYENGINE_INTERFACE)
else()
  message(VERBOSE "Skipping H2 tests because orbital rotation is not supported by complex build (QMC_COMPLEX=1)")
endif()
