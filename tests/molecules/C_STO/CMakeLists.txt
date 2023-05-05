if(NOT QMC_COMPLEX)

  # Test optimization of orbital rotation using carbon atom with STO wavefunction
  #  The error bounds are based on the spread of values from many runs of the optimizer,
  #   not a single run

  list(APPEND C_STO_OPT_SCALARS "totenergy" "-37.669 0.12") # total energy

  qmc_run_and_check(
    short-C_STO-opt-linear
    "${qmcpack_SOURCE_DIR}/tests/molecules/C_STO"
    carbon_sto
    opt.in.xml
    1
    16
    ${SUCCESS_STATUS_MP}
    3
    C_STO_OPT_SCALARS # OPT step 3
  )
else()
  message(VERBOSE "Skipping C_STO tests because linear optimizer is not supported by the complex build")
endif()

