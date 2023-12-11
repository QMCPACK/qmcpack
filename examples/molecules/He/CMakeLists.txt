# Test example

if(MIXED_PRECISION)
  set(MP_SUCCESS FALSE)
else()
  set(MP_SUCCESS TRUE)
endif()

if(NOT QMC_COMPLEX)

  list(APPEND HE_SIMPLE_VMC_SCALARS "totenergy" "-2.83 .006") # total energy

  qmc_run_and_check(
    example_He_simple
    "${qmcpack_SOURCE_DIR}/examples/molecules/He"
    He
    he_simple.xml
    1
    1
    TRUE
    0
    HE_SIMPLE_VMC_SCALARS)

  list(APPEND HE_SIMPLE_DMC_SCALARS "totenergy" "-2.902 .002") # total energy

  qmc_run_and_check(
    example_He_simple_dmc
    "${qmcpack_SOURCE_DIR}/examples/molecules/He"
    He
    he_simple_dmc.xml
    1
    1
    TRUE
    1
    HE_SIMPLE_DMC_SCALARS # series for DMC data
  )

  list(APPEND HE_SIMPLE_OPT_SCALARS "totenergy" "-2.88 .005") # total energy

  qmc_run_and_check(
    example_He_simple_opt
    "${qmcpack_SOURCE_DIR}/examples/molecules/He"
    He
    he_simple_opt.xml
    1
    1
    ${MP_SUCCESS}
    10
    HE_SIMPLE_OPT_SCALARS # final optimized series
  )

  qmc_run_and_check(
    example_He_example_wf
    "${qmcpack_SOURCE_DIR}/examples/molecules/He"
    He
    he_example_wf.xml
    1
    1
    ${MP_SUCCESS}
    10
    HE_SIMPLE_OPT_SCALARS # final optimized series
  )

  list(APPEND HE_FROM_GAMESS_SCALARS "totenergy" "-2.872 .003") # total energy

  qmc_run_and_check(
    example_He_from_gamess
    "${qmcpack_SOURCE_DIR}/examples/molecules/He"
    He
    he_from_gamess.xml
    1
    1
    TRUE
    0
    HE_FROM_GAMESS_SCALARS)

  # Not included for now because it fails spectacularly a small percentage of
  # the time (2/100 runs).  Values at each optimization stage seem to be okay,
  # but the final VMC is far off.
  #  LIST(APPEND HE_BSPLINE_OPT_SCALARS "totenergy" "-2.879 .003") # total energy
  #
  #  QMC_RUN_AND_CHECK(example_He_bspline_opt
  #                    "${qmcpack_SOURCE_DIR}/examples/molecules/He"
  #                    He
  #                    he_bspline_jastrow.xml
  #                    1 1
  #                    TRUE
  #                    10 HE_BSPLINE_OPT_SCALARS # series for VMC on optimized value
  #                    )
endif()
