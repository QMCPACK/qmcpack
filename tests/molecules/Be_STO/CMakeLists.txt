if(NOT QMC_COMPLEX)

  # From Bunge, Barrientos, and Bunge, Atomic Data and Nuclear Data Tables 53, 113-162 (1993)
  # Restricted HF results using STO basis functions
  list(APPEND BE_STO_SHORT LocalEnergy -14.57302313 0.0075)
  list(APPEND BE_STO_SHORT Kinetic 14.57302427 0.23)
  list(APPEND BE_STO_SHORT LocalPotential -29.14604740 0.23)

  qmc_run_and_check_custom_scalar(
    BASE_NAME
    short-Be_STO-vmc
    BASE_DIR
    "${qmcpack_SOURCE_DIR}/tests/molecules/Be_STO"
    PREFIX
    be_vmc
    INPUT_FILE
    be_sto_vmc.in.xml
    PROCS
    1
    THREADS
    16
    SERIES
    0
    SCALAR_VALUES
    BE_STO_SHORT)

  qmc_run_and_check_custom_scalar(
    BASE_NAME
    short-Be_STO-vmc-nodrift
    BASE_DIR
    "${qmcpack_SOURCE_DIR}/tests/molecules/Be_STO"
    PREFIX
    be_vmc_nodrift
    INPUT_FILE
    be_sto_nodrift_vmc.in.xml
    PROCS
    1
    THREADS
    16
    SERIES
    0
    SCALAR_VALUES
    BE_STO_SHORT)

  qmc_run_and_check_custom_scalar(
    BASE_NAME
    short-Be_STO-vmc-legacy
    BASE_DIR
    "${qmcpack_SOURCE_DIR}/tests/molecules/Be_STO"
    PREFIX
    be_vmc_legacy
    INPUT_FILE
    be_sto_vmc_legacy.in.xml
    PROCS
    1
    THREADS
    16
    SERIES
    0
    SCALAR_VALUES
    BE_STO_SHORT)


  if(NOT QMC_MIXED_PRECISION)
    # Test orbital rotation optimization.  Energy is HF value.  Error comes from
    #  visual inspection of an ensemble of runs.
    list(APPEND BE_STO_OPT_SHORT LocalEnergy -14.57302313 0.05)

    qmc_run_and_check_custom_scalar(
      BASE_NAME
      short-Be_STO-orb_opt
      BASE_DIR
      "${qmcpack_SOURCE_DIR}/tests/molecules/Be_STO"
      PREFIX
      be_opt
      INPUT_FILE
      be_opt.in.xml
      PROCS
      1
      THREADS
      16
      SERIES
      4
      SCALAR_VALUES
      BE_STO_OPT_SHORT)
  else()
    message(VERBOSE "Skipping Be_STO orb_opt tests because linear optimizer is not supported by mixed precision builds (QMC_MIXED_PRECISION=1)")
  endif()
else()
  message(VERBOSE "Skipping Be_STO tests because STO basis sets are not supported by complex build (QMC_COMPLEX=1)")
endif()
