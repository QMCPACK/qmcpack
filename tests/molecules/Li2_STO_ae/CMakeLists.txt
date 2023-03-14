if(NOT QMC_COMPLEX)
  #
  # Li2 molecular dimer STO tests, all electron
  #   Reproduces UNR 1993
  #

  # Short tests
  list(APPEND LI2STOSJS_VMC_SCALARS "totenergy" "-14.947214  0.000502") # total energy
  list(APPEND LI2STOSJS_VMC_SCALARS "variance" "  0.170618  0.001035") # energy variance
  list(APPEND LI2STOSJS_DMC_SCALARS "totenergy" "-14.989438  0.000643") # total energy
  list(APPEND LI2STOSJS_DMC_SCALARS "variance" "  0.151118  0.000312") # energy variance

  qmc_run_and_check(
    short-li2_sto-sj_dmc
    "${qmcpack_SOURCE_DIR}/tests/molecules/Li2_STO_ae"
    Li2.STO.short
    Li2.STO.short.in.xml
    1
    16
    TRUE
    0
    LI2STOSJS_VMC_SCALARS # VMC
    2
    LI2STOSJS_DMC_SCALARS # DMC
  )

  qmc_run_and_check(
    short-li2_sto-sj_dmc
    "${qmcpack_SOURCE_DIR}/tests/molecules/Li2_STO_ae"
    Li2.STO.short
    Li2.STO.short.in.xml
    4
    4
    TRUE
    0
    LI2STOSJS_VMC_SCALARS # VMC
    2
    LI2STOSJS_DMC_SCALARS # DMC
  )

  # Long tests
  list(APPEND LI2STOSJL_VMC_SCALARS "totenergy" "-14.947214  0.000166") # total energy
  list(APPEND LI2STOSJL_VMC_SCALARS "variance" "  0.170618  0.000342") # energy variance
  list(APPEND LI2STOSJL_DMC_SCALARS "totenergy" "-14.989438  0.000212") # total energy
  list(APPEND LI2STOSJL_DMC_SCALARS "variance" "  0.151118  0.000103") # energy variance

  qmc_run_and_check(
    long-li2_sto-sj_dmc
    "${qmcpack_SOURCE_DIR}/tests/molecules/Li2_STO_ae"
    Li2.STO.long
    Li2.STO.long.in.xml
    1
    16
    TRUE
    0
    LI2STOSJL_VMC_SCALARS # VMC
    2
    LI2STOSJL_DMC_SCALARS # DMC
  )

else()
  message(VERBOSE "Skipping Li2_STO_ae tests because the basis sets are not supported by complex build (QMC_COMPLEX=1)")
endif()
