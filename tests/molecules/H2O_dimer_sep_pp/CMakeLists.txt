if(NOT QMC_COMPLEX)

  # 2 water molecules far away test for Tmoves.
  # short tests
  list(APPEND H2OS2_DMCLA_SCALARS "totenergy" "-34.550288  0.0015") # total energy
  list(APPEND H2OS2_DMCLA_SCALARS "nonlocalecp" "3.7922 0.012") # non local ecp

  qmc_run_and_check(
    short-H2O_dimer_sep_pp-j3_dmc_la
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    H2O
    short-H2O-LA.xml
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
    H2O
    short-H2O-TM0.xml
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
    H2O
    short-H2O-TM1.xml
    1
    16
    TRUE
    2
    H2OS2_DMCTM1_SCALARS)

  list(APPEND H2OS2_DMCTM3_SCALARS "totenergy" "-34.560980  0.0015") # total energy
  list(APPEND H2OS2_DMCTM3_SCALARS "nonlocalecp" "3.174425 0.012") # non local ecp

  qmc_run_and_check(
    short-H2O_dimer_sep_pp-j3_dmc_tm3
    "${qmcpack_SOURCE_DIR}/tests/molecules/H2O_dimer_sep_pp"
    H2O
    short-H2O-TM3.xml
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
    H2O
    long-H2O-LA.xml
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
    H2O
    long-H2O-TM0.xml
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
    H2O
    long-H2O-TM1.xml
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
    H2O
    long-H2O-TM3.xml
    1
    16
    TRUE
    2
    LONG_H2OS2_DMCTM3_SCALARS)

else()
  message(
    VERBOSE
    "Skipping H2O_dimer_sep_pp tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
endif()
