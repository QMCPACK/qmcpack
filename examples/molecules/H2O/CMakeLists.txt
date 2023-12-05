# Test example

if(NOT QMC_COMPLEX)

  list(APPEND H2O_SCALARS "totenergy" "-17.26138 .006") # total energy

  qmc_run_and_check(
    example_H2O
    "${qmcpack_SOURCE_DIR}/examples/molecules/H2O"
    H2O
    simple-H2O.xml
    1
    1
    TRUE
    2
    H2O_SCALARS # check the DMC data
  )
endif()
