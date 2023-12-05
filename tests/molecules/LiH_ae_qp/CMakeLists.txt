if(NOT QMC_COMPLEX)
  #
  # LiH molecule gaussian tests, all electron using QP
  # Also check results for different number of mpi tasks and threads keeping total constant
  # Energy from Quantum Package: E=-7.987323645629950
  #
  #   "kinetic" "7.991334  0.000065") # kinetic energy
  #   "totenergy" "-7.9873140 0.0000086 ") # total energy
  #   "eeenergy" "3.488827 0.000011") # e-e energy
  #   "samples" "79626240000 0.0") # samples

  list(APPEND LIH_SHORT "kinetic" "7.991334 0.01") # kinetic energy
  list(APPEND LIH_SHORT "totenergy" "-7.9873140 0.0009") # total energy
  list(APPEND LIH_SHORT "eeenergy" "3.488827 0.00088") # e-e energy
  list(APPEND LIH_SHORT "samples" "9600000 0.0") # samples
  qmc_run_and_check(
    short-LiH_ae_qp-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_qp"
    vmc_short_noj
    vmc_short_noj.in.xml
    1
    16
    TRUE
    0
    LIH_SHORT # VMC
  )

  list(APPEND LIH_LONG "kinetic" "7.991334 0.01") # kinetic energy
  list(APPEND LIH_LONG "totenergy" "-7.9873140 0.0002") # total energy
  list(APPEND LIH_LONG "eeenergy" "3.488827 0.00026") # e-e energy
  list(APPEND LIH_LONG "samples" "96000000 0.0") # samples
  qmc_run_and_check(
    long-LiH_ae_qp-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_qp"
    vmc_long_noj
    vmc_long_noj.in.xml
    1
    16
    TRUE
    0
    LIH_LONG # VMC
  )

else()
  message(VERBOSE
          "Skipping LiH_ae tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
endif()
