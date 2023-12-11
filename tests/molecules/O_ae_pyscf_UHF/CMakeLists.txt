if(NOT QMC_COMPLEX)
  #
  # O triplet UHF molecule gaussian tests, all electron using pyscf
  # Also check results for different number of mpi tasks and threads keeping total constant
  # Energy from Pyscf: E=-74.7923237196089
  #
  #   "kinetic" "74.7451 0.0018") # kinetic energy
  #   "totenergy" "-74.793488 0.000038 ") # total energy
  #   "eeenergy" "28.483604 0.000065") # e-e energy
  #   "samples" "98304000000 0.0") # samples

  list(APPEND O_UHF_SHORT "kinetic" "74.7451  0.39") # kinetic energy
  list(APPEND O_UHF_SHORT "totenergy" "-74.793488 0.0075") # total energy
  list(APPEND O_UHF_SHORT "eeenergy" "28.483604 0.011") # e-e energy
  list(APPEND O_UHF_SHORT "samples" "9600000 0.0") # samples
  qmc_run_and_check(
    short-O_ae_uhf_pyscf-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/O_ae_pyscf_UHF"
    vmc_short_noj
    vmc_short_noj.in.xml
    1
    16
    TRUE
    0
    O_UHF_SHORT # VMC
  )

  list(APPEND O_UHF_LONG "kinetic" "74.7451 0.057") # kinetic energy
  list(APPEND O_UHF_LONG "totenergy" "-74.793488 0.0010") # total energy
  list(APPEND O_UHF_LONG "eeenergy" "28.483604 0.0018") # e-e energy
  list(APPEND O_UHF_LONG "samples" "96000000 0.0") # samples
  qmc_run_and_check(
    long-O_ae_uhf_pyscf-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/O_ae_pyscf_UHF"
    vmc_long_noj
    vmc_long_noj.in.xml
    1
    16
    TRUE
    0
    O_UHF_LONG # VMC
  )

else()
  message(VERBOSE
          "Skipping O_ae_uhf tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
endif()
