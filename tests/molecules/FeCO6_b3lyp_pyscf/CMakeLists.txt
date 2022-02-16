if(NOT QMC_CUDA)
  if(NOT QMC_COMPLEX)
    #
    # [Fe(CO)6]2+ molecule gaussian tests, with BFD ECP using pyscf

    # Ref for ref VMC test
    #   "kinetic" "163.8420 0.0022") # kinetic energy
    #   "totenergy" "-250.43499 0.00034") # total energy
    #   "eeenergy" "568.6012 0.0017") # e-e energy
    #   "samples" "655360000 0.0") # samples

    list(APPEND FeCO6_PYSCF_SHORT "kinetic" "163.8420 0.27") # kinetic energy
    list(APPEND FeCO6_PYSCF_SHORT "eeenergy" "568.6012 0.22") # e-e energy
    list(APPEND FeCO6_PYSCF_SHORT "totenergy" "-250.43499 0.035") # total energy
    list(APPEND FeCO6_PYSCF_SHORT "samples" "48000 0.0") # samples
    qmc_run_and_check(
      short-FeCO6_pyscf-vmc_noj
      "${qmcpack_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_pyscf"
      vmc_short_noj
      vmc_short_noj.in.xml
      1
      16
      TRUE
      0
      FeCO6_PYSCF_SHORT # VMC
    )

    list(APPEND FeCO6_PYSCF_LONG "kinetic" "163.8420 0.076") # kinetic energy
    list(APPEND FeCO6_PYSCF_LONG "eeenergy" "568.6012 0.063") # e-e energy
    list(APPEND FeCO6_PYSCF_LONG "totenergy" "-250.43499 0.011") # total energy
    list(APPEND FeCO6_PYSCF_LONG "samples" "480000 0.0") # samples
    qmc_run_and_check(
      long-FeCO6_pyscf-vmc_noj
      "${qmcpack_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_pyscf"
      vmc_long_noj
      vmc_long_noj.in.xml
      1
      16
      TRUE
      0
      FeCO6_PYSCF_LONG # VMC
    )

if(QMC_MIXED_PRECISION)
    list(APPEND DET_FeCO6_PYSCF_SHORT "kinetic" "205.29692267 0.003") # kinetic energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "eeenergy" "592.73730675 0.001") # e-e energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "totenergy" "-244.72443901 0.001") # total energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "potential" "-450.02136168 0.001") # potential energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "ionion" "422.43049005 0.0001") # i-i energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "localecp" "-1542.74871820 0.003") # localECP
    list(APPEND DET_FeCO6_PYSCF_SHORT "nonlocalecp" "77.55955951 0.002") # nonlocalECP
    list(APPEND DET_FeCO6_PYSCF_SHORT "samples" "9 0.0") # samples
else()
    list(APPEND DET_FeCO6_PYSCF_SHORT "kinetic" "157.44334976 0.000001") # kinetic energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "eeenergy" "564.18503814 0.000001") # e-e energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "totenergy" "-249.23219138 0.000001") # total energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "potential" "-406.67554113 0.000001") # potential energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "ionion" "422.43049172 0.000001") # i-i energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "localecp" "-1420.56055990 0.000001") # localECP
    list(APPEND DET_FeCO6_PYSCF_SHORT "nonlocalecp" "27.26948891 0.000001") # nonlocalECP
    list(APPEND DET_FeCO6_PYSCF_SHORT "samples" "9 0.0") # samples
endif()

    qmc_run_and_check(
      deterministic-FeCO6_pyscf-vmc_noj
      "${qmcpack_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_pyscf"
      det_vmc_noj
      det_vmc_noj.in.xml
      1
      1
      TRUE
      0
      DET_FeCO6_PYSCF_SHORT # VMC
    )


if(QMC_MIXED_PRECISION)
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "kinetic" "205.29692267 0.003") # kinetic energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "eeenergy" "592.73730675 0.001") # e-e energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "totenergy" "-244.72443644 0.001") # total energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "potential" "-450.02136168 0.001") # potential energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "ionion" "422.43049027 0.0001") # i-i energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "localecp" "-1542.74871820 0.003") # localECP
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "nonlocalecp" "77.55955951 0.002") # nonlocalECP
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "samples" "9 0.0") # samples
else()
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "kinetic" "157.44334976 0.000001") # kinetic energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "eeenergy" "564.18503814 0.000001") # e-e energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "totenergy" "-249.23219138 0.000001") # total energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "potential" "-406.67554113 0.000001") # potential energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "ionion" "422.43049172 0.000001") # i-i energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "localecp" "-1420.56055990 0.000001") # localECP
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "nonlocalecp" "27.26948891 0.000001") # nonlocalECP
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "samples" "9 0.0") # samples
endif()

    qmc_run_and_check(
      deterministic-FeCO6_pyscf-vmcbatch_noj
      "${qmcpack_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_pyscf"
      det_vmcbatch_noj
      det_vmcbatch_noj.in.xml
      1
      1
      TRUE
      0
      DET_FeCO6_PYSCF_SHORT_BATCHED # VMC
    )

  else()
    message(VERBOSE
      "Skipping  FeCO6_b3lyp_gms tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
  endif()
else()
  message(VERBOSE
    "Skipping FeCO6_b3lyp_gms tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
endif()
