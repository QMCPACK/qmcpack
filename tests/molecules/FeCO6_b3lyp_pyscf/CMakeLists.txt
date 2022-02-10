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
    list(APPEND DET_FeCO6_PYSCF_SHORT "kinetic" "158.69060050 0.00001") # kinetic energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "eeenergy" "567.89193761 0.00001") # e-e energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "totenergy" "-250.78482579 0.00001") # total energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "potential" "-409.47542629 0.00001") # potential energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "ionion" "422.43049027 0.00001") # i-i energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "localecp" "-1436.69967390 0.00001") # localECP
    list(APPEND DET_FeCO6_PYSCF_SHORT "nonlocalecp" "36.90181976 0.00001") # nonlocalECP
    list(APPEND DET_FeCO6_PYSCF_SHORT "samples" "144 0.0") # samples
else()
    list(APPEND DET_FeCO6_PYSCF_SHORT "kinetic" "165.56751508 0.000001") # kinetic energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "eeenergy" "569.46176293 0.000001") # e-e energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "totenergy" "-250.41330418 0.000001") # total energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "potential" "-415.98081925 0.000001") # potential energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "ionion" "422.43049172 0.000001") # i-i energy
    list(APPEND DET_FeCO6_PYSCF_SHORT "localecp" "-1438.65412200 0.000001") # localECP
    list(APPEND DET_FeCO6_PYSCF_SHORT "nonlocalecp" "30.78104810 0.000001") # nonlocalECP
    list(APPEND DET_FeCO6_PYSCF_SHORT "samples" "144 0.0") # samples
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
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "kinetic" "160.05153240 0.00001") # kinetic energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "eeenergy" "565.49285459 0.00001") # e-e energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "totenergy" "-248.30906741 0.00001") # total energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "potential" "-408.36060104 0.00001") # potential energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "ionion" "422.43049027 0.00001") # i-i energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "localecp" "-1440.43089300 0.00001") # localECP
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "nonlocalecp" "44.14694711 0.00001") # nonlocalECP
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "samples" "144 0.0") # samples
else()
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "kinetic" "161.33316595 0.000001") # kinetic energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "eeenergy" "569.38678149 0.000001") # e-e energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "totenergy" "-249.54304908 0.000001") # total energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "potential" "-410.87621503 0.000001") # potential energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "ionion" "422.43049172 0.000001") # i-i energy
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "localecp" "-1446.30226520 0.000001") # localECP
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "nonlocalecp" "43.60877699 0.000001") # nonlocalECP
    list(APPEND DET_FeCO6_PYSCF_SHORT_BATCHED "samples" "144 0.0") # samples
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
