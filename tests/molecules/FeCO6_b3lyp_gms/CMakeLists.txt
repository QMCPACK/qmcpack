if(NOT QMC_CUDA)
  if(NOT QMC_COMPLEX)
    #
    # [Fe(CO)6]2+ molecule gaussian tests, with BFD ECP using Gamess
    # Ref for ref  VMC test
    #  "kinetic"   "163.8407 0.0020") # kinetic energy
    #  "totenergy" "-250.43439 0.00032") # total energy
    #  "eeenergy"  "568.6009 0.0016") # e-e energy
    #  "samples"   "625868800 0.0") # samples

    list(APPEND FeCO6_GMS_SHORT "kinetic" "163.8407 0.23") # kinetic energy
    list(APPEND FeCO6_GMS_SHORT "eeenergy" "568.6009 0.18") # e-e energy
    list(APPEND FeCO6_GMS_SHORT "totenergy" "-250.43439 0.035") # total energy
    list(APPEND FeCO6_GMS_SHORT "samples" "48000 0.0") # samples
    qmc_run_and_check(
      short-FeCO6_gms-vmc_noj
      "${qmcpack_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_gms"
      vmc_short_noj
      vmc_short_noj.in.xml
      1
      16
      TRUE
      0
      FeCO6_GMS_SHORT # VMC
    )

    list(APPEND FeCO6_GMS_LONG "kinetic" "163.8407 0.088") # kinetic energy
    list(APPEND FeCO6_GMS_LONG "eeenergy" "568.6009 0.067") # e-e energy
    list(APPEND FeCO6_GMS_LONG "totenergy" "-250.43439 0.012") # total energy
    list(APPEND FeCO6_GMS_LONG "samples" "480000 0.0") # samples
    qmc_run_and_check(
      long-FeCO6_gms-vmc_noj
      "${qmcpack_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_gms"
      vmc_long_noj
      vmc_long_noj.in.xml
      1
      16
      TRUE
      0
      FeCO6_GMS_LONG # VMC
    )

if(QMC_MIXED_PRECISION)
    list(APPEND DET_FeCO6_GMS_SHORT "kinetic" "152.19971454 0.00001") # kinetic energy
    list(APPEND DET_FeCO6_GMS_SHORT "eeenergy" "563.06806164 0.00001") # e-e energy
    list(APPEND DET_FeCO6_GMS_SHORT "totenergy" "-251.28136179 0.00001") # total energy
    list(APPEND DET_FeCO6_GMS_SHORT "potential" "-403.48107633 0.00001") # potential energy
    list(APPEND DET_FeCO6_GMS_SHORT "ionion" "422.43052077 0.00001") # i-i energy
    list(APPEND DET_FeCO6_GMS_SHORT "localecp" "-1428.90698010 0.00001") # localECP
    list(APPEND DET_FeCO6_GMS_SHORT "nonlocalecp" "39.92732130 0.00001") # nonlocalECP
    list(APPEND DET_FeCO6_GMS_SHORT "samples" "144 0.0") # samples
else()
    list(APPEND DET_FeCO6_GMS_SHORT "kinetic" "162.27545500 0.000001") # kinetic energy
    list(APPEND DET_FeCO6_GMS_SHORT "eeenergy" "567.72642153 0.000001") # e-e energy
    list(APPEND DET_FeCO6_GMS_SHORT "totenergy" "-250.65585860 0.000001") # total energy
    list(APPEND DET_FeCO6_GMS_SHORT "potential" "-412.93131360 0.000001") # potential energy
    list(APPEND DET_FeCO6_GMS_SHORT "ionion" "422.43052231 0.000001") # i-i energy
    list(APPEND DET_FeCO6_GMS_SHORT "localecp" "-1439.95821000 0.000001") # localECP
    list(APPEND DET_FeCO6_GMS_SHORT "nonlocalecp" "36.86995255 0.000001") # nonlocalECP
    list(APPEND DET_FeCO6_GMS_SHORT "samples" "144 0.0") # samples
endif()

    qmc_run_and_check(
      deterministic-FeCO6_gms-vmc_noj
      "${qmcpack_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_gms"
      det_vmc_noj
      det_vmc_noj.in.xml
      1
      1
      TRUE
      0
      DET_FeCO6_GMS_SHORT # VMC
    )

if(QMC_MIXED_PRECISION)
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "kinetic" "152.71785045 0.00001") # kinetic energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "eeenergy" "566.17219082 0.00001") # e-e energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "totenergy" "-249.15602048 0.00001") # total energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "potential" "-401.87387108 0.00001") # potential energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "ionion" "422.43052077 0.00001") # i-i energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "localecp" "-1435.23597170 0.00001") # localECP
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "nonlocalecp" "44.75938902 0.00001") # nonlocalECP
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "samples" "144 0.0") # samples
else()
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "kinetic" "155.60906974 0.000001") # kinetic energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "eeenergy" "562.72044322 0.000001") # e-e energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "totenergy" "-251.54937078 0.000001") # total energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "potential" "-407.15844052 0.000001") # potential energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "ionion" "422.43052231 0.000001") # i-i energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "localecp" "-1430.23770790 0.000001") # localECP
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "nonlocalecp" "37.92830183 0.000001") # nonlocalECP
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "samples" "144 0.0") # samples
endif()

    qmc_run_and_check(
      deterministic-FeCO6_gms-vmcbatch_noj
      "${qmcpack_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_gms"
      det_vmcbatch_noj
      det_vmcbatch_noj.in.xml
      1
      1
      TRUE
      0
      DET_FeCO6_GMS_SHORT_BATCHED # VMC
    )

  else()
    message(VERBOSE
      "Skipping  FeCO6_b3lyp_gms tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
  endif()
else()
  message(VERBOSE
    "Skipping FeCO6_b3lyp_gms tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
endif()
