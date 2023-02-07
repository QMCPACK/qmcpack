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
    list(APPEND DET_FeCO6_GMS_SHORT "kinetic" "139.26843221 0.003") # kinetic energy
    list(APPEND DET_FeCO6_GMS_SHORT "eeenergy" "568.76484836 0.001") # e-e energy
    list(APPEND DET_FeCO6_GMS_SHORT "totenergy" "-250.11452519 0.001") # total energy
    list(APPEND DET_FeCO6_GMS_SHORT "potential" "-389.38295740 0.001") # potential energy
    list(APPEND DET_FeCO6_GMS_SHORT "ionion" "422.43052077 0.0001") # i-i energy
    list(APPEND DET_FeCO6_GMS_SHORT "localecp" "-1405.32291650 0.001") # localECP
    list(APPEND DET_FeCO6_GMS_SHORT "nonlocalecp" "24.74459001 0.002") # nonlocalECP
    list(APPEND DET_FeCO6_GMS_SHORT "samples" "9 0.0") # samples
  else()
    list(APPEND DET_FeCO6_GMS_SHORT "kinetic" "163.17377427 0.000001") # kinetic energy
    list(APPEND DET_FeCO6_GMS_SHORT "eeenergy" "561.38303067 0.000001") # e-e energy
    list(APPEND DET_FeCO6_GMS_SHORT "totenergy" "-252.63553526 0.000001") # total energy
    list(APPEND DET_FeCO6_GMS_SHORT "potential" "-415.80930953 0.000001") # potential energy
    list(APPEND DET_FeCO6_GMS_SHORT "ionion" "422.43052231 0.000001") # i-i energy
    list(APPEND DET_FeCO6_GMS_SHORT "localecp" "-1432.09655360 0.000001") # localECP
    list(APPEND DET_FeCO6_GMS_SHORT "nonlocalecp" "32.47369109 0.000001") # nonlocalECP
    list(APPEND DET_FeCO6_GMS_SHORT "samples" "9 0.0") # samples
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
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "kinetic" "134.10184342 0.003") # kinetic energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "eeenergy" "552.72060729 0.002") # e-e energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "totenergy" "-251.29185486 0.001") # total energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "potential" "-385.39369691 0.002") # potential energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "ionion" "422.43052077 0.0001") # i-i energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "localecp" "-1391.12666390 0.003") # localECP
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "nonlocalecp" "30.58183893 0.002") # nonlocalECP
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "samples" "9 0.0") # samples
  else()
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "kinetic" "163.17377427 0.000001") # kinetic energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "eeenergy" "561.38303067 0.000001") # e-e energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "totenergy" "-252.63553526 0.000001") # total energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "potential" "-415.80930953 0.000001") # potential energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "ionion" "422.43052231 0.000001") # i-i energy
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "localecp" "-1432.09655360 0.000001") # localECP
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "nonlocalecp" "32.47369109 0.000001") # nonlocalECP
    list(APPEND DET_FeCO6_GMS_SHORT_BATCHED "samples" "9 0.0") # samples
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
  message(
    VERBOSE
    "Skipping  FeCO6_b3lyp_gms tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
endif()
