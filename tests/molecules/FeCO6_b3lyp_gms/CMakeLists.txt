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

  else()
    message_verbose(
      "Skipping  FeCO6_b3lyp_gms tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
  endif()
else()
  message_verbose(
    "Skipping FeCO6_b3lyp_gms tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
endif()
