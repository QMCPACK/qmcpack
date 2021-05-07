
IF (NOT QMC_CUDA)
  IF (NOT QMC_COMPLEX)
#
# [Fe(CO)6]2+ molecule gaussian tests, with BFD ECP using pyscf

# Ref for ref VMC test
#   "kinetic" "163.8420 0.0022") # kinetic energy
#   "totenergy" "-250.43499 0.00034") # total energy
#   "eeenergy" "568.6012 0.0017") # e-e energy
#   "samples" "655360000 0.0") # samples

  LIST(APPEND FeCO6_PYSCF_SHORT "kinetic" "163.8420 0.27") # kinetic energy
  LIST(APPEND FeCO6_PYSCF_SHORT "eeenergy" "568.6012 0.22") # e-e energy
  LIST(APPEND FeCO6_PYSCF_SHORT "totenergy" "-250.43499 0.035") # total energy
  LIST(APPEND FeCO6_PYSCF_SHORT "samples" "48000 0.0") # samples
  QMC_RUN_AND_CHECK(short-FeCO6_pyscf-vmc_noj
                    "${qmcpack_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_pyscf"
                    vmc_short_noj
                    vmc_short_noj.in.xml
                    1 16
                    TRUE
                    0 FeCO6_PYSCF_SHORT # VMC
                    )

  LIST(APPEND FeCO6_PYSCF_LONG "kinetic" "163.8420 0.076") # kinetic energy
  LIST(APPEND FeCO6_PYSCF_LONG "eeenergy" "568.6012 0.063") # e-e energy
  LIST(APPEND FeCO6_PYSCF_LONG "totenergy" "-250.43499 0.011") # total energy
  LIST(APPEND FeCO6_PYSCF_LONG "samples" "480000 0.0") # samples
  QMC_RUN_AND_CHECK(long-FeCO6_pyscf-vmc_noj
                    "${qmcpack_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_pyscf"
                    vmc_long_noj
                    vmc_long_noj.in.xml
                    1 16
                    TRUE
                    0 FeCO6_PYSCF_LONG # VMC
                    )

  ELSE()
    MESSAGE_VERBOSE("Skipping  FeCO6_b3lyp_gms tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
  ENDIF()
ELSE()
    MESSAGE_VERBOSE("Skipping FeCO6_b3lyp_gms tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

