
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
  IF (ENABLE_SOA)
#
# [Fe(CO)6]2+ molecule gaussian tests, with BFD ECP
# Also check results for different number of mpi tasks and threads keeping total constant

# Ref for long VMC test
  LIST(APPEND FeCO6_GMS_SCALARS "kinetic" "163.8420 0.0022") # kinetic energy
  LIST(APPEND FeCO6_GMS_SCALARS "totenergy" "-250.43499 0.00034") # total energy
  LIST(APPEND FeCO6_GMS_SCALARS "eeenergy" "568.6012 0.0017") # e-e energy
  LIST(APPEND FeCO6_GMS_SCALARS "samples" "16000 0.0") # samples
#  LIST(APPEND LIH_SCALARS "flux" "0.0 0.1") # sampling check, should be zero

  QMC_RUN_AND_CHECK(short-FeCO6_pyscf-vmc_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_pyscf"
                    qmc_short_noj
                    qmc_short_noj.in.xml
                    1 16
                    TRUE
                    0 FeCO6_GMS_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(long-FeCO6_pyscf-vmc_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_pyscf"
                    qmc_long_noj
                    qmc_long_noj.in.xml
                    1 16
                    TRUE
                    0 FeCO6_GMS_SCALARS # VMC
                    )

  ELSE()
    MESSAGE("Skipping  FeCO6_b3lyp_pyscf tests because Spherical representaion not supported by AoS (ENABLE_SOA=0)")
  ENDIF()
 ELSE()
    MESSAGE("Skipping  FeCO6_b3lyp_gms tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping FeCO6_b3lyp_gms tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

