
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
#
# [Fe(CO)6]2+ molecule gaussian tests, with BFD ECP
# Also check results for different number of mpi tasks and threads keeping total constant

# Ref for ref  VMC test
  LIST(APPEND FeCO6_GMS_SCALARS "kinetic" "163.8407 0.0020") # kinetic energy
  LIST(APPEND FeCO6_GMS_SCALARS "totenergy" "-250.43439 0.00032") # total energy
  LIST(APPEND FeCO6_GMS_SCALARS "eeenergy" "568.6009 0.0016") # e-e energy
  LIST(APPEND FeCO6_GMS_SCALARS "samples" "625868800 0.0") # samples
#  LIST(APPEND LIH_SCALARS "flux" "0.0 0.1") # sampling check, should be zero

  LIST(APPEND FeCO6_GMS_SHORT "kinetic" "163.8407 0.090") # kinetic energy
  LIST(APPEND FeCO6_GMS_SHORT "eeenergy" "568.6009 0.052") # e-e energy
  LIST(APPEND FeCO6_GMS_SHORT "totenergy" "-250.43439 0.0122") # total energy
  LIST(APPEND FeCO6_GMS_SHORT "samples" "48000 0.0") # samples
  QMC_RUN_AND_CHECK(short-FeCO6_gms-vmc_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_gms"
                    vmc_short_noj
                    vmc_short_noj.in.xml
                    1 16
                    TRUE
                    0 FeCO6_GMS_SHORT # VMC
                    )

  LIST(APPEND FeCO6_GMS_LONG "kinetic" "163.8407 0.045") # kinetic energy
  LIST(APPEND FeCO6_GMS_LONG "eeenergy" "568.6009 0.025") # e-e energy
  LIST(APPEND FeCO6_GMS_LONG "totenergy" "-250.43439 0.006") # total energy
  LIST(APPEND FeCO6_GMS_LONG "samples" "480000 0.0") # samples
  QMC_RUN_AND_CHECK(long-FeCO6_gms-vmc_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_gms"
                    vmc_long_noj
                    vmc_long_noj.in.xml
                    1 16
                    TRUE
                    0 FeCO6_GMS_LONG # VMC
                    )

 ELSE()
    MESSAGE("Skipping  FeCO6_b3lyp_gms tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping FeCO6_b3lyp_gms tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

