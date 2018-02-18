
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
#
# LiH molecular dimer gaussian tests, all electron
# Also check results for different number of mpi tasks and threads keeping total constant
#
  LIST(APPEND LIH_SCALARS "kinetic" "163.30 4.24") # kinetic energy
  LIST(APPEND LIH_SCALARS "totenergy" "-250.415 4.24") # total energy
  LIST(APPEND LIH_SCALARS "eeenergy" "568.46 0.13") # e-e energy
  LIST(APPEND LIH_SCALARS "samples" "32000 0.0") # samples
#  LIST(APPEND LIH_SCALARS "flux" "0.0 0.1") # sampling check, should be zero

  QMC_RUN_AND_CHECK(long-FeCO6_gms-vmc_b3lyp_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_gms"
                    vmc_b3lyp_noj
                    vmc_b3lyp_noj.in.xml
                    16 1
                    TRUE
                    0 FeCO6_GMS_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(long-FeCO6_gms-vmc_b3lyp_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_gms"
                    vmc_b3lyp_noj
                    vmc_b3lyp_noj.in.xml
                    4 4
                    TRUE
                    0 FeCO6_GMS_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(long-FeCO6_gms-vmc_b3lyp_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_gms"
                    vmc_b3lyp_noj
                    vmc_b3lyp_noj.in.xml
                    1 16
                    TRUE
                    0 FeCO6_GMS_SCALARS # VMC
                    )
 ELSE()
    MESSAGE("Skipping  FeCO6_b3lyp_gms tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping FeCO6_b3lyp_gms tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

