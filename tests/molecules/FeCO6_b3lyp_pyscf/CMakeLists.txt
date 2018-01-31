
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
#
# LiH molecular dimer gaussian tests, all electron
# Also check results for different number of mpi tasks and threads keeping total constant
#
  LIST(APPEND LIH_SCALARS "kinetic" "165.48 0.36") # kinetic energy
  LIST(APPEND LIH_SCALARS "totenergy" "-250.8138 0.0015") # total energy
  LIST(APPEND LIH_SCALARS "eeenergy" "595.54 0.11") # e-e energy
  LIST(APPEND LIH_SCALARS "samples" "32000 0.0") # samples
#  LIST(APPEND LIH_SCALARS "flux" "0.0 0.1") # sampling check, should be zero

  QMC_RUN_AND_CHECK(long-FeCO6_pyscf-vmc_b3lyp_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_pyscf"
                    vmc_b3lyp_noj
                    vmc_b3lyp_noj.in.xml
                    16 1
                    TRUE
                    0 FeCO6_PYSCF_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(long-FeCO6_pyscf-vmc_b3lyp_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_pyscf"
                    vmc_b3lyp_noj
                    vmc_b3lyp_noj.in.xml
                    4 4
                    TRUE
                    0 FeCO6_PYSCF_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(long-FeCO6_pyscf-vmc_b3lyp_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/FeCO6_b3lyp_pyscf"
                    vmc_b3lyp_noj
                    vmc_b3lyp_noj.in.xml
                    1 16
                    TRUE
                    0 FeCO6_PYSCF_SCALARS # VMC
                    )
 ELSE()
    MESSAGE("Skipping  FeCO6_b3lyp_pyscf tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping FeCO6_b3lyp_pyscf tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

