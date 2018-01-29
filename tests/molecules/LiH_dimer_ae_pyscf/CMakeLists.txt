
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
#
# LiH molecular dimer gaussian tests, all electron
# Also check results for different number of mpi tasks and threads keeping total constant
# Energy from Gamess: E=-7.9873236457148  
#
  LIST(APPEND LIH_SCALARS "kinetic" "8.027 0.043") # kinetic energy
  LIST(APPEND LIH_SCALARS "totenergy" "-7.9886 0.0035") # total energy
  LIST(APPEND LIH_SCALARS "eeenergy" "3.4947 0.0043") # e-e energy
  LIST(APPEND LIH_SCALARS "samples" "1600000 0.0") # samples
#  LIST(APPEND LIH_SCALARS "flux" "0.0 0.1") # sampling check, should be zero

  QMC_RUN_AND_CHECK(short-LiH_dimer_ae_pyscf-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae_pyscf"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml
                    16 1
                    ${SOA_SUCCESS}
                    0 LIH_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-LiH_dimer_ae_pyscf-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae_pyscf"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml
                    4 4
                    ${SOA_SUCCESS}
                    0 LIH_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-LiH_dimer_ae_pyscf-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae_pyscf"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml
                    1 16
                    ${SOA_SUCCESS}
                    0 LIH_SCALARS # VMC
                    )
 ELSE()
    MESSAGE("Skipping LiH_dimer_ae tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping LiH_dimer_ae tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

