
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
#
# LiH molecular dimer gaussian tests, all electron
# Also check results for different number of mpi tasks and threads keeping total constant
# Energy from Quantum Package: E=-7.987323645629950
#
  LIST(APPEND LIH_SCALARS "kinetic" "8.051 0.033") # kinetic energy
  LIST(APPEND LIH_SCALARS "totenergy" "-7.9833 0.0026") # total energy
  LIST(APPEND LIH_SCALARS "eeenergy" "3.4899 0.0022") # e-e energy
  LIST(APPEND LIH_SCALARS "samples" "1600000 0.0") # samples
#  LIST(APPEND LIH_SCALARS "flux" "0.0 0.1") # sampling check, should be zero

  QMC_RUN_AND_CHECK(short-LiH_dimer_ae_QP-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae_qp"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml 
                    16 1
                    TRUE
                    0 LIH_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-LiH_dimer_ae_qp-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae_qp"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml 
                    4 4
                    TRUE
                    0 LIH_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-LiH_dimer_ae_qp-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae_qp"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml 
                    1 16
                    TRUE
                    0 LIH_SCALARS # VMC
                    )
 ELSE()
    MESSAGE("Skipping LiH_dimer_ae tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping LiH_dimer_ae tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

