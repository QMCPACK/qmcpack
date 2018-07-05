
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
#
# LiH molecular dimer gaussian tests, all electron
# Also check results for different number of mpi tasks and threads keeping total constant
# Energy from Quantum Package: E=-7.987323645629950
#
  LIST(APPEND LIH_SCALARS "kinetic" "7.991334  0.000065") # kinetic energy
  LIST(APPEND LIH_SCALARS "totenergy" "-7.9873140 0.0000086 ") # total energy
  LIST(APPEND LIH_SCALARS "eeenergy" "3.488827 0.000011") # e-e energy
  LIST(APPEND LIH_SCALARS "samples" "79626240000 0.0") # samples
#  LIST(APPEND LIH_SCALARS "flux" "0.0 0.1") # sampling check, should be zero

  LIST(APPEND LIH_SHORT "kinetic" "7.991334 0.006") # kinetic energy
  LIST(APPEND LIH_SHORT "totenergy" "-7.9873140 0.0008") # total energy
  LIST(APPEND LIH_SHORT "eeenergy" "3.488827 0.0008") # e-e energy
  LIST(APPEND LIH_SHORT "samples" "9600000 0.0") # samples
  QMC_RUN_AND_CHECK(short-LiH_dimer_ae_qp-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae_qp"
                    vmc_short_noj 
                    vmc_short_noj.in.xml
                    1 16
                    TRUE
                    0 LIH_SHORT # VMC
                    )

  LIST(APPEND LIH_LONG "kinetic" "7.991334 0.002") # kinetic energy
  LIST(APPEND LIH_LONG "totenergy" "-7.9873140 0.0003") # total energy
  LIST(APPEND LIH_LONG "eeenergy" "3.488827 0.0003") # e-e energy
  LIST(APPEND LIH_LONG "samples" "96000000 0.0") # samples
  QMC_RUN_AND_CHECK(long-LiH_dimer_ae_qp-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae_qp"
                    vmc_long_noj 
                    vmc_long_noj.in.xml 
                    1 16
                    TRUE
                    0 LIH_LONG # VMC
                    )


 ELSE()
    MESSAGE("Skipping LiH_dimer_ae tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping LiH_dimer_ae tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

