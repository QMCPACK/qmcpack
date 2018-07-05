
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
#
# LiH molecular dimer gaussian tests, all electron
# Also check results for different number of mpi tasks and threads keeping total constant
# Energy from Gamess: E=-7.9873236457
#
  LIST(APPEND LIH_SCALARS "kinetic" "7.991319   0.000066") # kinetic energy
  LIST(APPEND LIH_SCALARS "totenergy" "-7.9873131  0.0000088 ") # total energy
  LIST(APPEND LIH_SCALARS "eeenergy" "3.488823  0.000010") # e-e energy
  LIST(APPEND LIH_SCALARS "samples" "79134720000 0.0") # samples
#  LIST(APPEND LIH_SCALARS "flux" "0.0 0.1") # sampling check, should be zero

  LIST(APPEND LIH_SHORT "kinetic" "7.991319 0.006") # kinetic energy
  LIST(APPEND LIH_SHORT "totenergy" "-7.9873131 0.0006") # total energy
  LIST(APPEND LIH_SHORT "eeenergy" "3.488823 0.0007") # e-e energy
  LIST(APPEND LIH_SHORT "samples" "9600000 0.0") # samples
  QMC_RUN_AND_CHECK(short-LiH_dimer_ae_gms-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae_gms"
                    vmc_short_noj 
                    vmc_short_noj.in.xml
                    1 16
                    TRUE
                    0 LIH_SHORT # VMC
                    )

  LIST(APPEND LIH_LONG "kinetic" "7.991319 0.002") # kinetic energy
  LIST(APPEND LIH_LONG "totenergy" "-7.9873131 0.00025") # total energy
  LIST(APPEND LIH_LONG "eeenergy" "3.488823 0.0003") # e-e energy
  LIST(APPEND LIH_LONG "samples" "96000000 0.0") # samples
  QMC_RUN_AND_CHECK(long-LiH_dimer_ae_qp-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae_gms"
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

