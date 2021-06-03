
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
#
# LiH UHF molecule gaussian tests, all electron using pyscf
# Also check results for different number of mpi tasks and threads keeping total constant
# Energy from Pyscf: E=-7.98361527403875
#
#   "kinetic" "7.991344  0.000065") # kinetic energy
#   "totenergy" "-7.9873125  0.0000086 ") # total energy
#   "eeenergy" "3.488827 0.000010") # e-e energy
#   "samples" "80117760000 0.0") # samples

  LIST(APPEND LIH_UHF_SHORT "kinetic" "7.991344 0.015") # kinetic energy
  LIST(APPEND LIH_UHF_SHORT "totenergy" "-7.9873125 0.00062") # total energy
  LIST(APPEND LIH_UHF_SHORT "eeenergy" "3.488827 0.00089") # e-e energy
  LIST(APPEND LIH_UHF_SHORT "samples" "9600000 0.0") # samples
  QMC_RUN_AND_CHECK(short-LiH_ae_uhf_pyscf-vmc_hf_noj
                    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_pyscf_uhf"
                    vmc_short_noj
                    vmc_short_noj.in.xml
                    1 16
                    TRUE
                    0 LIH_UHF_SHORT # VMC
                    )

  LIST(APPEND LIH_UHF_LONG "kinetic" "7.991344 0.0022") # kinetic energy
  LIST(APPEND LIH_UHF_LONG "totenergy" "-7.9873125 0.000163") # total energy
  LIST(APPEND LIH_UHF_LONG "eeenergy" "3.488827 0.00033") # e-e energy
  LIST(APPEND LIH_UHF_LONG "samples" "96000000 0.0") # samples
  QMC_RUN_AND_CHECK(long-LiH_ae_uhf_pyscf-vmc_hf_noj
                    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_ae_pyscf_uhf"
                    vmc_long_noj
                    vmc_long_noj.in.xml
                    1 16
                    TRUE
                    0 LIH_UHF_LONG # VMC
                    )


 ELSE()
    MESSAGE_VERBOSE("Skipping LiH_ae_uhf tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
   MESSAGE_VERBOSE("Skipping LiH_ae_uhf tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

