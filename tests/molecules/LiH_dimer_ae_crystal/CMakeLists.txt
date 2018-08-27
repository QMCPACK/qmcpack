
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
#
# Energy from Pyscf: E=-7.9873228116622
#
#   "kinetic" "7.9915  0.0016") # kinetic energy
#   "totenergy" "-7.98737  0.00028") # total energy
#   "eeenergy" "3.488935 0.000089") # e-e energy
#   "samples" "1152000000 0.0") # samples

  LIST(APPEND LIH_SHORT "kinetic" "7.9946 0.0033") # kinetic energy
  LIST(APPEND LIH_SHORT "totenergy" "-7.98651 0.00028") # total energy
  LIST(APPEND LIH_SHORT "eeenergy" "3.48940 0.00032") # e-e energy
  LIST(APPEND LIH_SHORT "samples" "115200000 0.0") # samples
  QMC_RUN_AND_CHECK(short-LiH_dimer_ae_crystal-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae_crystal"
                    vmc_short_noj
                    vmc_short_noj.in.xml
                    1 16
                    ${SOA_SUCCESS}
                    0 LIH_SHORT # VMC
                    )

  LIST(APPEND LIH_LONG "kinetic" "7.99230 0.00083") # kinetic energy
  LIST(APPEND LIH_LONG "totenergy" "-7.987230 0.000080") # total energy
  LIST(APPEND LIH_LONG "eeenergy" "3.48882 0.00010") # e-e energy
  LIST(APPEND LIH_LONG "samples" "967680000 0.0") # samples
  QMC_RUN_AND_CHECK(long-LiH_dimer_ae_crystal-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_ae_crystal"
                    vmc_long_noj
                    vmc_long_noj.in.xml
                    1 16
                    ${SOA_SUCCESS}
                    0 LIH_LONG # VMC
                    )


 ELSE()
    MESSAGE("Skipping LiH_dimer_ae tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping LiH_dimer_ae tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

