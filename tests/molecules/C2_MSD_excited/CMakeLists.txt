IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)

#C2_MSD_Excited test targets the Multideterminant route of the code with excited states. The wavefunction data is  stored in a HDF5 Format where the multideterminants expansions are the same for the ground state and excited state (ext_level=1 and ext_level=2) but the coefficients are stored in two separate arrays. Both Ground State and Excited state should be ran. Since it is an all electron VMC calculation with no Jastrow, we should reproduce QP2 results.  
# Energy from QP2 Ground : E=-75.64507811 Ha
# Energy from QP2 Excited: E=-75.35872670 Ha
#
#   VMC REF - Ground State
#   "totenergy" " -75.64590 0.00036 ") # total energy
#   "kinetic" " 75.676 0.013") # kinetic energy
#   "eeenergy" " 40.37148 0.00070") # e-e energy
#   "samples" "450150400  0") # samples
#
#
#
#   VMC REF - Excited State
#   "totenergy" "-75.35987 0.00028  ") # total energy
#   "kinetic" " 75.632 0.010") # kinetic energy
#   "eeenergy" " 40.67007 0.00064") # e-e energy
#   "samples" "573440000  0") # samples

  LIST(APPEND C2_DZ_GROUND_SHORT "kinetic" "75.676 1.02") # kinetic energy
  LIST(APPEND C2_DZ_GROUND_SHORT "totenergy" "-75.64590 0.024") # total energy
  LIST(APPEND C2_DZ_GROUND_SHORT "eeenergy" "40.37148 0.041") # e-e energy
  LIST(APPEND C2_DZ_GROUND_SHORT "samples" "160000 0.0") # samples
  QMC_RUN_AND_CHECK(short-C2_msd-ground-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C2_MSD_excited"
                    C2_DZ_SHORT_GROUND_MSD 
                    vmc_short_msd-Ground.xml 
                    1 16
                    TRUE
                    1 C2_DZ_GROUND_SHORT # VMC
                    )

  LIST(APPEND C2_DZ_GROUND_LONG "kinetic" "75.676 0.49") # kinetic energy
  LIST(APPEND C2_DZ_GROUND_LONG "totenergy" "-75.64590 0.0079") # total energy
  LIST(APPEND C2_DZ_GROUND_LONG "eeenergy" "40.37148 0.015") # e-e energy
  LIST(APPEND C2_DZ_GROUND_LONG "samples" "160000 0.0") # samples
  QMC_RUN_AND_CHECK(long-C2_msd-ground-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C2_MSD_excited"
                    C2_DZ_LONG_GROUND_MSD 
                    vmc_long_msd-Ground.xml 
                    1 16
                    TRUE
                    1 C2_DZ_GROUND_LONG # VMC
                    )

  LIST(APPEND C2_DZ_EXCITED_SHORT "kinetic" "75.632 0.85") # kinetic energy
  LIST(APPEND C2_DZ_EXCITED_SHORT "totenergy" "-75.35987 0.030") # total energy
  LIST(APPEND C2_DZ_EXCITED_SHORT "eeenergy" "40.67007 0.038") # e-e energy
  LIST(APPEND C2_DZ_EXCITED_SHORT "samples" "160000 0.0") # samples
  QMC_RUN_AND_CHECK(short-C2_msd-excited-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C2_MSD_excited"
                    C2_DZ_SHORT_EXT2_MSD 
                    vmc_short_msd-Excited2.xml 
                    1 16
                    TRUE
                    1 C2_DZ_EXCITED_SHORT # VMC
                    )

  LIST(APPEND C2_DZ_EXCITED_LONG "kinetic" "75.632 0.34") # kinetic energy
  LIST(APPEND C2_DZ_EXCITED_LONG "totenergy" "-75.35987 0.0098") # total energy
  LIST(APPEND C2_DZ_EXCITED_LONG "eeenergy" "40.67007 0.013") # e-e energy
  LIST(APPEND C2_DZ_EXCITED_LONG "samples" "160000 0.0") # samples
  QMC_RUN_AND_CHECK(long-C2_msd-excited-vmc
                    "${CMAKE_SOURCE_DIR}/tests/molecules/C2_MSD_excited"
                    C2_DZ_LONG_EXT2_MSD 
                    vmc_Long_msd-Excited2.xml 
                    1 16
                    TRUE
                    1 C2_DZ_EXCITED_LONG # VMC
                    )

 ELSE()
    MESSAGE_VERBOSE("Skipping C2_MSD_Excited tests because Pure real WF are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE_VERBOSE("Skipping C2_MSD_Excited tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()
