if(NOT QMC_COMPLEX)

  # E_scf_QP=-150.95914112726956
  #  Reference from VMC runs
  list(APPEND C4_SD_AE_QP_SHORT "totenergy" "-150.95914112726956 0.074761412109")
  list(APPEND C4_SD_AE_QP_SHORT "variance" "34.741879231108 4.833758192799")
  list(APPEND C4_SD_AE_QP_SHORT "kinetic" "149.744668494694 3.0")
  list(APPEND C4_SD_AE_QP_SHORT "potential" "-300.705234559488 3.0")
  list(APPEND C4_SD_AE_QP_SHORT "eeenergy" "106.589716452578 0.093871118179")
  list(APPEND C4_SD_AE_QP_SHORT "ionion" "56.960598986000 0.001000000000")

  qmc_run_and_check(
    short-C4_sd-ae-qp-vmc
    "${qmcpack_SOURCE_DIR}/tests/molecules/C4_ae"
    qmc_short-SD-C4_AE_Mol_QP
    qmc_short-SD-C4_AE_Mol_QP.in.xml
    1
    16
    TRUE
    0
    C4_SD_AE_QP_SHORT # VMC
  )

  # Test case where number of orbitals is greater than the number of electrons
  qmc_run_and_check(
    short-C4_sd-ae-large-qp-vmc
    "${qmcpack_SOURCE_DIR}/tests/molecules/C4_ae"
    qmc_short-SD-C4_AE_Mol_large_QP
    qmc_short-SD-C4_AE_Mol_large_QP.in.xml
    1
    16
    TRUE
    0
    C4_SD_AE_QP_SHORT # VMC
  )

  list(APPEND C4_SD_AE_QP_LONG "totenergy" "-150.95914112726956 0.007790118606")
  list(APPEND C4_SD_AE_QP_LONG "variance" "34.741879231108 0.990741221699")
  list(APPEND C4_SD_AE_QP_LONG "kinetic" "149.744668494694 0.331776551096")
  list(APPEND C4_SD_AE_QP_LONG "potential" "-300.705234559488 0.330983245767")
  list(APPEND C4_SD_AE_QP_LONG "eeenergy" "106.589716452578 0.011179248515")
  list(APPEND C4_SD_AE_QP_LONG "ionion" "56.960598986000 0.001000000000")

  qmc_run_and_check(
    long-C4_sd-ae-qp-vmc
    "${qmcpack_SOURCE_DIR}/tests/molecules/C4_ae"
    qmc_long-SD-C4_AE_Mol_QP
    qmc_long-SD-C4_AE_Mol_QP.in.xml
    1
    16
    TRUE
    0
    C4_SD_AE_QP_LONG # VMC
  )

  #Nb_Det = 211
  #Ground:
  #E               =   -151.05936895837109

  #Excited:
  #E               =   -150.89066382584013

  list(APPEND C4_MD_AE_EXCITED_QP_SHORT "totenergy" "-150.89066382584013 0.062212972577")
  list(APPEND C4_MD_AE_EXCITED_QP_SHORT "variance" "34.255149340568 3.877696647688")
  list(APPEND C4_MD_AE_EXCITED_QP_SHORT "kinetic" "149.395366601073 2.289442218166")
  list(APPEND C4_MD_AE_EXCITED_QP_SHORT "potential" "-300.287675965416 2.302220898617")
  list(APPEND C4_MD_AE_EXCITED_QP_SHORT "eeenergy" "106.039594332889 0.113866937600")
  list(APPEND C4_MD_AE_EXCITED_QP_SHORT "ionion" "56.960598986000 0.001000000000")

  qmc_run_and_check(
    short-C4_md-ae-excited-qp-vmc
    "${qmcpack_SOURCE_DIR}/tests/molecules/C4_ae"
    qmc_short-MD-C4_AE_Mol_Excited_QP
    qmc_short-MD-C4_AE_Mol_Excited_QP.in.xml
    1
    16
    TRUE
    0
    C4_MD_AE_EXCITED_QP_SHORT # VMC
  )

  list(APPEND C4_MD_AE_GROUND_QP_SHORT "totenergy" "-151.05936895837109 0.067531501861")
  list(APPEND C4_MD_AE_GROUND_QP_SHORT "variance" "34.274239480694 1.882676106875")
  list(APPEND C4_MD_AE_GROUND_QP_SHORT "kinetic" "149.806496130140 2.275152959288")
  list(APPEND C4_MD_AE_GROUND_QP_SHORT "potential" "-300.866518836168 2.280577422331")
  list(APPEND C4_MD_AE_GROUND_QP_SHORT "eeenergy" "106.350153396650 0.098463197040")
  list(APPEND C4_MD_AE_GROUND_QP_SHORT "ionion" "56.960598986000 0.001000000000")

  qmc_run_and_check(
    short-C4_md-ae-ground-qp-vmc
    "${qmcpack_SOURCE_DIR}/tests/molecules/C4_ae"
    qmc_short-MD-C4_AE_Mol_Ground_QP
    qmc_short-MD-C4_AE_Mol_Ground_QP.in.xml
    1
    16
    TRUE
    0
    C4_MD_AE_GROUND_QP_SHORT # VMC
  )

  #--- Deterministic

  # Deterministic tests for 1 node & 1 thread
  ## FOR det_qmc_short-SD-C4_AE_Mol_QP.s000.scalar.dat
  #  LIST(APPEND DET_C4_SD_AE_QP "kinetic" "25.57801835 0.00001")
  #  LIST(APPEND DET_C4_SD_AE_QP "totenergy" "-147.14097018 0.00001")
  #  LIST(APPEND DET_C4_SD_AE_QP "potential" "-172.71898854 0.00001")
  #  LIST(APPEND DET_C4_SD_AE_QP "eeenergy" "101.53797646 0.00001")
  #  LIST(APPEND DET_C4_SD_AE_QP "ionion" "56.96059899 0.00001")
  #  LIST(APPEND DET_C4_SD_AE_QP "samples" "9 0.0")
  #  QMC_RUN_AND_CHECK(deterministic-C4_sd-ae-qp-vmc
  #                    "${qmcpack_SOURCE_DIR}/tests/molecules/C4_ae"
  #                    det_qmc_short-SD-C4_AE_Mol_QP
  #                    det_qmc_short-SD-C4_AE_Mol_QP.in.xml
  #                    1 1
  #                    TRUE
  #                    0 DET_C4_SD_AE_QP # VMC
  #                    )

  ## FOR det_qmc_short-MD-C4_AE_Mol_Ground_QP.s000.scalar.dat
  #  LIST(APPEND DET_C4_MD_AE_GROUND_QP "kinetic" "99.09812052 0.00001")
  #  LIST(APPEND DET_C4_MD_AE_GROUND_QP "totenergy" "-146.85148712 0.00001")
  #  LIST(APPEND DET_C4_MD_AE_GROUND_QP "potential" "-245.94960764 0.00001")
  #  LIST(APPEND DET_C4_MD_AE_GROUND_QP "eeenergy" "96.78479850  0.00001")
  #  LIST(APPEND DET_C4_MD_AE_GROUND_QP "ionion" "56.96059899 0.00001")
  #  LIST(APPEND DET_C4_MD_AE_GROUND_QP "samples" "9 0.0")
  #  QMC_RUN_AND_CHECK(deterministic-C4_md-ae-ground-qp-vmc
  #                    "${qmcpack_SOURCE_DIR}/tests/molecules/C4_ae"
  #                    det_qmc_short-MD-C4_AE_Mol_Ground_QP
  #                    det_qmc_short-MD-C4_AE_Mol_Ground_QP.in.xml
  #                    1 1
  #                    TRUE
  #                    0 DET_C4_MD_AE_GROUND_QP # VMC
  #                    )

  ## FOR det_qmc_short-MD-C4_AE_Mol_Excited_QP.s000.scalar.dat
  #  LIST(APPEND DET_C4_MD_AE_EXCITED_QP "kinetic" "25.92674520 0.00001")
  #  LIST(APPEND DET_C4_MD_AE_EXCITED_QP "totenergy" "-145.57864234 0.00001")
  #  LIST(APPEND DET_C4_MD_AE_EXCITED_QP "potential" "-171.50538754 0.00001")
  #  LIST(APPEND DET_C4_MD_AE_EXCITED_QP "eeenergy" "102.16100082 0.00001")
  #  LIST(APPEND DET_C4_MD_AE_EXCITED_QP "ionion" "56.96059899 0.00001")
  #  LIST(APPEND DET_C4_MD_AE_EXCITED_QP "samples" "9 0.0")
  #  QMC_RUN_AND_CHECK(deterministic-C4_md-ae-excited-qp-vmc
  #                    "${qmcpack_SOURCE_DIR}/tests/molecules/C4_ae"
  #                    det_qmc_short-MD-C4_AE_Mol_Excited_QP
  #                    det_qmc_short-MD-C4_AE_Mol_Excited_QP.in.xml
  #                    1 1
  #                    TRUE
  #                    0 DET_C4_MD_AE_EXCITED_QP # VMC
  #                     )

else()
  message(VERBOSE "Skipping C4_ae tests because Pure real WF are not supported by complex build (QMC_COMPLEX=1)")
endif()
