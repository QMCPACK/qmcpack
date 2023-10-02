if(NOT QMC_COMPLEX)
  # H Positron test with SD and MD wavefunction
  #   Wavefunction reference energies:
  #   --------------------------------
  #   SCF (SD)  :  -0.6660424124598098
  #   CI  (MD)  :  -0.6742027301356266

  list(APPEND H_POSITRON_SD_SCALARS "kinetic" "0.66113370 0.003563")
  list(APPEND H_POSITRON_SD_SCALARS "totenergy" "-0.66610792 0.001381")
  list(APPEND H_POSITRON_SD_SCALARS "eeenergy" "-0.11912710 0.001801")
  list(APPEND H_POSITRON_SD_SCALARS "potential" "-1.32724162 0.003844")
  list(APPEND H_POSITRON_SD_SCALARS "samples" "160000 0.0")

  qmc_run_and_check(
    short-H_positron-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/H_positron"
    qmc_short-SD-H_Ground
    qmc_short-SD-H_Ground.in.xml
    1
    16
    TRUE
    0
    H_POSITRON_SD_SCALARS # VMC
  )

  qmc_run_and_check(
    short-H_positron-vmcbatch_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/H_positron"
    qmc_short_vmcbatch-SD-H_Ground
    qmc_short_vmcbatch-SD-H_Ground.in.xml
    1
    16
    TRUE
    0
    H_POSITRON_SD_SCALARS # VMC
  )

  list(APPEND H_POSITRON_MD_SCALARS "kinetic" "0.66325047 0.003595")
  list(APPEND H_POSITRON_MD_SCALARS "totenergy" "-0.67429644 0.001355")
  list(APPEND H_POSITRON_MD_SCALARS "eeenergy" "-0.14894710 0.001747")
  list(APPEND H_POSITRON_MD_SCALARS "potential" "-1.33754691 0.003872")
  list(APPEND H_POSITRON_MD_SCALARS "samples" "160000 0.0")

  qmc_run_and_check(
    short-H_positron_MSD-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/H_positron"
    qmc_short-MD-H_Ground
    qmc_short-MD-H_Ground.in.xml
    1
    16
    TRUE
    0
    H_POSITRON_MD_SCALARS # VMC
  )

  qmc_run_and_check(
    short-H_positron_MSD-vmcbatch_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/H_positron"
    qmc_short_vmcbatch-MD-H_Ground
    qmc_short_vmcbatch-MD-H_Ground.in.xml
    1
    16
    TRUE
    0
    H_POSITRON_MD_SCALARS # VMC
  )

  # Deterministic test

  if(QMC_MIXED_PRECISION)
    list(APPEND DET_H_POSITRON_SD_SCALARS "kinetic" "0.29681403 0.0001")
    list(APPEND DET_H_POSITRON_SD_SCALARS "totenergy" "-0.73754906 0.0001")
    list(APPEND DET_H_POSITRON_SD_SCALARS "eeenergy" "-0.35638536 0.00001")
    list(APPEND DET_H_POSITRON_SD_SCALARS "potential" "-1.03436309 0.00001")
  else()
    list(APPEND DET_H_POSITRON_SD_SCALARS "kinetic" "0.29671714 0.000001")
    list(APPEND DET_H_POSITRON_SD_SCALARS "totenergy" "-0.73764449 0.000001")
    list(APPEND DET_H_POSITRON_SD_SCALARS "eeenergy" "-0.35638597 0.000001")
    list(APPEND DET_H_POSITRON_SD_SCALARS "potential" "-1.03436163 0.000001")
  endif()

  qmc_run_and_check(
    deterministic-H_positron-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/H_positron"
    det_qmc-SD-H_Ground
    det_qmc-SD-H_Ground.in.xml
    1
    1
    TRUE
    0
    DET_H_POSITRON_SD_SCALARS # VMC
  )

  if(QMC_MIXED_PRECISION)
    list(APPEND DET_H_POSITRON_VMCBATCH_SD_SCALARS "kinetic" "0.29681399 0.0003")
    list(APPEND DET_H_POSITRON_VMCBATCH_SD_SCALARS "totenergy" "-0.73754911 0.0003")
    list(APPEND DET_H_POSITRON_VMCBATCH_SD_SCALARS "eeenergy" "-0.35638538 0.00001")
    list(APPEND DET_H_POSITRON_VMCBATCH_SD_SCALARS "potential" "-1.03436308 0.00001")
  else()
    list(APPEND DET_H_POSITRON_VMCBATCH_SD_SCALARS "kinetic" "0.29671714 0.000001")
    list(APPEND DET_H_POSITRON_VMCBATCH_SD_SCALARS "totenergy" "-0.73764449 0.000001")
    list(APPEND DET_H_POSITRON_VMCBATCH_SD_SCALARS "eeenergy" "-0.35638597 0.000001")
    list(APPEND DET_H_POSITRON_VMCBATCH_SD_SCALARS "potential" "-1.03436163 0.000001")
  endif()

  qmc_run_and_check(
    deterministic-H_positron-vmcbatch_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/H_positron"
    det_qmc_vmcbatch-SD-H_Ground
    det_qmc_vmcbatch-SD-H_Ground.in.xml
    1
    1
    TRUE
    0
    DET_H_POSITRON_VMCBATCH_SD_SCALARS # VMC
  )


  if(QMC_MIXED_PRECISION)	 
    list(APPEND DET_H_POSITRON_MD_SCALARS "kinetic" "0.99113870 0.0002")
    list(APPEND DET_H_POSITRON_MD_SCALARS "totenergy" "-0.61425842 0.0002")
    list(APPEND DET_H_POSITRON_MD_SCALARS "eeenergy" "-0.11788154 0.00001")
    list(APPEND DET_H_POSITRON_MD_SCALARS "potential" "-1.60539712 0.00001")
  else()
    list(APPEND DET_H_POSITRON_MD_SCALARS "kinetic" "0.99113311 0.000001")
    list(APPEND DET_H_POSITRON_MD_SCALARS "totenergy" "-0.61426433 0.000001")
    list(APPEND DET_H_POSITRON_MD_SCALARS "eeenergy" "-0.11788167 0.000001")
    list(APPEND DET_H_POSITRON_MD_SCALARS "potential" "-1.60539744 0.000001")
  endif()    

  qmc_run_and_check(
    deterministic-H_positron_MSD-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/H_positron"
    det_qmc-MD-H_Ground
    det_qmc-MD-H_Ground.in.xml
    1
    1
    TRUE
    0
    DET_H_POSITRON_MD_SCALARS # VMC
  )

  if(QMC_MIXED_PRECISION)
    list(APPEND DET_H_POSITRON_VMCBATCH_MD_SCALARS "kinetic" "0.99113875 0.0002")
    list(APPEND DET_H_POSITRON_VMCBATCH_MD_SCALARS "totenergy" "-0.61425823 0.0002")
    list(APPEND DET_H_POSITRON_VMCBATCH_MD_SCALARS "eeenergy" "-0.11788158 0.00001")
    list(APPEND DET_H_POSITRON_VMCBATCH_MD_SCALARS "potential" "-1.60539697 0.00001")
  else()
    list(APPEND DET_H_POSITRON_VMCBATCH_MD_SCALARS "kinetic" "0.99113311 0.000001")
    list(APPEND DET_H_POSITRON_VMCBATCH_MD_SCALARS "totenergy" "-0.61426433 0.000001")
    list(APPEND DET_H_POSITRON_VMCBATCH_MD_SCALARS "eeenergy" "-0.11788167 0.000001")
    list(APPEND DET_H_POSITRON_VMCBATCH_MD_SCALARS "potential" "-1.60539744 0.000001")
  endif()

  qmc_run_and_check(
    deterministic-H_positron_MSD-vmcbatch_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/H_positron"
    det_qmc_vmcbatch-MD-H_Ground
    det_qmc_vmcbatch-MD-H_Ground.in.xml
    1
    1
    TRUE
    0
    DET_H_POSITRON_VMCBATCH_MD_SCALARS # VMC
  )

else()
  message(VERBOSE "Skipping H_positron test because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
endif()
