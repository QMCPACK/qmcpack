if(NOT QMC_CUDA)
  if(NOT QMC_COMPLEX)

    #C2_PP test targets the Multideterminant route of the code. The wavefunction data is  stored in a XML format AND HDF5 Format. Both files should be run
    # Energy from QP: E=-10.9015995458600 Ha
    #
    #   "kinetic" " 7.74018 0.00065") # kinetic energy
    #   "totenergy" " -11.060274 0.000055 ") # total energy
    #   "eeenergy" " 12.14369 0.00036") # e-e energy
    #   "samples" "67174400 0") # samples

    list(APPEND C2_PP_MSDJ_SHORT "kinetic" "7.74018 0.012") # kinetic energy
    list(APPEND C2_PP_MSDJ_SHORT "totenergy" "-11.060274 0.00099") # total energy
    list(APPEND C2_PP_MSDJ_SHORT "eeenergy" "12.14369 0.0052") # e-e energy
    list(APPEND C2_PP_MSDJ_SHORT "samples" "160000 0.0") # samples
    qmc_run_and_check(
      short-C2_pp_msdj-H5-vmc
      "${qmcpack_SOURCE_DIR}/tests/molecules/C2_pp"
      C2_TZ_PP_SHORT_H5_MSDJ
      vmc_short_pp_msdj-H5.xml
      1
      16
      TRUE
      1
      C2_PP_MSDJ_SHORT # VMC
    )

    qmc_run_and_check(
      short-C2_pp_msdj-H5-msd_no_precompute-vmc
      "${qmcpack_SOURCE_DIR}/tests/molecules/C2_pp"
      C2_TZ_PP_SHORT_H5_MSDJ
      vmc_short_pp_msd_no_precompute_jas-H5.xml
      1
      16
      TRUE
      1
      C2_PP_MSDJ_SHORT # VMC
    )

    qmc_run_and_check(
      short-C2_pp_msdj-XML-vmc
      "${qmcpack_SOURCE_DIR}/tests/molecules/C2_pp"
      C2_TZ_PP_SHORT_XML_MSDJ
      vmc_short_pp_msdj-XML.xml
      1
      16
      TRUE
      1
      C2_PP_MSDJ_SHORT # VMC
    )

    list(APPEND C2_PP_MSDJ_LONG "kinetic" "7.74018 0.0034") # kinetic energy
    list(APPEND C2_PP_MSDJ_LONG "totenergy" "-11.060274 0.00037") # total energy
    list(APPEND C2_PP_MSDJ_LONG "eeenergy" "12.14369 0.0020") # e-e energy
    list(APPEND C2_PP_MSDJ_LONG "samples" "1280000 0.0") # samples
    qmc_run_and_check(
      long-C2_pp_msdj-H5-vmc
      "${qmcpack_SOURCE_DIR}/tests/molecules/C2_pp"
      C2_TZ_PP_LONG_H5_MSDJ
      vmc_long_pp_msdj-H5.xml
      1
      16
      TRUE
      1
      C2_PP_MSDJ_LONG # VMC
    )

    qmc_run_and_check(
      long-C2_pp_msdj-XML-vmc
      "${qmcpack_SOURCE_DIR}/tests/molecules/C2_pp"
      C2_TZ_PP_LONG_XML_MSDJ
      vmc_long_pp_msdj-XML.xml
      1
      16
      TRUE
      1
      C2_PP_MSDJ_LONG # VMC
    )

  else()
    message(VERBOSE
      "Skipping C2_pp tests because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
  endif()
else()
  message(VERBOSE "Skipping C2_pp tests because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
endif()
