if(NOT QMC_CUDA)

  list(APPEND HE2_SP_SCALARS "totenergy" "-2.32676861 0.0863")
  list(APPEND HE2_SP_SCALARS "kinetic" "16.23198428 0.1117")
  list(APPEND HE2_SP_SCALARS "potential" "-18.55875289 0.1136")
  list(APPEND HE2_SP_SCALARS "eeenergy" "4.59673345 0.0174")
  list(APPEND HE2_SP_SCALARS "ionion" "0.86420131 0.0001")
  list(APPEND HE2_SP_SCALARS "samples" "16000 0.0")

  qmc_run_and_check(
    short-he2_1x1x1_ae-vmc_gaussian-sp-sph_noj
    "${qmcpack_SOURCE_DIR}/tests/solids/he2_1x1x1_ae"
    he2_sp_sph
    he2_sp_sph.qmc.in-wfnoj.xml
    1
    16
    TRUE
    0
    HE2_SP_SCALARS # VMC
  )

  list(APPEND HE2_SD_SCALARS "totenergy" "-2.32942007 0.0835")
  list(APPEND HE2_SD_SCALARS "kinetic" "16.21547981 0.1067")
  list(APPEND HE2_SD_SCALARS "potential" "-18.54489989 0.1189")
  list(APPEND HE2_SD_SCALARS "eeenergy" "4.59377532 0.0184")
  list(APPEND HE2_SD_SCALARS "ionion" "0.86420131 0.0001")
  list(APPEND HE2_SD_SCALARS "samples" "16000 0.0")

  qmc_run_and_check(
   short-he2_1x1x1_ae-vmc_gaussian-sd-sph_noj
   "${qmcpack_SOURCE_DIR}/tests/solids/he2_1x1x1_ae"
   he2_sd_sph
   he2_sd_sph.qmc.in-wfnoj.xml
   1
   16
   TRUE
   0
   HE2_SD_SCALARS # VMC
  )

  list(APPEND HE2_SF_SCALARS "totenergy" "-2.35068897 0.1249")
  list(APPEND HE2_SF_SCALARS "kinetic" "16.20083403 0.1568")
  list(APPEND HE2_SF_SCALARS "potential" "-18.55152300 0.1705")
  list(APPEND HE2_SF_SCALARS "eeenergy" "4.59192604 0.0261")
  list(APPEND HE2_SF_SCALARS "ionion" "0.86420131 0.0001")
  list(APPEND HE2_SF_SCALARS "samples" "16000 0.0")

  qmc_run_and_check(
   short-he2_1x1x1_ae-vmc_gaussian-sf-sph_noj
   "${qmcpack_SOURCE_DIR}/tests/solids/he2_1x1x1_ae"
   he2_sf_sph
   he2_sf_sph.qmc.in-wfnoj.xml
   1
   16
   TRUE
   0
   HE2_SF_SCALARS # VMC
  )

  list(APPEND HE2_SG_SCALARS "totenergy" "-2.30636188 0.1233")
  list(APPEND HE2_SG_SCALARS "kinetic" "16.21415978 0.1493")
  list(APPEND HE2_SG_SCALARS "potential" "-18.52052166 0.1655")
  list(APPEND HE2_SG_SCALARS "eeenergy" "4.59711595 0.0251")
  list(APPEND HE2_SG_SCALARS "ionion" "0.86420131 0.0001")
  list(APPEND HE2_SG_SCALARS "samples" "16000 0.0")

  qmc_run_and_check(
   short-he2_1x1x1_ae-vmc_gaussian-sg-sph_noj
   "${qmcpack_SOURCE_DIR}/tests/solids/he2_1x1x1_ae"
   he2_sg_sph
   he2_sg_sph.qmc.in-wfnoj.xml
   1
   16
   TRUE
   0
   HE2_SG_SCALARS # VMC
  )

  list(APPEND HE2_SH_SCALARS "totenergy" "-2.31529851 0.0837")
  list(APPEND HE2_SH_SCALARS "kinetic" "16.21254352 0.1095")
  list(APPEND HE2_SH_SCALARS "potential" "-18.52784203 0.1163")
  list(APPEND HE2_SH_SCALARS "eeenergy" "4.59749106 0.0174")
  list(APPEND HE2_SH_SCALARS "ionion" " 0.86420131 0.0001")
  list(APPEND HE2_SH_SCALARS "samples" "16000 0.0")

  qmc_run_and_check(
   short-he2_1x1x1_ae-vmc_gaussian-sh-sph_noj
   "${qmcpack_SOURCE_DIR}/tests/solids/he2_1x1x1_ae"
   he2_sh_sph
   he2_sh_sph.qmc.in-wfnoj.xml
   1
   16
   TRUE
   0
   HE2_SH_SCALARS # VMC
  )

else()
  message(VERBOSE "Skipping Periodic LCAO as not supported by CUDA build (QMC_CUDA=1)")
endif()
