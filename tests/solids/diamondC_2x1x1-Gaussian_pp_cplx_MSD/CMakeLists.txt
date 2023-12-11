if(QMC_COMPLEX)

  # VMC long run with no Jastrows
  # Compare directly to sCI Energy of: E=-21.43955772
  #Data generated using utils/make_ref_data.sh  diamondC_2x1x1-Gaussian_pp_cplx_MSD_LONG qmc_ref/qmc_long.s000.scalar.dat  qmc_ref/qmc_ref_long.s000.scalar.dat
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_LONG_SCALARS "totenergy" "-21.43955772     0.016631586347")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_LONG_SCALARS "variance" "3.940110880232 0.143095003305")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_LONG_SCALARS "kinetic" "19.632464437269 0.067207073331")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_LONG_SCALARS "potential" "-41.072258549452 0.064597662591")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_LONG_SCALARS "eeenergy" "-4.745704876978 0.023930352973")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_LONG_SCALARS "localecp" "-12.898536216985 0.093725944996")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_LONG_SCALARS "nonlocalecp" "2.123309301169 0.027315413153")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_LONG_SCALARS "ionion" "-25.551326757000 0.100000000000")
  qmc_run_and_check(
    long-diamondC_2x1x1-Gaussian_pp_cplx_MSD
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_cplx_MSD"
    vmc_long
    vmc_long.in.xml
    1
    16
    TRUE
    0
    diamondC_2x1x1-Gaussian_pp_cplx_MSD_LONG_SCALARS #VMC
  )

  # VMC short run with no Jastrows
  # Compare directly to sCI Energy of: E=-21.43955772
  #Data generated using utils/make_ref_data.sh  diamondC_2x1x1-Gaussian_pp_cplx_MSD_SHORT qmc_ref/qmc_short.s000.scalar.dat  qmc_ref/qmc_ref_long.s000.scalar.dat
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_SHORT_SCALARS "totenergy" "-21.43955772     0.044296322552")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_SHORT_SCALARS "variance" "3.940110880232 0.290453940902")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_SHORT_SCALARS "kinetic" "19.632464437269 0.214155100497")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_SHORT_SCALARS "potential" "-41.072258549452 0.210512631857")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_SHORT_SCALARS "eeenergy" "-4.745704876978 0.052545086243")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_SHORT_SCALARS "localecp" "-12.898536216985 0.256538864675")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_SHORT_SCALARS "nonlocalecp" "2.123309301169 0.064408195985")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_SHORT_SCALARS "ionion" "-25.551326757000 0.220000000000")
  qmc_run_and_check(
    short-diamondC_2x1x1-Gaussian_pp_cplx_MSD
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_cplx_MSD"
    vmc_short
    vmc_short.in.xml
    1
    16
    TRUE
    0
    diamondC_2x1x1-Gaussian_pp_cplx_MSD_SHORT_SCALARS #VMC
  )

else()
  message(VERBOSE "Skipping Complex Periodic LCAO  if Complex code not build (QMC_COMPLEX=0)")
endif()
