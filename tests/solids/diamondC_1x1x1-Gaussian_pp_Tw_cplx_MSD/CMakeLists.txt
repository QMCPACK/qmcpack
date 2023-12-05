if(QMC_COMPLEX)

  # VMC long run with no Jastrows
  # Compare directly to sCI Energy of: E=-11.65519310
  #Data generated using utils/make_ref_data.sh  diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG qmc_ref/qmc_long.s000.scalar.dat  qmc_ref/qmc_ref_long.s000.scalar.dat
  list(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "totenergy" "-11.65519310     0.006077601667")
  list(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "variance" "1.539759784521 0.027855789336")
  list(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "kinetic" "7.857471129807 0.019397520629")
  list(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "potential" "-19.513384734772 0.020252648227")
  list(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "eeenergy" "-2.835644692415 0.007764684716")
  list(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "localecp" "-5.360598387703 0.029116710254")
  list(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "nonlocalecp" "1.458525910786 0.010167449064")
  list(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "ionion" "-12.775667565000 0.100000000000")
  qmc_run_and_check(
    long-diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD"
    vmc_long
    vmc_long.in.xml
    1
    16
    TRUE
    0
    diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS #VMC
  )

  # VMC short run with no Jastrows
  # Compare directly to sCI Energy of: E=-11.65519310
  #Data generated using utils/make_ref_data.sh  diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT qmc_ref/qmc_short.s000.scalar.dat  qmc_ref/qmc_ref_long.s000.scalar.dat
  list(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "totenergy" "-11.65519310     0.011805237730")
  list(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "variance" "1.539759784521 0.211329760575")
  list(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "kinetic" "7.857471129807 0.053834125981")
  list(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "potential" "-19.513384734772 0.054910387369")
  list(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "eeenergy" "-2.835644692415 0.017280154413")
  list(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "localecp" "-5.360598387703 0.096994263787")
  list(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "nonlocalecp" "1.458525910786 0.029846773207")
  list(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "ionion" "-12.775667565000 0.100000000000")
  qmc_run_and_check(
    short-diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD"
    vmc_short
    vmc_short.in.xml
    1
    16
    TRUE
    0
    diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS #VMC
  )

else()
  message(VERBOSE "Skipping Complex Periodic LCAO  if Complex code not build (QMC_COMPLEX=0)")
endif()
