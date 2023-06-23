if(QMC_COMPLEX)

  # VMC long run with no Jastrows
  # Compare directly to Hartree Fock Energy of: E=-21.30438154915548
  #Data generated using utils/make_ref_data.sh  diamondC_2x1x1-Gaussian_pp_cplx_LONG qmc_ref/qmc_long.s000.scalar.dat  qmc_ref/qmc_ref_long.s000.scalar.dat
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_LONG_SCALARS "totenergy" "-21.30438154915548 0.021337507878")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_LONG_SCALARS "variance" "4.388589388914 0.103220760086")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_LONG_SCALARS "kinetic" "19.545637226792 0.077931426068")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_LONG_SCALARS "potential" "-40.850779026183 0.081804233229")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_LONG_SCALARS "eeenergy" "-4.537323794961 0.028697376864")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_LONG_SCALARS "localecp" "-12.875979090005 0.104598406850")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_LONG_SCALARS "nonlocalecp" "2.113850615398 0.028971593529")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_LONG_SCALARS "ionion" "-25.551326757000 0.020000000000")
  qmc_run_and_check(
    long-diamondC_2x1x1-Gaussian_pp_Tw_cplx
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_cplx"
    vmc_long
    vmc_long.in.xml
    1
    16
    TRUE
    0
    diamondC_2x1x1-Gaussian_pp_cplx_LONG_SCALARS #VMC
  )

  # VMC short run with no Jastrows
  # Compare directly to Hartree Fock Energy of: E=-21.30438154915548
  #Data generated using utils/make_ref_data.sh  diamondC_2x1x1-Gaussian_pp_cplx_SHORT qmc_ref/qmc_short.s000.scalar.dat  qmc_ref/qmc_ref_long.s000.scalar.dat
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_SHORT_SCALARS "totenergy" "-21.30438154915548 0.045341262663")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_SHORT_SCALARS "variance" "4.388589388914 0.25")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_SHORT_SCALARS "kinetic" "19.545637226792 0.133872779266")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_SHORT_SCALARS "potential" "-40.850779026183 0.130602718601")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_SHORT_SCALARS "eeenergy" "-4.537323794961 0.058408305166")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_SHORT_SCALARS "localecp" "-12.875979090005 0.173005314408")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_SHORT_SCALARS "nonlocalecp" "2.113850615398 0.059863222645")
  list(APPEND diamondC_2x1x1-Gaussian_pp_cplx_SHORT_SCALARS "ionion" "-25.551326757000 0.020000000000")
  qmc_run_and_check(
    short-diamondC_2x1x1-Gaussian_pp_Tw_cplx
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_cplx"
    vmc_short
    vmc_short.in.xml
    1
    16
    TRUE
    0
    diamondC_2x1x1-Gaussian_pp_cplx_SHORT_SCALARS #VMC
  )

else()
  message(VERBOSE "Skipping Complex Periodic LCAO  if Complex code not build (QMC_COMPLEX=0)")
endif()
