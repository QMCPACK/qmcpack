
 IF (NOT QMC_CUDA)
   IF (QMC_COMPLEX)

# VMC long run with no Jastrows 
# Compare directly to Hartree Fock Energy of: E=-11.187187769284947
#Data generated using utils/make_ref_data.sh  diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG qmc_ref/qmc_long.s000.scalar.dat  qmc_ref/qmc_ref_long.s000.scalar.dat 
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "totenergy" "-11.187187769284947 0.006185568972")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "variance" "1.893849107587 0.031952937899")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "kinetic" "8.992601869357 0.031209738681")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "potential" "-20.180009436360 0.022888123934")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "eeenergy" "-2.603439510064 0.011081783142")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "localecp" "-5.969614858813 0.043988176158")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "nonlocalecp" "1.168712497974 0.012309208368")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "ionion" "-12.775667565000 0.020000000000")
  QMC_RUN_AND_CHECK(long-diamondC_1x1x1-Gaussian_pp_Tw_cplx
                    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_Tw_cplx"
                    vmc_long  
                    vmc_long.in.xml
                    1 16
                    TRUE
                    0 diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS #VMC 
                    )

# VMC short run with no Jastrows 
# Compare directly to Hartree Fock Energy of: E=-11.187187769284947
#Data generated using utils/make_ref_data.sh  diamondC_1x1x1-Gaussian_pp_Tw_cplx_SHORT qmc_ref/qmc_short.s000.scalar.dat  qmc_ref/qmc_ref_long.s000.scalar.dat 
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_SHORT_SCALARS "totenergy" "-11.187187769284947 0.015799832919")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_SHORT_SCALARS "variance" "1.893849107587 0.087708561871")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_SHORT_SCALARS "kinetic" "8.992601869357 0.048503358207")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_SHORT_SCALARS "potential" "-20.180009436360 0.051663837931")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_SHORT_SCALARS "eeenergy" "-2.603439510064 0.018396250955")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_SHORT_SCALARS "localecp" "-5.969614858813 0.068868300966")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_SHORT_SCALARS "nonlocalecp" "1.168712497974 0.025981779181")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_SHORT_SCALARS "ionion" "-12.775667565000 0.020000000000")
  QMC_RUN_AND_CHECK(short-diamondC_1x1x1-Gaussian_pp_Tw_cplx
                    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_Tw_cplx"
                    vmc_short 
                    vmc_short.in.xml
                    1 16
                    TRUE
                    0 diamondC_1x1x1-Gaussian_pp_Tw_cplx_SHORT_SCALARS #VMC 
                    )




    ELSE()
      MESSAGE_VERBOSE("Skipping Complex Periodic LCAO  if Complex code not build (QMC_COMPLEX=0)")
    ENDIF()
  ELSE()
    MESSAGE_VERBOSE("Skipping Periodic LCAO as not supported by CUDA build (QMC_CUDA=1)")
  ENDIF()
