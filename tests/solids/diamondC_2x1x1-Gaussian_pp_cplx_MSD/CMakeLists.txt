
IF(ENABLE_SOA)
 IF (NOT QMC_CUDA)
   IF (QMC_COMPLEX)

# VMC long run with no Jastrows 
# Compare directly to sCI Energy of: E=-21.43955772
#Data generated using utils/make_ref_data.sh  diamondC_2x1x1-Gaussian_pp_cplx_MSD_LONG qmc_ref/qmc_long.s000.scalar.dat  qmc_ref/qmc_ref_long.s000.scalar.dat 
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "totenergy" "-11.76861606 0.007202242728")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "variance" "1.000602519743 0.045687701277")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "kinetic" "8.081902038337 0.046347290312")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "potential" "-19.850035624680 0.047014871698")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "eeenergy" "-3.014169096687 0.014078543711")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "localecp" "-5.496998282285 0.067595358823")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "nonlocalecp" "1.436799319733 0.024004326954")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_LONG_SCALARS "ionion" "-12.775667565000 0.100000000000")
  QMC_RUN_AND_CHECK(long-diamondC_2x1x1-Gaussian_pp_cplx_MSD
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_cplx_MSD"
                    vmc_long  
                    vmc_long.in.xml
                    1 16
                    TRUE
                    0 diamondC_2x1x1-Gaussian_pp_cplx_MSD_LONG_SCALARS #VMC 
                    )

# VMC short run with no Jastrows 
# Compare directly to sCI Energy of: E=-21.43955772
#Data generated using utils/make_ref_data.sh  diamondC_2x1x1-Gaussian_pp_cplx_MSD_SHORT qmc_ref/qmc_short.s000.scalar.dat  qmc_ref/qmc_ref_long.s000.scalar.dat 
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "totenergy" "-11.76861606 0.013373782570")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "variance" "1.000602519743 0.063924769247")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "kinetic" "8.081902038337 0.102927629801")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "potential" "-19.850035624680 0.105507379338")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "eeenergy" "-3.014169096687 0.026116076315")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "localecp" "-5.496998282285 0.140482720505")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "nonlocalecp" "1.436799319733 0.046041270355")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_MSD_SHORT_SCALARS "ionion" "-12.775667565000 0.200000000000")
  QMC_RUN_AND_CHECK(short-diamondC_2x1x1-Gaussian_pp_cplx_MSD
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_cplx_MSD"
                    vmc_short 
                    vmc_short.in.xml
                    1 16
                    TRUE
                    0 diamondC_2x1x1-Gaussian_pp_cplx_MSD_SHORT_SCALARS #VMC 
                    )




    ELSE()
      MESSAGE_VERBOSE("Skipping Complex Periodic LCAO  if Complex code not build (QMC_COMPLEX=0)")
    ENDIF()
  ELSE()
    MESSAGE_VERBOSE("Skipping Periodic LCAO as not supported by CUDA build (QMC_CUDA=1)")
  ENDIF()
ELSE()
  MESSAGE_VERBOSE("Skipping Periodic LCAO as not supported by AoS build (ENABLE_SOA=0)")
ENDIF()
