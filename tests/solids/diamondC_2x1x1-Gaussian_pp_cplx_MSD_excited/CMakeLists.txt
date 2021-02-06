
 IF (NOT QMC_CUDA)
   IF (QMC_COMPLEX)







      

#Data generated using utils/make_ref_data.sh  diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG qmc_ref/vmc_long_excited.s000.scalar.dat qmc_ref/qmc_ref_excited_long.s000.scalar.dat
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS "totenergy" "-21.043422517714 0.008715737370")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS "variance" "1.449654318091 0.023049744979")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS "kinetic" "21.005794394605 0.052657189151")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS "potential" "-42.049216912354 0.053155315732")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS "eeenergy" "-5.428377875848 0.014377274802")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS "localecp" "-13.053874507313 0.067010768286")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS "nonlocalecp" "1.984362532886 0.019826769332")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS "ionion" "-25.551327062000 0.001000000000")
  QMC_RUN_AND_CHECK(long-diamondC_2x1x1-Gaussian_pp_cplx_excited_MSD
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_cplx_MSD_excited"                    
                    vmc_long_excited
                    vmc_long_excited.in.xml 
                    1 16
                    TRUE
                    0 diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_LONG_SCALARS #VMC                              
                    )
      


    ELSE()
      MESSAGE_VERBOSE("Skipping Complex Periodic LCAO  if Complex code not build (QMC_COMPLEX=0)")
    ENDIF()
  ELSE()
    MESSAGE_VERBOSE("Skipping Periodic LCAO as not supported by CUDA build (QMC_CUDA=1)")
  ENDIF()
