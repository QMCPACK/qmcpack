
IF(ENABLE_SOA)
 IF (NOT QMC_CUDA)
   IF (NOT QMC_COMPLEX)

# VMC long run with no Jastrows 
# Compare directly to Hartree Fock Energy of: E=-11.187187769284947
#Data generated using utils/make_ref_data.sh  diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG qmc-ref/qmc_long.s000.scalar.dat  qmc-ref/qmc_ref_long.s000.scalar.dat 
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "totenergy" "-11.187187769284947 0.006185568972")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "variance" "1.893849107587 0.031952937899")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "kinetic" "8.992601869357 0.031209738681")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "potential" "-20.180009436360 0.022888123934")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "eeenergy" "-2.603439510064 0.011081783142")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "localecp" "-5.969614858813 0.043988176158")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "nonlocalecp" "1.168712497974 0.012309208368")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "ionion" "-12.775667565000 0.000000000000")

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1-Gaussian_pp_Tw_cplx
                    "${CMAKE_SOURCE_DIR}/tests/solids/Al-1x1x1-Gaussian_pp_Tw_cplx"
                    vmc_long  
                    vmc_long.in.xml
                    1 16
                    TRUE
                    0 diamondC_1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS #VMC 
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
