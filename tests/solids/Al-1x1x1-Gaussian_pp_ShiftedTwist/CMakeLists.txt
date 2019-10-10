
IF(ENABLE_SOA)
 IF (NOT QMC_CUDA)
   IF (NOT QMC_COMPLEX)

# VMC long run with no Jastrows 
# Compare directly to Hartree Fock Energy of: E=-8.221102284516801
#Data generated using utils/make_ref_data.sh  Al-1x1x1-Gaussian_pp_Tw_cplx_LONG qmc_ref/qmc_ref_short.s000.scalar.dat  qmc_ref/qmc_ref_long.s000.scalar.dat 

LIST(APPEND Al-1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "totenergy" "-8.221102284516801 0.010033175665")
LIST(APPEND Al-1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "variance" "0.925547990401 0.017872039850")
LIST(APPEND Al-1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "kinetic" "3.405591892337 0.007761089867")
LIST(APPEND Al-1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "potential" "-11.626261269528 0.012650810074")
LIST(APPEND Al-1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "eeenergy" "-3.102527469122 0.009855539320")
LIST(APPEND Al-1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "localecp" "0.656663039125 0.020622464088")
LIST(APPEND Al-1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "nonlocalecp" "1.632289072030 0.019602324152")
LIST(APPEND Al-1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS "ionion" "-10.812685912000 0.000000000000")
  QMC_RUN_AND_CHECK(long-Al-1x1x1-Gaussian-pp-cplx
                    "${CMAKE_SOURCE_DIR}/tests/solids/Al-1x1x1-Gaussian_pp_ShiftedTwist"
                    vmc_long 
                    vmc_long.in.xml
                    1 16
                    TRUE
                    0 Al-1x1x1-Gaussian_pp_Tw_cplx_LONG_SCALARS #VMC 
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
