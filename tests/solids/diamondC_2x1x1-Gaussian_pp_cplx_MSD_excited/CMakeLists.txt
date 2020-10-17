
 IF (NOT QMC_CUDA)
   IF (QMC_COMPLEX)

# VMC short run for 1st excited state with no Jastrows 
# Compare directly to sCI Energy of: E=-20.6068860259079
#  LocalEnergy           =         -20.55017 +/-          0.00093
#  Variance              =             4.582 +/-            0.017
#  Kinetic               =           21.0864 +/-           0.0026
#  LocalPotential        =          -41.6366 +/-           0.0027
#  ElecElec              =           -4.4376 +/-           0.0011
#  LocalECP              =          -13.7431 +/-           0.0033
#  NonLocalECP           =            2.0955 +/-           0.0012
#  IonIon                =            -25.55 +/-             0.00
#  TotalSamples          =          12226560 +/-                0
#Data generated using utils/make_ref_data.sh  diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_SHORT qmc_ref/vmc_short_excited_noj.s000.scalar.dat qmc_ref/qmc_ref_excited_long_noj.s000.scalar.dat
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_SHORT_SCALARS "totenergy" "-20.550276480536 0.039294710408")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_SHORT_SCALARS "variance" "4.581714552801 0.137321481608")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_SHORT_SCALARS "kinetic" "21.086544715418 0.177460056409")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_SHORT_SCALARS "potential" "-41.636821195959 0.185935091698")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_SHORT_SCALARS "eeenergy" "-4.437697375480 0.048410609989")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_SHORT_SCALARS "localecp" "-13.743308395587 0.225541823772")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_SHORT_SCALARS "nonlocalecp" "2.095511637186 0.061381665845")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_SHORT_SCALARS "ionion" "-25.551327062000 0.001000000000")
  QMC_RUN_AND_CHECK(short-diamondC_2x1x1-Gaussian_pp_cplx_excited_noj_MSD
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_cplx_MSD_excited"
                    vmc_short_excited_noj
                    vmc_short_excited_noj.in.xml
                    1 16
                    TRUE
                    0 diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_SHORT_SCALARS #VMC 
                    )

# VMC long run for 1st excited state with no Jastrows 
# Compare directly to sCI Energy of: E=-20.6068860259079
#  LocalEnergy           =         -20.55017 +/-          0.00093
#  Variance              =             4.582 +/-            0.017
#  Kinetic               =           21.0864 +/-           0.0026
#  LocalPotential        =          -41.6366 +/-           0.0027
#  ElecElec              =           -4.4376 +/-           0.0011
#  LocalECP              =          -13.7431 +/-           0.0033
#  NonLocalECP           =            2.0955 +/-           0.0012
#  IonIon                =            -25.55 +/-             0.00
#  TotalSamples          =          12226560 +/-                0
#Data generated using utils/make_ref_data.sh  diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_LONG qmc_ref/vmc_long_excited_noj.s000.scalar.dat qmc_ref/qmc_ref_excited_long_noj.s000.scalar.dat
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_LONG_SCALARS "totenergy" "-20.550276480536 0.017932198313")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_LONG_SCALARS "variance" "4.581714552801 0.153381113169")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_LONG_SCALARS "kinetic" "21.086544715418 0.059079832112")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_LONG_SCALARS "potential" "-41.636821195959 0.061592666021")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_LONG_SCALARS "eeenergy" "-4.437697375480 0.020549201829")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_LONG_SCALARS "localecp" "-13.743308395587 0.077797562873")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_LONG_SCALARS "nonlocalecp" "2.095511637186 0.024858616101")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_LONG_SCALARS "ionion" "-25.551327062000 0.001000000000")
  QMC_RUN_AND_CHECK(long-diamondC_2x1x1-Gaussian_pp_cplx_excited_noj_MSD
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_cplx_MSD_excited"
                    vmc_long_excited_noj
                    vmc_long_excited_noj.in.xml
                    1 16
                    TRUE
                    0 diamondC_2x1x1-Gaussian_pp_cplx_MSD_Excited_NOJ_LONG_SCALARS #VMC 
                    )


# VMC short run for ground state with Jastrows 
# Compare directly to sCI Energy of: E=-21.3957193782765
#Data generated using utils/make_ref_data.sh  diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT qmc_ref/vmc_short_ground.s000.scalar.dat qmc_ref/qmc_ref_ground_long.s000.scalar.dat
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS "totenergy" "-21.821112814361 0.022782931094")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS "variance" "1.051140170460 0.052716149007")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS "kinetic" "20.090944056549 0.145500025188")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS "potential" "-41.912056870903 0.152689438716")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS "eeenergy" "-5.471182382216 0.048835207119")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS "localecp" "-12.996437953778 0.214186877447")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS "nonlocalecp" "2.106890527097 0.064007725383")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS "ionion" "-25.551327062000 0.001000000000")
  QMC_RUN_AND_CHECK(short-diamondC_2x1x1-Gaussian_pp_cplx_ground_MSD                                             
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_cplx_MSD_excited"                    
                    vmc_short_ground
                    vmc_short_ground.in.xml 
                    1 16
                    TRUE
                    0 diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_SHORT_SCALARS #VMC                              
                    )
      

# VMC long  run for ground state with Jastrows 
# Compare directly to sCI Energy of: E=-21.3957193782765
#Data generated using utils/make_ref_data.sh  diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG qmc_ref/vmc_long_ground.s000.scalar.dat qmc_ref/qmc_ref_ground_long.s000.scalar.dat
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS "totenergy" "-21.821112814361 0.006665434144")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS "variance" "1.051140170460 0.024813291267")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS "kinetic" "20.090944056549 0.066235207833")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS "potential" "-41.912056870903 0.068349719327")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS "eeenergy" "-5.471182382216 0.016467750812")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS "localecp" "-12.996437953778 0.076427146609")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS "nonlocalecp" "2.106890527097 0.022662692528")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS "ionion" "-25.551327062000 0.001000000000")
  QMC_RUN_AND_CHECK(long-diamondC_2x1x1-Gaussian_pp_cplx_ground_MSD
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_cplx_MSD_excited"                    
                    vmc_long_ground
                    vmc_long_ground.in.xml 
                   1 16
                    TRUE
                    0 diamondC_2x1x1-Gaussian_pp_cplx_MSD_GROUND_LONG_SCALARS #VMC                              
                    )
      

# VMC short run for 1st excited state with Jastrows 
# Compare directly to sCI Energy of: E=-20.6068860259079
#  LocalEnergy           =         -21.04339 +/-          0.00045
#  Variance              =           1.44953 +/-          0.00094
#  Kinetic               =           21.0062 +/-           0.0027
#  LocalPotential        =          -42.0496 +/-           0.0028
#  ElecElec              =          -5.42851 +/-          0.00076
#  LocalECP              =          -13.0542 +/-           0.0031
#  NonLocalECP           =            1.9844 +/-           0.0010
#  IonIon                =            -25.55 +/-             0.00
#  TotalSamples          =           9912320 +/-                0
#Data generated using utils/make_ref_data.sh  diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT qmc_ref/vmc_short_excited.s000.scalar.dat qmc_ref/qmc_ref_excited_long.s000.scalar.dat
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS "totenergy" "-21.043422517714 0.026265551021")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS "variance" "1.449654318091 0.052594664103")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS "kinetic" "21.005794394605 0.133579143250")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS "potential" "-42.049216912354 0.141199415636")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS "eeenergy" "-5.428377875848 0.066531858421")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS "localecp" "-13.053874507313 0.175835920111")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS "nonlocalecp" "1.984362532886 0.060268377421")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS "ionion" "-25.551327062000 0.001000000000")
  QMC_RUN_AND_CHECK(short-diamondC_2x1x1-Gaussian_pp_cplx_excited_MSD                                             
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_cplx_MSD_excited"                    
                    vmc_short_excited
                    vmc_short_excited.in.xml 
                    1 16
                    TRUE
                    0 diamondC_2x1x1-Gaussian_pp_cplx_MSD_EXCITED_SHORT_SCALARS #VMC                              
                    )
      

# VMC long  run for 1st excited state with Jastrows 
# Compare directly to sCI Energy of: E=-20.6068860259079
#  LocalEnergy           =         -21.04339 +/-          0.00045
#  Variance              =           1.44953 +/-          0.00094
#  Kinetic               =           21.0062 +/-           0.0027
#  LocalPotential        =          -42.0496 +/-           0.0028
#  ElecElec              =          -5.42851 +/-          0.00076
#  LocalECP              =          -13.0542 +/-           0.0031
#  NonLocalECP           =            1.9844 +/-           0.0010
#  IonIon                =            -25.55 +/-             0.00
#  TotalSamples          =           9912320 +/-                0
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
