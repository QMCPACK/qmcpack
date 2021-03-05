
 IF (NOT QMC_CUDA)

# VMC long run with no Jastrows 
# Compare directly to sCI Energy of: E=-10.39616312
LIST(APPEND diamondC_1x1x1-Gaussian_pp_MSD_LONG_SCALARS "totenergy" "-10.39616312 0.03453829914")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_MSD_LONG_SCALARS "kinetic" "11.17885390 0.1357441899")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_MSD_LONG_SCALARS "potential" "-21.59361551 0.1469985177")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_MSD_LONG_SCALARS "eeenergy" "-2.53195870 0.04706895858")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_MSD_LONG_SCALARS "localecp" "-7.14187724 0.17959301094")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_MSD_LONG_SCALARS "nonlocalecp" "0.85588797 0.0446549163")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_MSD_LONG_SCALARS "ionion" "-12.77566754 0.000001")
  QMC_RUN_AND_CHECK(long-diamondC_1x1x1-Gaussian_pp_MSD
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_MSD"
                    vmc_long  
                    vmc_long.in.xml
                    1 16
                    TRUE
                    0 diamondC_1x1x1-Gaussian_pp_MSD_LONG_SCALARS #VMC 
                    )

# VMC short run with no Jastrows 
# Compare directly to sCI Energy of: E=--10.39616312 
LIST(APPEND diamondC_1x1x1-Gaussian_pp_MSD_SHORT_SCALARS "totenergy" "-10.39616312 0.0772197745")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_MSD_SHORT_SCALARS "kinetic" "11.17885390 0.3034931075")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_MSD_SHORT_SCALARS "potential" "-21.59361551 0.3286552225")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_MSD_SHORT_SCALARS "eeenergy" "-2.53195870 0.1052354765")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_MSD_SHORT_SCALARS "localecp" "-7.14187724 0.4015290895")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_MSD_SHORT_SCALARS "nonlocalecp" "0.85588797 0.0998382275")
LIST(APPEND diamondC_1x1x1-Gaussian_pp_MSD_SHORT_SCALARS "ionion" "-12.77566754 0.000001")
  QMC_RUN_AND_CHECK(short-diamondC_1x1x1-Gaussian_pp_MSD
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_MSD"
                    vmc_short 
                    vmc_short.in.xml
                    1 16
                    TRUE
                    0 diamondC_1x1x1-Gaussian_pp_MSD_SHORT_SCALARS #VMC 
                    )


#Deterministic test
IF (QMC_MIXED_PRECISION)
 #VMC
 LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "totenergy" "-10.62293537 0.00002")
 LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "kinetic" "12.00640791 0.00002")
 LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "potential" "-22.62934328 0.00001")
 LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "eeenergy" "-2.50996032 0.000002")
 LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "localecp" "-7.67975015 0.00001")
 LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "nonlocalecp" "0.33603735 0.000005")
 LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "ionion" "-12.77567016 0.000001")
 IF (QMC_COMPLEX)
  #DMC
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-10.99915119 0.00002")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "11.94458085 0.00002")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-22.94373204 0.000003")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-2.97216487 0.000002")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-6.84338204 0.000002")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "-0.35251498 0.000001")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-12.77567016 0.000001")
 ELSE()
  #DMC
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-10.63742696 0.00002")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "12.07373931 0.00002")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-22.71116627 0.000003")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-2.72638743 0.000002")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-7.91852665 0.000002")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "0.70941796 0.000001")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-12.77567016 0.000001")
 ENDIF()
ELSE()
 #VMC
 LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "totenergy" "-10.62293307 0.000001")
 LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "kinetic" "12.00640819 0.000001")
 LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "potential" "-22.62934126 0.000001")
 LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "eeenergy" "-2.50995954 0.000001")
 LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "localecp" "-7.67975160 0.000001")
 LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "nonlocalecp" "0.33603731 0.000001")
 LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "ionion" "-12.77566743 0.000001")
 IF (QMC_COMPLEX)
  #DMC
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-10.99915141 0.000001")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "11.94457971 0.000001")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-22.94373111 0.000001")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-2.97216524 0.000001")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-6.84338340 0.000001")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "-0.35251504 0.000001")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-12.77566743 0.000001")
  ELSE()
  #DMC
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-10.63742384 0.000001")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "12.07372470 0.000001")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-22.71114854 0.000001")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-2.72638698 0.000001")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-7.91851937 0.000001")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "0.70942524 0.000001")
  LIST(APPEND det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-12.77566743 0.000001")
  ENDIF()
ENDIF()

  QMC_RUN_AND_CHECK(deterministic-diamondC_1x1x1-Gaussian_pp_MSD
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp_MSD"
                    det_short_vmc_dmc
                    det_short_vmc_dmc.in.xml
                    1 1
                    TRUE
                    0 det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS #VMC
                    1 det_diamondC_1x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS #DMC
                    )

  ELSE()
    MESSAGE_VERBOSE("Skipping Periodic LCAO as not supported by CUDA build (QMC_CUDA=1)")
  ENDIF()
