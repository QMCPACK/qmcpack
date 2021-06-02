
IF (NOT QMC_CUDA)

# VMC long run with no Jastrows 
# Compare directly to sCI Energy of: E=-21.34158422, T=19.99861447, and V=-41.34019869
# The other references are from VMC reference runs.
LIST(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "totenergy" "-21.34158422 0.0743588132")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "kinetic" "19.99861447 0.0655952484")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "potential" "-41.34019869 0.0743588132")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "eeenergy" "-4.48998595 0.0218346836")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "localecp" "-13.06854246 0.0827662144")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "nonlocalecp" "2.05740964 0.0246666124")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "ionion" "-25.55132719 0.000001")
  QMC_RUN_AND_CHECK(long-diamondC_2x1x1-Gaussian_pp_MSD
                    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_MSD"
                    vmc_long  
                    vmc_long.in.xml
                    1 16
                    TRUE
                    0 diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS #VMC 
                    )

# VMC short run with no Jastrows 
# Compare directly to sCI Energy of: E=-21.34158422, T=19.99861447, and V=-41.34019869
# The other references are from VMC reference runs.
LIST(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "totenergy" "-21.34158422 0.234973849712")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "kinetic" "19.99861447 0.207280984944")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "potential" "-41.34019869 0.234973849712")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "eeenergy" "-4.48998595 0.068997600176")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "localecp" "-13.06854246 0.077946495184")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "nonlocalecp" "2.05740964 0.077946495184")
LIST(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "ionion" "-25.55132719 0.000001")
  QMC_RUN_AND_CHECK(short-diamondC_2x1x1-Gaussian_pp_MSD
                    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_MSD"
                    vmc_short 
                    vmc_short.in.xml
                    1 16
                    TRUE
                    0 diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS #VMC 
                    )

# Deterministic test
IF (QMC_MIXED_PRECISION)
 IF (QMC_COMPLEX)
  # VMC
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "totenergy" "-21.36457073 0.0005")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "kinetic" "18.75860361 0.0005")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "potential" "-40.12317435 0.0001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "eeenergy" "-4.77761915 0.00003")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "localecp" "-10.90122112 0.0002")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "nonlocalecp" "1.10699937 0.0001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "ionion" "-25.55133345 0.000002")
  #DMC
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-21.69556787 0.0005")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "18.36514441 0.0005")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-40.06071229 0.0001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-5.32649930 0.00003")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-10.70569492 0.0002")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "1.52281538 0.0001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-25.55133345 0.000002")
 ELSE()
  # VMC
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "totenergy" "-21.00670927 0.0005")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "kinetic" "20.84452705 0.0005")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "potential" "-41.85123632 0.0001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "eeenergy" "-4.24910344 0.00003")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "localecp" "-13.65564795 0.0002")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "nonlocalecp" "1.60484851 0.0001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "ionion" "-25.55133345 0.000002")
  #DMC
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-21.62353122 0.0005")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "20.66684133 0.0005")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-42.29037255 0.0001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-4.94421717 0.00003")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-13.24215237 0.0002")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "1.44733044 0.0001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-25.55133345 0.000002")
 ENDIF()
ELSE()
 IF (QMC_COMPLEX)
  # VMC
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "totenergy" "-21.36456878 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "kinetic" "18.75860076 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "potential" "-40.12316954 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "eeenergy" "-4.77761896 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "localecp" "-10.90122054 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "nonlocalecp" "1.10699716 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "ionion" "-25.55132719 0.000001")
  #DMC
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-21.69556672 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "18.36514388 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-40.06071060 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-5.32649923 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-10.70569785 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "1.52281367 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-25.55132719 0.000001")
 ELSE()
  # VMC
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "totenergy" "-21.03869868 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "kinetic" "20.61455762 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "potential" "-41.65325630 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "eeenergy" "-4.40963853 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "localecp" "-13.17860152 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "nonlocalecp" "1.48631094 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "ionion" "-25.55132719 0.000001")
  #DMC
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-21.21832528 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "18.98312444 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-40.20144972 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-4.75574225 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-10.90139086 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "1.00701057 0.000001")
  LIST(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-25.55132719 0.000001")
 ENDIF()
ENDIF()

  QMC_RUN_AND_CHECK(deterministic-diamondC_2x1x1-Gaussian_pp_MSD
                    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_MSD"
                    det_vmc_dmc
                    det_vmc_dmc.in.xml
                    1 1
                    TRUE
                    0 det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS #VMC
                    1 det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS #DMC
                    )


ELSE()
  MESSAGE_VERBOSE("Skipping Periodic LCAO as not supported by CUDA build (QMC_CUDA=1)")
ENDIF()
