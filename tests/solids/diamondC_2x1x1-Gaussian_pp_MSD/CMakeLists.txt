# VMC long run with no Jastrows
# Compare directly to sCI Energy of: E=-21.34158422, T=19.99861447, and V=-41.34019869
# The other references are from VMC reference runs.
list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "totenergy" "-21.34158422 0.0743588132")
list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "kinetic" "19.99861447 0.0655952484")
list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "potential" "-41.34019869 0.0743588132")
list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "eeenergy" "-4.48998595 0.0218346836")
list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "localecp" "-13.06854246 0.0827662144")
list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "nonlocalecp" "2.05740964 0.0246666124")
list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "ionion" "-25.55132719 0.000003")
qmc_run_and_check(
  long-diamondC_2x1x1-Gaussian_pp_MSD
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_MSD"
  vmc_long
  vmc_long.in.xml
  1
  16
  TRUE
  0
  diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS #VMC
)

# VMC short run with no Jastrows
# Compare directly to sCI Energy of: E=-21.34158422, T=19.99861447, and V=-41.34019869
# The other references are from VMC reference runs.
list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "totenergy" "-21.34158422 0.234973849712")
list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "kinetic" "19.99861447 0.207280984944")
list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "potential" "-41.34019869 0.234973849712")
list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "eeenergy" "-4.48998595 0.1")
list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "localecp" "-13.06854246 0.2")
list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "nonlocalecp" "2.05740964 0.077946495184")
list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "ionion" "-25.55132719 0.000003")
qmc_run_and_check(
  short-diamondC_2x1x1-Gaussian_pp_MSD
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_MSD"
  vmc_short
  vmc_short.in.xml
  1
  16
  TRUE
  0
  diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS #VMC
)

# Deterministic test
if(QMC_MIXED_PRECISION)
  if(QMC_COMPLEX)
    # VMC
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "totenergy" "-21.36457073 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "kinetic" "18.75860361 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "potential" "-40.12317435 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "eeenergy" "-4.77761915 0.00003")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "localecp" "-10.90122112 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "nonlocalecp" "1.10699937 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "ionion" "-25.55132413 0.000002")
    #DMC
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-21.69556787 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "18.36514441 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-40.06071229 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-5.32649930 0.00003")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-10.70569492 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "1.52281538 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-25.55132413 0.000002")
  else()
    # VMC
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "totenergy" "-21.00670927 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "kinetic" "20.84452705 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "potential" "-41.85123632 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "eeenergy" "-4.24910344 0.00003")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "localecp" "-13.65564795 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "nonlocalecp" "1.60484851 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "ionion" "-25.55132413 0.000002")
    #DMC
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-21.62353122 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "20.66684133 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-42.29037255 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-4.94421717 0.00003")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-13.24215237 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "1.44733044 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-25.55132413 0.000002")
  endif()
else()
  if(QMC_COMPLEX)
    # VMC
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "totenergy" "-21.36458168 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "kinetic" "18.75860076 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "potential" "-40.12318244 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "eeenergy" "-4.77761896 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "localecp" "-10.90122054 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "nonlocalecp" "1.10698394 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "ionion" "-25.55132719 0.000001")
    #DMC
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-21.69562984 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "18.36514388 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-40.06077372 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-5.32649923 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-10.70569785 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "1.52275023 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-25.55132719 0.000001")
  else()
    # VMC
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "totenergy" "-21.03866915 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "kinetic" "20.61455762 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "potential" "-41.65322677 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "eeenergy" "-4.40963853 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "localecp" "-13.17860152 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "nonlocalecp" "1.48634016 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "ionion" "-25.55132719 0.000001")
    #DMC
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-21.21829172 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "18.98312444 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-40.20141616 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-4.75574225 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-10.90139086 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "1.00704382 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-25.55132719 0.000001")
  endif()
endif()

qmc_run_and_check(
  deterministic-diamondC_2x1x1-Gaussian_pp_MSD
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_MSD"
  det_vmc_dmc
  det_vmc_dmc.in.xml
  1
  1
  TRUE
  0
  det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS #VMC
  1
  det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS #DMC
)

if(QMC_MIXED_PRECISION)
  if(QMC_COMPLEX)
    # VMC without drift
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "totenergy"
         "-21.03266507 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "kinetic"
         "21.53043728 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "potential"
         "-42.56310227 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "eeenergy"
         "-4.20408649 0.00003")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "ionion"
         "-25.55132413 0.000002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "localecp"
         "-14.93105730 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "nonlocalecp"
         "2.12337515 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "samples" "144 0.0")
    # VMC with drift
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "totenergy"
         "-21.79244047 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "kinetic"
         "20.34265653 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "potential"
         "-42.13509704 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "eeenergy"
         "-4.99765058 0.00003")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "ionion"
         "-25.55132413 0.000002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "localecp"
         "-14.00169909 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "nonlocalecp"
         "2.41558627 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "samples" "144 0.0")
    # DMC
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "totenergy"
         "-21.77438963 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "kinetic"
         "21.72919297 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "potential"
         "-43.50358262 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "eeenergy"
         "-4.47590829 0.00003")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "ionion"
         "-25.55132413 0.000002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "localecp"
         "-15.03973198 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "nonlocalecp"
         "1.56339129 0.0001")
  else()
    # VMC without drift
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "totenergy"
         "-21.07292465 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "kinetic"
         "20.30465866 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "potential"
         "-41.37758333 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "eeenergy"
         "-4.43122622 0.00003")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "ionion"
         "-25.55132413 0.000002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "localecp"
         "-12.89019745 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "nonlocalecp"
         "1.49517397 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "samples" "144.00000000 0.0")
    # VMC with drift
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "totenergy"
         "-21.01757324 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "kinetic"
         "20.30786316 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "potential"
         "-41.32543639 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "eeenergy"
         "-4.34588506 0.00003")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "ionion"
         "-25.55132413 0.000002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "localecp"
         "-13.94038870 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "nonlocalecp"
         "2.51217101 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "samples" "144.00000000 0.0")
    # DMC
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "totenergy"
         "-20.84853817 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "kinetic"
         "20.55342353 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "potential"
         "-41.40196171 0.0001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "eeenergy"
         "-4.03628785 0.00003")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "ionion"
         "-25.55132413 0.000002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "localecp"
         "-14.32128624 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "nonlocalecp"
         "2.50694602 0.0001")
  endif()
else()
  if(QMC_COMPLEX)
    # VMC without drift
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "totenergy"
         "-21.03263547 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "kinetic"
         "21.53043534 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "potential"
         "-42.56307081 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "eeenergy"
         "-4.20408808 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "ionion"
         "-25.55132706 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "localecp"
         "-14.93105741 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "nonlocalecp"
         "2.12340158 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "samples" "144 0.0")
    # VMC with drift
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "totenergy"
         "-21.79249332 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "kinetic"
         "20.34264322 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "potential"
         "-42.13513654 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "eeenergy"
         "-4.99765165 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "ionion"
         "-25.55132706 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "localecp"
         "-14.00168724 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "nonlocalecp"
         "2.41552925 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "samples" "144 0.0")
    # DMC
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "totenergy"
         "-21.77441682 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "kinetic"
         "21.72916062 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "potential"
         "-43.50357863 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "eeenergy"
         "-4.47591187 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "ionion"
         "-25.55132706 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "localecp"
         "-15.03968077 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "nonlocalecp"
         "1.56334141 0.000001")
  else()
    # VMC without drift
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "totenergy"
         "-20.52919488 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "kinetic"
         "20.07128892 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "potential"
         "-40.60048380 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "eeenergy"
         "-3.93794601 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "ionion"
         "-25.55132719 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "localecp"
         "-12.96301894 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "nonlocalecp"
         "1.85180834 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "samples" "144.00000000 0.0")
    # VMC with drift
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "totenergy"
         "-21.16501162 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "kinetic"
         "19.67417373 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "potential"
         "-40.83918534 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "eeenergy"
         "-4.53265658 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "ionion"
         "-25.55132719 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "localecp"
         "-13.36631725 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "nonlocalecp"
         "2.61111536 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "samples" "144.00000000 0.0")
    # DMC
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "totenergy"
         "-20.79231965 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "kinetic"
         "21.05924365 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "potential"
         "-41.85156470 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "eeenergy"
         "-3.83828756 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "ionion"
         "-25.55132719 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "localecp"
         "-16.19033691 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "nonlocalecp"
         "3.72838694 0.000001")
  endif()
endif()

qmc_run_and_check(
  deterministic-diamondC_2x1x1-gaussian_pp_MSD-vmcbatch-dmcbatch-mwalkers_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_MSD"
  det_qmc_vmcbatch_dmcbatch_mwalkers
  det_qmc_vmcbatch_dmcbatch_mwalkers.in.xml
  1
  4
  TRUE
  0
  det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS # VMC without drift
  1
  det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS # VMC with drift
  2
  det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS # DMC
)

if(QMC_MIXED_PRECISION)
  if(QMC_COMPLEX)
    # VMC without drift
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "totenergy"
         "-21.03294851 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "kinetic"
         "20.65585656 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "potential"
         "-41.68880500 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "eeenergy"
         "-4.40762368 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "ionion"
         "-25.55132413 0.000002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "localecp"
         "-13.87483391 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "nonlocalecp"
         "2.14498622 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "samples" "144 0.0")
    # VMC with drift
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "totenergy"
         "-21.51006228 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "kinetic"
         "19.78928388 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "potential"
         "-41.29934624 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "eeenergy"
         "-4.87340890 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "ionion"
         "-25.55132413 0.000002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "localecp"
         "-13.00551457 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "nonlocalecp"
         "2.13091087 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "samples" "144 0.0")
    # DMC
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "totenergy"
         "-21.60474948 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "kinetic"
         "21.43358030 0.0010")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "potential"
         "-43.03832971 0.0010")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "eeenergy"
         "-4.62211747 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "ionion"
         "-25.55132413 0.000002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "localecp"
         "-14.01500084 0.0004")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "nonlocalecp"
         "1.15012223 0.0010")
  else()
    # VMC without drift
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "totenergy"
         "-20.92212995 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "kinetic"
         "19.65122803 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "potential"
         "-40.57335794 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "eeenergy"
         "-4.44192720 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "ionion"
         "-25.55132413 0.000002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "localecp"
         "-11.69255184 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "nonlocalecp"
         "1.11245474 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "samples" "144.00000000 0.0")
    # VMC with drift
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "totenergy"
         "-20.81739152 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "kinetic"
         "20.27356775 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "potential"
         "-41.09095922 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "eeenergy"
         "-4.23387512 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "ionion"
         "-25.55132413 0.000002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "localecp"
         "-13.74506329 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "nonlocalecp"
         "2.43931282 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "samples" "144.00000000 0.0")
    # DMC
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "totenergy"
         "-20.52433054 0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "kinetic"
         "20.96842887 0.0010")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "potential"
         "-41.49275929 0.0010")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "eeenergy"
         "-3.69096631 0.0002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "ionion"
         "-25.55132413 0.000002")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "localecp"
         "-15.85799109  0.0005")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "nonlocalecp"
         "3.60753175 0.0010")
  endif()
else()
  if(QMC_COMPLEX)
    # VMC without drift
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "totenergy"
         "-21.03298214 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "kinetic"
         "20.65585352 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "potential"
         "-41.68883566 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "eeenergy"
         "-4.40762168 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "ionion"
         "-25.55132706 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "localecp"
         "-13.87483540 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "nonlocalecp"
         "2.14494833 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "samples" "144 0.0")
    # VMC with drift
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "totenergy"
         "-21.51010377 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "kinetic"
         "19.78926663 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "potential"
         "-41.29937040 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "eeenergy"
         "-4.87341474 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "ionion"
         "-25.55132706 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "localecp"
         "-13.00549989 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "nonlocalecp"
         "2.13087113 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "samples" "144 0.0")
    # DMC
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "totenergy"
         "-21.60481550 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "kinetic"
         "21.43354420 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "potential"
         "-43.03836124 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "eeenergy"
         "-4.62212905 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "ionion"
         "-25.55132706 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "localecp"
         "-14.01493383 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "nonlocalecp"
         "1.15003007 0.000001")
  else()
    # VMC without drift
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "totenergy"
         "-20.84291775 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "kinetic"
         "19.68818617 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "potential"
         "-40.53110391 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "eeenergy"
         "-4.27768409 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "ionion"
         "-25.55132706 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "localecp"
         "-11.80294323 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "nonlocalecp"
         "1.10085031 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "samples" "144.00000000 0.0")
    # VMC with drift
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "totenergy"
         "-20.97512384 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "kinetic"
         "19.90145613 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "potential"
         "-40.87657996 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "eeenergy"
         "-4.46662861 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "ionion"
         "-25.55132706 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "localecp"
         "-13.57462590 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "nonlocalecp"
         "2.71600145 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "samples" "144.00000000 0.0")
    # DMC
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "totenergy"
         "-21.06958081 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "kinetic"
         "21.41422718 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "potential"
         "-42.48381031 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "eeenergy"
         "-4.24343589 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "ionion"
         "-25.55132706 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "localecp"
         "-16.56049089 0.000001")
    list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "nonlocalecp"
         "3.87144211 0.000001")
  endif()
endif()

qmc_run_and_check(
  deterministic-diamondC_2x1x1-gaussian_pp_MSD-vmcbatch-dmcbatch-mwalkers_sdj
  "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_MSD"
  det_qmc_vmcbatch_dmcbatch_mwalkers
  det_qmc_vmcbatch_dmcbatch_mwalkers.in.xml
  2
  2
  TRUE
  0
  det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS # VMC without drift
  1
  det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS # VMC with drift
  2
  det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS # DMC
)
