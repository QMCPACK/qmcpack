if(NOT QMC_CUDA)

  # VMC long run with no Jastrows
  # Compare directly to sCI Energy of: E=-21.34158422, T=19.99861447, and V=-41.34019869
  # The other references are from VMC reference runs.
  list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "totenergy" "-21.34158422 0.0743588132")
  list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "kinetic" "19.99861447 0.0655952484")
  list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "potential" "-41.34019869 0.0743588132")
  list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "eeenergy" "-4.48998595 0.0218346836")
  list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "localecp" "-13.06854246 0.0827662144")
  list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "nonlocalecp" "2.05740964 0.0246666124")
  list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_LONG_SCALARS "ionion" "-25.55132719 0.000001")
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
  list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "eeenergy" "-4.48998595 0.068997600176")
  list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "localecp" "-13.06854246 0.077946495184")
  list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "nonlocalecp" "2.05740964 0.077946495184")
  list(APPEND diamondC_2x1x1-Gaussian_pp_MSD_SHORT_SCALARS "ionion" "-25.55132719 0.000001")
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
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "ionion" "-25.55133345 0.000002")
      #DMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-21.69556787 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "18.36514441 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-40.06071229 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-5.32649930 0.00003")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-10.70569492 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "1.52281538 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-25.55133345 0.000002")
    else()
      # VMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "totenergy" "-21.00670927 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "kinetic" "20.84452705 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "potential" "-41.85123632 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "eeenergy" "-4.24910344 0.00003")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "localecp" "-13.65564795 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "nonlocalecp" "1.60484851 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "ionion" "-25.55133345 0.000002")
      #DMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-21.62353122 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "20.66684133 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-42.29037255 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-4.94421717 0.00003")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-13.24215237 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "1.44733044 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-25.55133345 0.000002")
    endif()
  else()
    if(QMC_COMPLEX)
      # VMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "totenergy" "-21.36456878 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "kinetic" "18.75860076 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "potential" "-40.12316954 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "eeenergy" "-4.77761896 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "localecp" "-10.90122054 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "nonlocalecp" "1.10699716 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "ionion" "-25.55132719 0.000001")
      #DMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-21.69556672 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "18.36514388 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-40.06071060 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-5.32649923 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-10.70569785 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "1.52281367 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "ionion" "-25.55132719 0.000001")
    else()
      # VMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "totenergy" "-21.03869868 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "kinetic" "20.61455762 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "potential" "-41.65325630 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "eeenergy" "-4.40963853 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "localecp" "-13.17860152 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "nonlocalecp" "1.48631094 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_VMC_SCALARS "ionion" "-25.55132719 0.000001")
      #DMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "totenergy" "-21.21832528 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "kinetic" "18.98312444 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "potential" "-40.20144972 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "eeenergy" "-4.75574225 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "localecp" "-10.90139086 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_SHORT_DMC_SCALARS "nonlocalecp" "1.00701057 0.000001")
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
           "-25.55133363 0.000002")
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
           "-25.55133363 0.000002")
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
           "-25.55133363 0.000002")
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
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "localecp"
           "-12.89019745 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "nonlocalecp"
           "1.49517397 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "samples"
           "144.00000000 0.0")
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
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "localecp"
           "-13.94038870 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "nonlocalecp"
           "2.51217101 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "samples"
           "144.00000000 0.0")
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
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "localecp"
           "-14.32128624 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "nonlocalecp"
           "2.50694602 0.0001")
    endif()
  else()
    if(QMC_COMPLEX)
      # VMC without drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "totenergy"
           "-21.03266250 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "kinetic"
           "21.53043534 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "potential"
           "-42.56309784 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "eeenergy"
           "-4.20408808 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "localecp"
           "-14.93105741 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "nonlocalecp"
           "2.12337471 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "samples" "144 0.0")
      # VMC with drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "totenergy"
           "-21.79244251 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "kinetic"
           "20.34264322 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "potential"
           "-42.13508573 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "eeenergy"
           "-4.99765165 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "localecp"
           "-14.00168724 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "nonlocalecp"
           "2.41558022 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "samples" "144 0.0")
      # DMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "totenergy"
           "-21.77438877 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "kinetic"
           "21.72916062 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "potential"
           "-43.50354939 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "eeenergy"
           "-4.47591187 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "localecp"
           "-15.03968077 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "nonlocalecp"
           "1.56337031 0.000001")
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
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_1_SCALARS "samples"
           "144.00000000 0.0")
      # VMC with drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "totenergy"
           "-21.16496822 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "kinetic"
           "19.67417373 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "potential"
           "-40.83914194 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "eeenergy"
           "-4.53265658 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "ionion"
           "-25.55132719 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "localecp"
           "-13.36631725 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "nonlocalecp"
           "2.61115908 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_2_SCALARS "samples"
           "144.00000000 0.0")
      # DMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "totenergy"
           "-20.79227061 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "kinetic"
           "21.05924365 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "potential"
           "-41.85151426 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "eeenergy"
           "-3.83828756 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "ionion"
           "-25.55132719 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "localecp"
           "-16.19033691 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t4_3_SCALARS "nonlocalecp"
           "3.72843740 0.000001")
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
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "totenergy"
           "-21.44268115 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "kinetic"
           "17.50276721 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "potential"
           "-38.94544861 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "eeenergy"
           "-5.08527221 0.00003")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "ionion"
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "localecp"
           "-9.83228444 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "nonlocalecp"
           "1.52344166 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "samples" "144 0.0")
      # VMC with drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "totenergy"
           "-21.28433738 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "kinetic"
           "21.27983264 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "potential"
           "-42.56417009 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "eeenergy"
           "-4.68529565 0.00003")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "ionion"
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "localecp"
           "-15.31087459 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "nonlocalecp"
           "2.98333379 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "samples" "144 0.0")
      # DMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "totenergy"
           "-22.02240256 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "kinetic"
           "22.26955921 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "potential"
           "-44.29196173 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "eeenergy"
           "-4.92804498 0.00003")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "ionion"
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "localecp"
           "-17.22154454 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "nonlocalecp"
           "3.40896143 0.0001")
    else()
      # VMC without drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "totenergy"
           "-21.24377773 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "kinetic"
           "18.95261961 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "potential"
           "-40.19639737 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "eeenergy"
           "-5.04910503 0.00003")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "ionion"
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "localecp"
           "-10.87790041 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "nonlocalecp"
           "1.28194171 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "samples"
           "144.00000000 0.0")
      # VMC with drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "totenergy"
           "-21.71080331 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "kinetic"
           "20.92941038 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "potential"
           "-42.64021376 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "eeenergy"
           "-4.64881369 0.00003")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "ionion"
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "localecp"
           "-14.32409006 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "nonlocalecp"
           "1.88402361 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "samples"
           "144.00000000 0.0")
      # DMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "totenergy"
           "-21.60271630 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "kinetic"
           "21.47736273 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "potential"
           "-43.08007903 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "eeenergy"
           "-4.62934502 0.00003")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "ionion"
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "localecp"
           "-14.42788710 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "nonlocalecp"
           "1.52848673 0.0001")
    endif()
  else()
    if(QMC_COMPLEX)
      # VMC without drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "totenergy"
           "-21.44267742 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "kinetic"
           "17.50276467 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "potential"
           "-38.94544209 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "eeenergy"
           "-5.08527024 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "localecp"
           "-9.83228352 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "nonlocalecp"
           "1.52343873 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "samples" "144 0.0")
      # VMC with drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "totenergy"
           "-21.28433581 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "kinetic"
           "21.27982716 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "potential"
           "-42.56416297 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "eeenergy"
           "-4.68529657 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "localecp"
           "-15.31087672 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "nonlocalecp"
           "2.98333738 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "samples" "144 0.0")
      # DMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "totenergy"
           "-22.02239785 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "kinetic"
           "22.26956045 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "potential"
           "-44.29195830 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "eeenergy"
           "-4.92804100 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "localecp"
           "-17.22155131 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "nonlocalecp"
           "3.40896107 0.000001")
    else()
      # VMC without drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "totenergy"
           "-21.30267144 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "kinetic"
           "19.93254641 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "potential"
           "-41.23521784 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "eeenergy"
           "-5.06956602 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "ionion"
           "-25.55132719 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "localecp"
           "-12.19151812 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "nonlocalecp"
           "1.57719350 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS "samples"
           "144.00000000 0.0")
      # VMC with drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "totenergy"
           "-21.50423449 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "kinetic"
           "20.70369649 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "potential"
           "-42.20793098 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "eeenergy"
           "-4.69046880 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "ionion"
           "-25.55132719 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "localecp"
           "-13.26890009 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "nonlocalecp"
           "1.30276510 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS "samples"
           "144.00000000 0.0")
      # DMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "totenergy"
           "-21.64073241 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "kinetic"
           "20.93567485 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "potential"
           "-42.57640726 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "eeenergy"
           "-4.90231612 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "ionion"
           "-25.55132719 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "localecp"
           "-13.86140151 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS "nonlocalecp"
           "1.73863757 0.000001")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1-gaussian_pp_MSD-vmcbatch-dmcbatch-mwalkers_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_MSD"
    det_qmc_vmcbatch_dmcbatch_mwalkers
    det_qmc_vmcbatch_dmcbatch_mwalkers.in.xml
    1
    16
    TRUE
    0
    det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_1_SCALARS # VMC without drift
    1
    det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_2_SCALARS # VMC with drift
    2
    det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r1-t16_3_SCALARS # DMC
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
           "-25.55133363 0.000002")
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
           "-25.55133363 0.000002")
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
           "-25.55133363 0.000002")
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
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "localecp"
           "-11.69255184 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "nonlocalecp"
           "1.11245474 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "samples"
           "144.00000000 0.0")
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
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "localecp"
           "-13.74506329 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "nonlocalecp"
           "2.43931282 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "samples"
           "144.00000000 0.0")
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
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "localecp"
           "-15.85799109  0.0004")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "nonlocalecp"
           "3.60753175 0.0010")
    endif()
  else()
    if(QMC_COMPLEX)
      # VMC without drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "totenergy"
           "-21.03294400 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "kinetic"
           "20.65585352 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "potential"
           "-41.68879752 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "eeenergy"
           "-4.40762168 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "localecp"
           "-13.87483540 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "nonlocalecp"
           "2.14498663 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "samples" "144 0.0")
      # VMC with drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "totenergy"
           "-21.51006989 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "kinetic"
           "19.78926663 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "potential"
           "-41.29933652 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "eeenergy"
           "-4.87341474 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "localecp"
           "-13.00549989 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "nonlocalecp"
           "2.13090517 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "samples" "144 0.0")
      # DMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "totenergy"
           "-21.60475684 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "kinetic"
           "21.43356346 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "potential"
           "-43.03832030 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "eeenergy"
           "-4.62212422 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "localecp"
           "-14.01496528 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "nonlocalecp"
           "1.15009626 0.000001")
    else()
      # VMC without drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "totenergy"
           "-20.84292853 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "kinetic"
           "19.68818617 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "potential"
           "-40.53111470 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "eeenergy"
           "-4.27768409 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "localecp"
           "-11.80294323 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "nonlocalecp"
           "1.10083968 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_1_SCALARS "samples"
           "144.00000000 0.0")
      # VMC with drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "totenergy"
           "-20.97505088 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "kinetic"
           "19.90145613 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "potential"
           "-40.87650700 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "eeenergy"
           "-4.46662861 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "localecp"
           "-13.57462590 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "nonlocalecp"
           "2.71607457 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_2_SCALARS "samples"
           "144.00000000 0.0")
      # DMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "totenergy"
           "-21.06954729 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "kinetic"
           "21.41422718 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "potential"
           "-42.48377447 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "eeenergy"
           "-4.24343589 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "localecp"
           "-16.56048685 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t2_3_SCALARS "nonlocalecp"
           "3.87147533 0.000001")
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

  if(QMC_MIXED_PRECISION)
    if(QMC_COMPLEX)
      # VMC without drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "totenergy"
           "-21.72802929 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "kinetic"
           "20.66308784 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "potential"
           "-42.39111731 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "eeenergy"
           "-4.90724094 0.00003")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "ionion"
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "localecp"
           "-13.85965673 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "nonlocalecp"
           "1.92711400 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "samples" "144 0.0")
      # VMC with drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "totenergy"
           "-21.59728865 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "kinetic"
           "21.27139123 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "potential"
           "-42.86867993 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "eeenergy"
           "-4.45806315 0.00003")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "ionion"
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "localecp"
           "-14.75288075 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "nonlocalecp"
           "1.89359760 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "samples" "144 0.0")
      # DMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "totenergy"
           "-21.66458378 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "kinetic"
           "20.18933860 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "potential"
           "-41.85392239 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "eeenergy"
           "-4.81842941 0.00003")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "ionion"
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "localecp"
           "-12.98244950 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "nonlocalecp"
           "1.49829015 0.0001")
    else()
      # VMC without drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "totenergy"
           "-21.26278722 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "kinetic"
           "20.70065712 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "potential"
           "-41.96344439 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "eeenergy"
           "-4.45926415 0.00003")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "ionion"
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "localecp"
           "-13.74208188 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "nonlocalecp"
           "1.78923528 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "samples"
           "144.00000000 0.0")
      # VMC with drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "totenergy"
           "-21.44219502 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "kinetic"
           "18.72733136 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "potential"
           "-40.16952630 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "eeenergy"
           "-4.93674832 0.00003")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "ionion"
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "localecp"
           "-11.60307190 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "nonlocalecp"
           "1.92162757 0.0003")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "samples"
           "144.00000000 0.0")
      # DMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "totenergy"
           "-21.74222240 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "kinetic"
           "19.01632430 0.0005")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "potential"
           "-40.75854664 0.0001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "eeenergy"
           "-5.59636306 0.00003")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "ionion"
           "-25.55133363 0.000002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "localecp"
           "-10.52715232 0.0002")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "nonlocalecp"
           "0.91630238 0.0001")
    endif()
  else()
    if(QMC_COMPLEX)
      # VMC without drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "totenergy"
           "-21.58289265 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "kinetic"
           "20.38781180 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "potential"
           "-41.97070445 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "eeenergy"
           "-4.79329398 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "localecp"
           "-13.37215305 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "nonlocalecp"
           "1.74606964 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "samples" "144 0.0")
      # VMC with drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "totenergy"
           "-21.76136769 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "kinetic"
           "21.68396377 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "potential"
           "-43.44533145 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "eeenergy"
           "-4.61354584 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "localecp"
           "-14.61321213 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "nonlocalecp"
           "1.33275358 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "samples" "144 0.0")
      # DMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "totenergy"
           "-21.27076431 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "kinetic"
           "20.20677813 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "potential"
           "-41.47754244 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "eeenergy"
           "-4.42449658 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "localecp"
           "-12.46418943 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "nonlocalecp"
           "0.96247062 0.000001")
    else()
      # VMC without drift
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "totenergy"
           "-21.03372308 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "kinetic"
           "20.68179752 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "potential"
           "-41.71552059 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "eeenergy"
           "-4.05542459 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "localecp"
           "-14.26281752 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "nonlocalecp"
           "2.15404858 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS "samples"
           "144.00000000 0.0")
      # VMC with drif
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "totenergy"
           "-21.54032047 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "kinetic"
           "19.92752205 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "potential"
           "-41.46784252 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "eeenergy"
           "-4.99217592 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "localecp"
           "-12.74451933 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "nonlocalecp"
           "1.82017979 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS "samples"
           "144.00000000 0.0")
      # DMC
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "totenergy"
           "-21.61183432 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "kinetic"
           "20.63430258 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "potential"
           "-42.24613689 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "eeenergy"
           "-5.37797272 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "ionion"
           "-25.55132706 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "localecp"
           "-12.33231856 0.000001")
      list(APPEND det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS "nonlocalecp"
           "1.01548145 0.000001")
    endif()
  endif()

  qmc_run_and_check(
    deterministic-diamondC_2x1x1-gaussian_pp_MSD-vmcbatch-dmcbatch-mwalkers_sdj
    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_2x1x1-Gaussian_pp_MSD"
    det_qmc_vmcbatch_dmcbatch_mwalkers
    det_qmc_vmcbatch_dmcbatch_mwalkers.in.xml
    2
    3
    TRUE
    0
    det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_1_SCALARS # VMC without drift
    1
    det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_2_SCALARS # VMC with drift
    2
    det_diamondC_2x1x1-Gaussian_pp_MSD_VMC_DMC_BATCH_MWALKERS-r2-t3_3_SCALARS # DMC
  )

else()
  message(VERBOSE "Skipping Periodic LCAO as not supported by CUDA build (QMC_CUDA=1)")
endif()
