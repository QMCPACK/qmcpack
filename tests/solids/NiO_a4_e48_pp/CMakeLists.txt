if(NOT QMC_DATA)
  message_verbose("QMC_DATA not set. NiO_a4_e48_pp tests not added.")
elseif(NOT EXISTS ${QMC_DATA}/NiO)
  message("NiO directory under QMC_DATA does not exist. NiO_a4_e48_pp tests not added.")
else()

  set(H5_FILE NiO-fcc-supertwist111-supershift000-S1.h5)
  set(H5_FULL_PATH "${QMC_DATA}/NiO/${H5_FILE}")

  if(EXISTS ${H5_FULL_PATH})
    #symlink h5 file
    maybe_symlink(${H5_FULL_PATH} ${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp/${H5_FILE})

    #  LocalEnergy           =         -370.7944 +/-           0.0011
    #  Variance              =            12.836 +/-            0.052
    #  Kinetic               =           230.471 +/-            0.018
    #  LocalPotential        =          -601.265 +/-            0.018
    #  ElecElec              =           78.4679 +/-           0.0075
    #  LocalECP              =          -389.958 +/-            0.018
    #  NonLocalECP           =           -50.477 +/-            0.011
    #  IonIon                =           -239.30 +/-             0.00
    list(APPEND NIO_A4_E48_SCALARS "totenergy" "-370.7944 0.024")
    list(APPEND NIO_A4_E48_SCALARS "kinetic" "230.471 0.46")
    list(APPEND NIO_A4_E48_SCALARS "samples" "32000 0.0")

    qmc_run_and_check(
      short-NiO_a4_e48_pp-vmc_sdj
      "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
      NiO-fcc-S1-vmc
      NiO-vmc-short.in.xml
      1
      16
      TRUE
      1
      NIO_A4_E48_SCALARS # VMC
    )

    #  LocalEnergy           =  -371.320817 +/- 0.000430   13.153990 +/- 0.007166
    #  LocalEnergy           =  -371.288648 +/- 0.000530   13.352353 +/- 0.006475
    #  LocalEnergy           =  -371.283531 +/- 0.000519   13.590739 +/- 0.008039
    #  LocalEnergy           =  -371.283912 +/- 0.000493   13.588468 +/- 0.007902
    list(APPEND NIO_A4_E48_DMC_SCALARS "totenergy" "-371.320817 0.023")
    list(APPEND NIO_A4_E48_DMC_V0_SCALARS "totenergy" "-371.288648 0.018")
    list(APPEND NIO_A4_E48_DMC_V1_SCALARS "totenergy" "-371.283531 0.018")
    list(APPEND NIO_A4_E48_DMC_V3_SCALARS "totenergy" "-371.283912 0.018")

    qmc_run_and_check(
      short-NiO_a4_e48_pp-dmc_sdj
      "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
      NiO-fcc-S1-vmc-dmc
      NiO-vmc-dmc-short.in.xml
      1
      16
      TRUE
      2
      NIO_A4_E48_DMC_SCALARS # DMC no Tmove
      3
      NIO_A4_E48_DMC_V0_SCALARS # DMC Tmove v0
    )

    qmc_run_and_check(
      short-NiO_a4_e48_pp-delayed_update-dmc_sdj
      "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
      NiO-fcc-S1-vmc-dmc
      NiO-delayedupdate-vmc-dmc-short.in.xml
      1
      16
      TRUE
      2
      NIO_A4_E48_DMC_SCALARS # DMC no Tmove
      3
      NIO_A4_E48_DMC_V0_SCALARS # DMC Tmove v0
    )

    if(NOT QMC_CUDA)
      qmc_run_and_check(
        short-NiO_a4_e48_pp-dmc-TMv1v3_sdj
        "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
        NiO-fcc-S1-vmc-dmc
        NiO-vmc-dmc-TMv1v3-short.in.xml
        1
        16
        TRUE
        2
        NIO_A4_E48_DMC_V1_SCALARS # DMC Tmove v1
        3
        NIO_A4_E48_DMC_V3_SCALARS # DMC Tmove v3
      )

      qmc_run_and_check(
        short-NiO_a4_e48_pp-delayed_update-dmc-TMv1v3_sdj
        "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
        NiO-fcc-S1-vmc-dmc
        NiO-delayedupdate-vmc-dmc-TMv1v3-short.in.xml
        1
        16
        TRUE
        2
        NIO_A4_E48_DMC_V1_SCALARS # DMC Tmove v1
        3
        NIO_A4_E48_DMC_V3_SCALARS # DMC Tmove v3
      )
    else()
      message_verbose("Skipping NiO_a4_e48 Tmove v1/v3 tests on GPU")
    endif()

    #  LocalEnergy           =        -370.88498 +/-          0.00088
    #  Variance              =            11.118 +/-            0.031
    #  Kinetic               =           230.555 +/-            0.020
    #  LocalPotential        =          -601.440 +/-            0.020
    #  ElecElec              =           78.3624 +/-           0.0077
    #  LocalECP              =          -389.971 +/-            0.022
    #  NonLocalECP           =           -50.534 +/-            0.012
    #  IonIon                =           -239.30 +/-             0.00
    list(APPEND NIO_A4_E48_J3_SCALARS "totenergy" "-370.88498 0.023")
    list(APPEND NIO_A4_E48_J3_SCALARS "kinetic" "230.555 0.46")
    list(APPEND NIO_A4_E48_J3_SCALARS "samples" "32000 0.0")

    if(NOT QMC_CUDA)
      qmc_run_and_check(
        short-NiO_a4_e48_pp-vmc_sdj3
        "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
        NiO-fcc-S1-vmc
        NiO-j3-vmc-short.in.xml
        1
        16
        TRUE
        1
        NIO_A4_E48_J3_SCALARS # VMC
      )

      # Testing batched algorithm in NLPP.
      # Hitting all the routines using VirtualParticleSet.
      qmc_run_and_check(
        short-NiO_a4_e48-batched_pp-vmc_sdj3
        "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
        NiO-fcc-S1-vmc
        NiO-j3-batched-vmc-short.in.xml
        1
        16
        TRUE
        1
        NIO_A4_E48_J3_SCALARS # VMC
      )

    else()
      message_verbose("Skipping NiO_a4_e48 J3 tests for CUDA build (QMC_CUDA=1)")
    endif()

    if(NOT QMC_CUDA)
      list(APPEND NIO_A4_E48_L2_SCALARS "totenergy" "-371.119855 0.019547")
      list(APPEND NIO_A4_E48_L2_SCALARS "samples" "128000 0.0")

      qmc_run_and_check(
        short-NiO_a4_e48_pp-vmc_sdj_L2
        "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
        NiO-L2-short
        NiO-L2-short.in.xml
        16
        1
        TRUE
        1
        NIO_A4_E48_L2_SCALARS # VMC
      )

      qmc_run_and_check(
        short-NiO_a4_e48_pp-vmc_sdj_L2
        "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
        NiO-L2-short
        NiO-L2-short.in.xml
        4
        4
        TRUE
        1
        NIO_A4_E48_L2_SCALARS # VMC
      )
    else()
      message_verbose("Skipping NiO_a4_e48 L2 tests for CUDA build (QMC_CUDA=1)")
    endif()

    # Hybridrep is not implemented in legacy CUDA but should be error trapped. J3 also not implemented.    
    qmc_run_and_check(
      short-NiO_a4_e48-hybridrep-pp-vmc_sdj3
      "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
      NiO-fcc-S1-vmc
      NiO-hybridrep-j3-vmc-short.in.xml
      1
      16
      ${SUCCESS_STATUS_CUDA}
      1
      NIO_A4_E48_J3_SCALARS # VMC
    )

    # Testing batched algorithm in NLPP.
    # Hitting all the routines using VirtualParticleSet.
    qmc_run_and_check(
      short-NiO_a4_e48-hybridrep-batched_pp-vmc_sdj3
      "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
      NiO-fcc-S1-vmc
      NiO-hybridrep-j3-batched-vmc-short.in.xml
      1
      16
      ${SUCCESS_STATUS_CUDA}
      1
      NIO_A4_E48_J3_SCALARS # VMC
    )

    #
    # Long test
    #

    list(APPEND LONG_NIO_A4_E48_SCALARS "totenergy" "-370.7944 0.0070")
    list(APPEND LONG_NIO_A4_E48_SCALARS "samples" "320000 0.0")

    qmc_run_and_check(
      long-NiO_a4_e48_pp-vmc_sdj
      "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
      NiO-fcc-S1-vmc
      NiO-vmc-long.in.xml
      1
      16
      TRUE
      1
      LONG_NIO_A4_E48_SCALARS # VMC
    )

    list(APPEND LONG_NIO_A4_E48_J3_SCALARS "totenergy" "-370.88498 0.0066")
    list(APPEND LONG_NIO_A4_E48_J3_SCALARS "samples" "320000 0.0")

    qmc_run_and_check(
      long-NiO_a4_e48_pp-vmc_sdj3
      "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
      NiO-fcc-S1-vmc
      NiO-j3-vmc-long.in.xml
      1
      16
      TRUE
      1
      LONG_NIO_A4_E48_J3_SCALARS # VMC
    )

    if(NOT QMC_CUDA)
      list(APPEND LONG_NIO_A4_E48_L2_SCALARS "totenergy" "-371.119855 0.006451")
      list(APPEND LONG_NIO_A4_E48_L2_SCALARS "samples" "1280000 0.0")

      qmc_run_and_check(
        long-NiO_a4_e48_pp-vmc_sdj_L2
        "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
        NiO-L2-long
        NiO-L2-long.in.xml
        16
        1
        TRUE
        1
        LONG_NIO_A4_E48_L2_SCALARS # VMC
      )

      qmc_run_and_check(
        long-NiO_a4_e48_pp-vmc_sdj_L2
        "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
        NiO-L2-long
        NiO-L2-long.in.xml
        4
        4
        TRUE
        1
        LONG_NIO_A4_E48_L2_SCALARS # VMC
      )
    else()
      message_verbose("Skipping NiO_a4_e48 L2 tests on GPU")
    endif()

    #Determinsitic test

    if (QMC_CUDA)
     if (QMC_MIXED_PRECISION)
      if (QMC_COMPLEX)
       list(APPEND DET_NIO_A4_E48_SCALARS "totenergy" "-372.80830529 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "kinetic" "261.64612933 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "potential" "-634.45443462 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "eeenergy" "98.86604184 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "ionion" "-239.29814209 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "localecp" "-421.63781023 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "nonlocalecp" "-72.38452415 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "samples" "9 0.0")
      else()
       list(APPEND DET_NIO_A4_E48_SCALARS "totenergy" "-372.82709433 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "kinetic" "261.68320335 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "potential" "-634.51029768 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "eeenergy" "98.89647167 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "ionion" "-239.29814209 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "localecp" "-421.67137916 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "nonlocalecp" "-72.43724810 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "samples" "9 0.0")
      endif()
     else()
      if (QMC_COMPLEX)
       list(APPEND DET_NIO_A4_E48_SCALARS "totenergy" "-368.30237217 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "kinetic" "237.67396162 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "potential" "-605.97633379 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "eeenergy" "80.27766416 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "ionion" "-239.29814209 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "localecp" "-400.09667060 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "nonlocalecp" "-46.85918525 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "samples" "9 0.0")
      else()
       list(APPEND DET_NIO_A4_E48_SCALARS "totenergy" "-368.30237167 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "kinetic" "237.67180484 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "potential" "-605.97417652 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "eeenergy" "80.27750941 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "ionion" "-239.29814209 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "localecp" "-400.09547457 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "nonlocalecp" "-46.85806927 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "samples" "9 0.0")
      endif()
     endif()
    else()
     if (QMC_MIXED_PRECISION)
      if (QMC_COMPLEX)
       list(APPEND DET_NIO_A4_E48_SCALARS "totenergy" "-370.62466304 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "kinetic" "230.18358429 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "potential" "-600.80824733 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "eeenergy" "79.61912395 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "ionion" "-239.29823649 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "localecp" "-413.61387086 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "nonlocalecp" "-27.51526394 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "samples" "9 0.0")
      else()
       list(APPEND DET_NIO_A4_E48_SCALARS "totenergy" "-370.61130240 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "kinetic" "230.25580904 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "potential" "-600.86711144 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "eeenergy" "79.66720353 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "ionion" "-239.29823649 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "localecp" "-413.73287206 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "nonlocalecp" "-27.50320642 0.00005")
       list(APPEND DET_NIO_A4_E48_SCALARS "samples" "9 0.0")
      endif()
     else()
      if (QMC_COMPLEX)
       list(APPEND DET_NIO_A4_E48_SCALARS "totenergy" "-371.11389777 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "kinetic" "207.41466400 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "potential" "-578.52856177 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "eeenergy" "76.75959629 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "ionion" "-239.29814066 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "localecp" "-368.76230477 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "nonlocalecp" "-47.22771264 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "samples" "9 0.0")
      else()
       list(APPEND DET_NIO_A4_E48_SCALARS "totenergy" "-371.11474959 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "kinetic" "207.34064557 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "potential" "-578.45539516 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "eeenergy" "76.82246146 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "ionion" "-239.29814066 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "localecp" "-368.77476518 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "nonlocalecp" "-47.20495078 0.000001")
       list(APPEND DET_NIO_A4_E48_SCALARS "samples" "9 0.0")
      endif() 
     endif()    
    endif()

    qmc_run_and_check(
      deterministic-NiO_a4_e48_pp-vmc_sdj
      "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
      det_NiO-fcc-S1-vmc
      det_NiO-vmc.in.xml
      1
      1
      TRUE
      1
      DET_NIO_A4_E48_SCALARS # VMC
      )   

    if (NOT QMC_CUDA)
     if (QMC_MIXED_PRECISION)
      if (QMC_COMPLEX)
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "totenergy" "-371.10852051 0.00005")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "kinetic" "207.37844321 0.00005")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "potential" "-578.48696313 0.00005")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "eeenergy" "076.78619943 0.00005")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "ionion" "-239.29823649 0.00005")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "localecp" "-368.76812902 0.00005")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "nonlocalecp" "-47.20679705 0.00005")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "samples" "9 0.0")
      else()
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "totenergy" "-371.11052450 0.00005")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "kinetic" "207.47431382 0.00005")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "potential" "-578.58484291 0.00005")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "eeenergy" "76.70489285 0.00005")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "ionion" "-239.29823649 0.00005")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "localecp" "-368.75195719 0.00005")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "nonlocalecp" "-47.23954208 0.00005")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "samples" "9 0.0")
      endif()
     else()
      if (QMC_COMPLEX)
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "totenergy" "-371.11389777 0.000001")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "kinetic" "207.41466400 0.000001")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "potential" "-578.52856177 0.000001")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "eeenergy" "76.75959629 0.000001")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "ionion" "-239.29814066 0.000001")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "localecp" "-368.76230477 0.000001")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "nonlocalecp" "-47.22771264 0.000001")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "samples" "9 0.0")
      else()
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "totenergy" "-371.11474959 0.000001")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "kinetic" "207.34064557 0.000001")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "potential" "-578.45539516 0.000001")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "eeenergy" "76.82246146 0.000001")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "ionion" "-239.29814066 0.000001")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "localecp" "-368.77476518 0.000001")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "nonlocalecp" "-47.20495078 0.000001")
       list(APPEND DET_NIO_BATCHED_A4_E48_SCALARS "samples" "9 0.0")
      endif()
     endif()

     qmc_run_and_check(
       deterministic-NiO_a4_e48_pp-vmc_sdj_batched
       "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
       det_NiO-batched-fcc-S1-vmc 
       det_NiO-batched-vmc.in.xml
       1
       1
       TRUE
       1
       DET_NIO_BATCHED_A4_E48_SCALARS # VMC
       )
     endif()

  else()
    message_verbose("NiO_a4_e48_pp not added because the corresponding h5 file not found: ${H5_FULL_PATH}")
  endif()

endif()
