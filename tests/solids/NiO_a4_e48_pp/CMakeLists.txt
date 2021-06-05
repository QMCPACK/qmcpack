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
      message("Skipping NiO_a4_e48 Tmove v1/v3 tests on GPU")
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
      message("Skipping NiO_a4_e48 J3 tests on GPU")
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
      message("Skipping NiO_a4_e48 L2 tests on GPU")
    endif()

    qmc_run_and_check(
      short-NiO_a4_e48-hybridrep-pp-vmc_sdj3
      "${qmcpack_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
      NiO-fcc-S1-vmc
      NiO-hybridrep-j3-vmc-short.in.xml
      1
      16
      TRUE
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
      TRUE
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
      message("Skipping NiO_a4_e48 L2 tests on GPU")
    endif()

  else()
    message("NiO_a4_e48_pp not added because the corresponding h5 file not found: ${H5_FULL_PATH}")
  endif()

endif()
