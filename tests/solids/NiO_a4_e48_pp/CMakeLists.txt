IF (NOT QMC_DATA)
  MESSAGE("QMC_DATA not set. NiO_a4_e48_pp tests not added.")
ELSEIF(NOT EXISTS ${QMC_DATA}/NiO)
  MESSAGE("NiO directory under QMC_DATA does not exist. NiO_a4_e48_pp tests not added.")
ELSE()

SET(H5_FILE NiO-fcc-supertwist111-supershift000-S1.h5)
SET(H5_FULL_PATH "${QMC_DATA}/NiO/${H5_FILE}")

IF (EXISTS ${H5_FULL_PATH})
  #symlink h5 file
  MAYBE_SYMLINK(${H5_FULL_PATH} ${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp/${H5_FILE})


#  LocalEnergy           =         -370.7944 +/-           0.0011
#  Variance              =            12.836 +/-            0.052
#  Kinetic               =           230.471 +/-            0.018
#  LocalPotential        =          -601.265 +/-            0.018
#  ElecElec              =           78.4679 +/-           0.0075
#  LocalECP              =          -389.958 +/-            0.018
#  NonLocalECP           =           -50.477 +/-            0.011
#  IonIon                =           -239.30 +/-             0.00
  LIST(APPEND NIO_A4_E48_SCALARS "totenergy" "-370.7944 0.024")
  LIST(APPEND NIO_A4_E48_SCALARS "kinetic" "230.471 0.46")
  LIST(APPEND NIO_A4_E48_SCALARS "samples" "32000 0.0")

  QMC_RUN_AND_CHECK(short-NiO_a4_e48_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-fcc-S1-vmc
                    NiO-vmc-short.in.xml
                    1 16
                    TRUE
                    1 NIO_A4_E48_SCALARS # VMC
                    )

#  LocalEnergy           =  -371.322168 +/- 0.000843   13.223478 +/- 0.062204
#  LocalEnergy           =  -371.288408 +/- 0.000773   13.344812 +/- 0.011638
#  LocalEnergy           =  -371.283854 +/- 0.001185   13.587116 +/- 0.019469
#  LocalEnergy           =  -371.284416 +/- 0.000750   13.613860 +/- 0.013540
  LIST(APPEND NIO_A4_E48_DMC_SCALARS "totenergy" "-371.322168 0.023")
  LIST(APPEND NIO_A4_E48_DMC_V0_SCALARS "totenergy" "-371.288408 0.018")
  LIST(APPEND NIO_A4_E48_DMC_V1_SCALARS "totenergy" "-371.283854 0.018")
  LIST(APPEND NIO_A4_E48_DMC_V3_SCALARS "totenergy" "-371.284416 0.018")

  QMC_RUN_AND_CHECK(short-NiO_a4_e48_pp-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-fcc-S1-vmc-dmc
                    NiO-vmc-dmc-short.in.xml
                    1 16
                    TRUE
                    2 NIO_A4_E48_DMC_SCALARS # DMC no Tmove
                    3 NIO_A4_E48_DMC_V0_SCALARS # DMC no Tmove
                    4 NIO_A4_E48_DMC_V1_SCALARS # DMC no Tmove
                    5 NIO_A4_E48_DMC_V3_SCALARS # DMC no Tmove
                    )

  QMC_RUN_AND_CHECK(short-NiO_a4_e48_pp-delayed_update-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-fcc-S1-vmc-dmc
                    NiO-delayedupdate-vmc-dmc-short.in.xml
                    1 16
                    TRUE
                    2 NIO_A4_E48_DMC_SCALARS # DMC no Tmove
                    3 NIO_A4_E48_DMC_V0_SCALARS # DMC no Tmove
                    4 NIO_A4_E48_DMC_V1_SCALARS # DMC no Tmove
                    5 NIO_A4_E48_DMC_V3_SCALARS # DMC no Tmove
                    )

#  LocalEnergy           =        -370.88498 +/-          0.00088
#  Variance              =            11.118 +/-            0.031
#  Kinetic               =           230.555 +/-            0.020
#  LocalPotential        =          -601.440 +/-            0.020
#  ElecElec              =           78.3624 +/-           0.0077
#  LocalECP              =          -389.971 +/-            0.022
#  NonLocalECP           =           -50.534 +/-            0.012
#  IonIon                =           -239.30 +/-             0.00
  LIST(APPEND NIO_A4_E48_J3_SCALARS "totenergy" "-370.88498 0.023")
  LIST(APPEND NIO_A4_E48_J3_SCALARS "kinetic" "230.555 0.46")
  LIST(APPEND NIO_A4_E48_J3_SCALARS "samples" "32000 0.0")

  QMC_RUN_AND_CHECK(short-NiO_a4_e48_pp-vmc_sdj3
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-fcc-S1-vmc
                    NiO-j3-vmc-short.in.xml
                    1 16
                    TRUE
                    1 NIO_A4_E48_J3_SCALARS # VMC
                    )

  # Testing batched algorithm in NLPP.
  # Hitting all the routines using VirtualParticleSet.
  QMC_RUN_AND_CHECK(short-NiO_a4_e48-batched_pp-vmc_sdj3
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-fcc-S1-vmc
                    NiO-j3-batched-vmc-short.in.xml
                    1 16
                    TRUE
                    1 NIO_A4_E48_J3_SCALARS # VMC
                    )



IF(NOT QMC_CUDA)
  LIST(APPEND NIO_A4_E48_L2_SCALARS "totenergy" "-371.119855 0.019547")
  LIST(APPEND NIO_A4_E48_L2_SCALARS "samples" "128000 0.0")

  QMC_RUN_AND_CHECK(short-NiO_a4_e48_pp-vmc_sdj_L2
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-L2-short
                    NiO-L2-short.in.xml
                    16 1
                    TRUE
                    1 NIO_A4_E48_L2_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-NiO_a4_e48_pp-vmc_sdj_L2
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-L2-short
                    NiO-L2-short.in.xml
                    4 4
                    TRUE
                    1 NIO_A4_E48_L2_SCALARS # VMC
                    )
ELSE()
  MESSAGE("Skipping NiO_a4_e48-l2 tests on gpu")
ENDIF()


IF(ENABLE_SOA OR QMC_CUDA)
  QMC_RUN_AND_CHECK(short-NiO_a4_e48-hybridrep-pp-vmc_sdj3
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-fcc-S1-vmc
                    NiO-hybridrep-j3-vmc-short.in.xml
                    1 16
                    TRUE
                    1 NIO_A4_E48_J3_SCALARS # VMC
                    )

  # Testing batched algorithm in NLPP.
  # Hitting all the routines using VirtualParticleSet.
  QMC_RUN_AND_CHECK(short-NiO_a4_e48-hybridrep-batched_pp-vmc_sdj3
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-fcc-S1-vmc
                    NiO-hybridrep-j3-batched-vmc-short.in.xml
                    1 16
                    TRUE
                    1 NIO_A4_E48_J3_SCALARS # VMC
                    )


ELSE()
  MESSAGE("Skipping NiO_a4_e48-hybridrep tests because neither SOA nor CUDA is enabled")
ENDIF()
#
# Long test
#


  LIST(APPEND LONG_NIO_A4_E48_SCALARS "totenergy" "-370.7944 0.0070")
  LIST(APPEND LONG_NIO_A4_E48_SCALARS "samples" "320000 0.0")

  QMC_RUN_AND_CHECK(long-NiO_a4_e48_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-fcc-S1-vmc
                    NiO-vmc-long.in.xml
                    1 16
                    TRUE
                    1 LONG_NIO_A4_E48_SCALARS # VMC
                    )

  LIST(APPEND LONG_NIO_A4_E48_J3_SCALARS "totenergy" "-370.88498 0.0066")
  LIST(APPEND LONG_NIO_A4_E48_J3_SCALARS "samples" "320000 0.0")

  QMC_RUN_AND_CHECK(long-NiO_a4_e48_pp-vmc_sdj3
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-fcc-S1-vmc
                    NiO-j3-vmc-long.in.xml
                    1 16
                    TRUE
                    1 LONG_NIO_A4_E48_J3_SCALARS # VMC
                    )


IF(NOT QMC_CUDA)
  LIST(APPEND LONG_NIO_A4_E48_L2_SCALARS "totenergy" "-371.119855 0.006451")
  LIST(APPEND LONG_NIO_A4_E48_L2_SCALARS "samples" "1280000 0.0")

  QMC_RUN_AND_CHECK(long-NiO_a4_e48_pp-vmc_sdj_L2
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-L2-long
                    NiO-L2-long.in.xml
                    16 1
                    TRUE
                    1 LONG_NIO_A4_E48_L2_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(long-NiO_a4_e48_pp-vmc_sdj_L2
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-L2-long
                    NiO-L2-long.in.xml
                    4 4
                    TRUE
                    1 LONG_NIO_A4_E48_L2_SCALARS # VMC
                    )
ELSE()
  MESSAGE("Skipping NiO_a4_e48-l2 tests on gpu")
ENDIF()


ELSE()
  MESSAGE("NiO_a4_e48_pp not added because the corresponding h5 file not found: ${H5_FULL_PATH}")
ENDIF()

ENDIF()
