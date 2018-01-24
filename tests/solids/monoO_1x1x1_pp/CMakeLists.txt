
  LIST(APPEND MONO_O_SCALARS "totenergy" "-31.776596 0.006")
  LIST(APPEND MONO_O_SCALARS "samples" "64000 0.0")
#  LIST(APPEND MONO_O_SCALARS "flux" "0.0 0.5")

  QMC_RUN_AND_CHECK(short-monoO_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 MONO_O_SCALARS # VMC
                    )

  LIST(APPEND MONO_O_J3_SCALARS "totenergy" "-31.913838 0.0025")
  LIST(APPEND MONO_O_J3_SCALARS "samples" "204800 0.0")
  QMC_RUN_AND_CHECK(short-monoO_1x1x1_pp-vmc_sdj3
                    "${CMAKE_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp"
                    qmc_j3_short
                    qmc_j3_short.in.xml
                    1 16
                    TRUE
                    0 MONO_O_J3_SCALARS # VMC
                    )

#
# Long test
#


  LIST(APPEND LONG_MONO_O_SCALARS "totenergy" "-31.776596 0.0002")
  LIST(APPEND LONG_MONO_O_SCALARS "samples" "61440000 0.0")
  LIST(APPEND LONG_MONO_O_SCALARS "flux" "0.0 0.05")

  QMC_RUN_AND_CHECK(long-monoO_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_MONO_O_SCALARS # VMC
                    )

  LIST(APPEND LONG_MONO_O_J3_SCALARS "totenergy" "-31.913838 0.0008")
  LIST(APPEND LONG_MONO_O_J3_SCALARS "samples" "2048000 0.0")
  LIST(APPEND LONG_MONO_O_J3_SCALARS "flux" "0.0 0.27")

  QMC_RUN_AND_CHECK(long-monoO_1x1x1_pp-vmc_sdj3
                    "${CMAKE_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp"
                    qmc_j3_long
                    qmc_j3_long.in.xml
                    1 16
                    TRUE
                    0 LONG_MONO_O_J3_SCALARS # VMC
                    )

