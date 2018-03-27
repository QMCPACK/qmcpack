IF (NOT QMC_DATA)
  MESSAGE("QMC_DATA not set. NiO_a4_e48_pp tests not added.")
ELSEIF(NOT EXISTS ${QMC_DATA}/NiO)
  MESSAGE("NiO directory under QMC_DATA does not exist. NiO_a4_e48_pp tests not added.")
ELSE()

SET(H5_FILE NiO-fcc-supertwist111-supershift000-S1.h5)
SET(H5_FULL_PATH "${QMC_DATA}/NiO/${H5_FILE}")

IF (EXISTS ${H5_FULL_PATH})
  #symlink h5 file
  MAYBE_SYMLINK(${H5_FULL_PATH} ${CMAKE_CURRENT_BINARY_DIR}/${H5_FILE})
  #symlink pp files
  MAYBE_SYMLINK(${CMAKE_SOURCE_DIR}/tests/performance/NiO/O.xml ${CMAKE_CURRENT_BINARY_DIR}/O.xml)
  MAYBE_SYMLINK(${CMAKE_SOURCE_DIR}/tests/performance/NiO/Ni.opt.xml ${CMAKE_CURRENT_BINARY_DIR}/Ni.opt.xml)

  LIST(APPEND NIO_A4_E48_SCALARS "totenergy" "-370.785413 0.03")
  LIST(APPEND NIO_A4_E48_SCALARS "samples" "32000 0.0")
#  LIST(APPEND NIO_A4_E48_SCALARS "flux" "0.0 0.5")

  QMC_RUN_AND_CHECK(short-NiO_a4_e48_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-fcc-S1-vmc
                    NiO-vmc-short.in.xml
                    1 16
                    TRUE
                    1 NIO_A4_E48_SCALARS # VMC
                    )

  LIST(APPEND NIO_A4_E48_J3_SCALARS "totenergy" "-370.886859 0.021")
  LIST(APPEND NIO_A4_E48_J3_SCALARS "samples" "32000 0.0")
  QMC_RUN_AND_CHECK(short-NiO_a4_e48_pp-vmc_sdj3
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-fcc-S1-vmc
                    NiO-j3-vmc-short.in.xml
                    1 16
                    TRUE
                    1 NIO_A4_E48_J3_SCALARS # VMC
                    )

#
# Long test
#


  LIST(APPEND LONG_NIO_A4_E48_SCALARS "totenergy" "-31.776596 0.0002")
  LIST(APPEND LONG_NIO_A4_E48_SCALARS "samples" "61440000 0.0")
  LIST(APPEND LONG_NIO_A4_E48_SCALARS "flux" "0.0 0.05")

  QMC_RUN_AND_CHECK(long-NiO_a4_e48_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-fcc-S1-vmc
                    NiO-vmc-long.in.xml
                    1 16
                    TRUE
                    1 LONG_NIO_A4_E48_SCALARS # VMC
                    )

  LIST(APPEND LONG_NIO_A4_E48_J3_SCALARS "totenergy" "-31.913838 0.0008")
  LIST(APPEND LONG_NIO_A4_E48_J3_SCALARS "samples" "2048000 0.0")
  LIST(APPEND LONG_NIO_A4_E48_J3_SCALARS "flux" "0.0 0.27")

  QMC_RUN_AND_CHECK(long-NiO_a4_e48_pp-vmc_sdj3
                    "${CMAKE_SOURCE_DIR}/tests/solids/NiO_a4_e48_pp"
                    NiO-fcc-S1-vmc
                    NiO-j3-vmc-long.in.xml
                    1 16
                    TRUE
                    1 LONG_NIO_A4_E48_J3_SCALARS # VMC
                    )

ELSE()
  MESSAGE("NiO_a4_e48_pp not added because the corresponding h5 file not found: ${H5_FULL_PATH}")
ENDIF()

ENDIF()
