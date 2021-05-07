IF (NOT QMC_DATA)
  MESSAGE("QMC_DATA not set. NiO_afqmc tests not added.")
ELSEIF(NOT EXISTS ${QMC_DATA}/NiO)
  MESSAGE("NiO directory under QMC_DATA does not exist. NiO_afqmc tests not added.")
ELSE()

SET(H5_FILE NiO_afm_fcidump.h5)
SET(H5_FULL_PATH "${QMC_DATA}/NiO/${H5_FILE}")

IF (EXISTS ${H5_FULL_PATH})
  #symlink h5 file
  MAYBE_SYMLINK(${H5_FULL_PATH} ${CMAKE_CURRENT_BINARY_DIR}/${H5_FILE})
  #symlink wfn file
  MAYBE_SYMLINK(${QMC_DATA}/NiO/NiO_afm_wfn.dat ${CMAKE_CURRENT_BINARY_DIR}/NiO_afm_wfn.dat)

  LIST(APPEND NIO_AFM_SCALARS "Eloc" "-379.2037 0.0282")

  QMC_RUN_AND_CHECK(long-NiO_afm-afqmc
                    "${qmcpack_SOURCE_DIR}/tests/solids/NiO_afqmc"
                    qmc_long_afm
                    qmc_long_afm.in.xml
                    16 1
                    TRUE
                    0 NIO_AFM_SCALARS # AFQMC
                    )


  MAYBE_SYMLINK(${QMC_DATA}/NiO/NiO_nm_choldump.h5 ${CMAKE_CURRENT_BINARY_DIR}/NiO_nm_choldump.h5)


  LIST(APPEND NIO_NM_SCALARS "Eloc" "-379.2024 0.0192")

  QMC_RUN_AND_CHECK(long-NiO_nm-afqmc
                    "${qmcpack_SOURCE_DIR}/tests/solids/NiO_afqmc"
                    qmc_long_nm
                    qmc_long_nm.in.xml
                    16 1
                    TRUE
                    0 NIO_NM_SCALARS # AFQMC
                    )


ELSE()
  MESSAGE("NiO_afqmc not added because the corresponding h5 file not found: ${H5_FULL_PATH}")
ENDIF()

ENDIF()
