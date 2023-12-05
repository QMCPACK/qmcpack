if(NOT QMC_DATA)
  message("QMC_DATA not set. NiO_afqmc tests not added.")
elseif(NOT EXISTS ${QMC_DATA}/NiO)
  message("NiO directory under QMC_DATA does not exist. NiO_afqmc tests not added.")
else()

  set(H5_FILE NiO_afm_fcidump.h5)
  set(H5_FULL_PATH "${QMC_DATA}/NiO/${H5_FILE}")

  if(EXISTS ${H5_FULL_PATH})
    #symlink h5 file
    maybe_symlink(${H5_FULL_PATH} ${CMAKE_CURRENT_BINARY_DIR}/${H5_FILE})
    #symlink wfn file
    maybe_symlink(${QMC_DATA}/NiO/NiO_afm_wfn.dat ${CMAKE_CURRENT_BINARY_DIR}/NiO_afm_wfn.dat)

    list(APPEND NIO_AFM_SCALARS "Eloc" "-379.2037 0.0282")

    qmc_run_and_check(
      long-NiO_afm-afqmc
      "${qmcpack_SOURCE_DIR}/tests/solids/NiO_afqmc"
      qmc_long_afm
      qmc_long_afm.in.xml
      16
      1
      TRUE
      0
      NIO_AFM_SCALARS # AFQMC
    )

    maybe_symlink(${QMC_DATA}/NiO/NiO_nm_choldump.h5 ${CMAKE_CURRENT_BINARY_DIR}/NiO_nm_choldump.h5)

    list(APPEND NIO_NM_SCALARS "Eloc" "-379.2024 0.0192")

    qmc_run_and_check(
      long-NiO_nm-afqmc
      "${qmcpack_SOURCE_DIR}/tests/solids/NiO_afqmc"
      qmc_long_nm
      qmc_long_nm.in.xml
      16
      1
      TRUE
      0
      NIO_NM_SCALARS # AFQMC
    )

  else()
    message("NiO_afqmc not added because the corresponding h5 file not found: ${H5_FULL_PATH}")
  endif()

endif()
