if(NOT QMC_COMPLEX)
  #
  # LiH molecule pseudopotential tests
  #
  list(APPEND LIH_PP_SCALARS "kinetic" "0.6354491487 0.0013")
  list(APPEND LIH_PP_SCALARS "totenergy" "-0.7507264125 0.00043")
  list(APPEND LIH_PP_SCALARS "eeenergy" "0.4814819907 0.00061")
  list(APPEND LIH_PP_SCALARS "potential" "-1.3861755612 0.0013")
  list(APPEND LIH_PP_SCALARS "samples" "1600000 0.0")
  #  LIST(APPEND LIH_PP_SCALARS "flux" "0.0 0.0025")

  qmc_run_and_check(
    short-LiH_pp-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_pp"
    vmc_hf_noj
    vmc_hf_noj.in.xml
    1
    16
    TRUE
    0
    LIH_PP_SCALARS # VMC
  )

  #
  # LiH molecule pseudopotential tests with Jastrow from XML
  #
  list(APPEND LIH_PP_XML_SCALARS "kinetic" "0.6714 0.0022")
  list(APPEND LIH_PP_XML_SCALARS "totenergy" "-0.78492 0.00015")
  list(APPEND LIH_PP_XML_SCALARS "eeenergy" "0.42605 0.00087")
  list(APPEND LIH_PP_XML_SCALARS "potential" "-1.4563 0.0022")
  list(APPEND LIH_PP_XML_SCALARS "samples" "320000 0.0")
  #  LIST(APPEND LIH_PP_SCALARS "flux" "-0.0082 0.0051")

  qmc_run_and_check(
    short-LiH_pp-vmc_hf_sdj_xml
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_pp"
    vmc_hf_sdj_xml
    vmc_hf_sdj_xml.xml
    1
    16
    TRUE
    0
    LIH_PP_XML_SCALARS # VMC
  )

  qmc_run_and_check(
    short-LiH_pp-vmc_hf_sdj_xml
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_pp"
    vmc_hf_sdj_xml
    vmc_hf_sdj_xml.xml
    16
    1
    TRUE
    0
    LIH_PP_XML_SCALARS # VMC
  )

  qmc_run_and_check(
    short-LiH_pp-vmc_hf_sdj_xml
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_pp"
    vmc_hf_sdj_xml
    vmc_hf_sdj_xml.xml
    4
    4
    TRUE
    0
    LIH_PP_XML_SCALARS # VMC
  )

  #
  # LiH molecule pseudopotential tests with Jastrow from HDF5
  #

  qmc_run_and_check(
    short-LiH_pp-vmc_hf_sdj_hdf5
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_pp"
    vmc_hf_sdj_hdf5
    vmc_hf_sdj_hdf5.xml
    1
    16
    TRUE
    0
    LIH_PP_XML_SCALARS # VMC
  )

  qmc_run_and_check(
    short-LiH_pp-vmc_hf_sdj_hdf5
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_pp"
    vmc_hf_sdj_hdf5
    vmc_hf_sdj_hdf5.xml
    16
    1
    TRUE
    0
    LIH_PP_XML_SCALARS # VMC
  )

  qmc_run_and_check(
    short-LiH_pp-vmc_hf_sdj_hdf5
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_pp"
    vmc_hf_sdj_hdf5
    vmc_hf_sdj_hdf5.xml
    4
    4
    TRUE
    0
    LIH_PP_XML_SCALARS # VMC
  )

  #
  # LiH molecule pseudopotential tests with Jastrow from OLD HDF5
  #

  qmc_run_and_check(
    short-LiH_pp-vmc_hf_sdj_old_hdf5
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_pp"
    vmc_hf_sdj_old-hdf5
    vmc_hf_sdj_old-hdf5.xml
    1
    16
    TRUE
    0
    LIH_PP_XML_SCALARS # VMC
  )

  #
  # Deterministic tests
  #

  # LiH molecule pseudopotential tests
  if(QMC_MIXED_PRECISION)
    list(APPEND DET_LIH_PP_SCALARS "kinetic" "0.14526597 0.00009298")
    list(APPEND DET_LIH_PP_SCALARS "totenergy" "-0.72279678 0.00009313")
    list(APPEND DET_LIH_PP_SCALARS "eeenergy" "0.44908009 0.000001")
    list(APPEND DET_LIH_PP_SCALARS "potential" "-0.86806275 0.000001")
    list(APPEND DET_LIH_PP_SCALARS "ionion" "0.33179336 0.000001")
    list(APPEND DET_LIH_PP_SCALARS "localecp" "-1.64893620 0.00000134")
    list(APPEND DET_LIH_PP_SCALARS "nonlocalecp" "0.00000000 0.000001")
    list(APPEND DET_LIH_PP_SCALARS "samples" "9.00000000 0.0")
    list(APPEND DET_LIH_PP_SCALARS "flux" "1.05207942 0.00018638")
  else()
    list(APPEND DET_LIH_PP_SCALARS "kinetic" "0.14529655 0.000001")
    list(APPEND DET_LIH_PP_SCALARS "totenergy" "-0.72276621 0.000001")
    list(APPEND DET_LIH_PP_SCALARS "eeenergy" "0.44908026 0.000001")
    list(APPEND DET_LIH_PP_SCALARS "potential" "-0.86806275 0.000001")
    list(APPEND DET_LIH_PP_SCALARS "ionion" "0.33179335 0.000001")
    list(APPEND DET_LIH_PP_SCALARS "localecp" "-1.64893636 0.000001")
    list(APPEND DET_LIH_PP_SCALARS "nonlocalecp" "0.00000000 0.000001")
    list(APPEND DET_LIH_PP_SCALARS "samples" "9.00000000 0.0")
    list(APPEND DET_LIH_PP_SCALARS "flux" "1.05201758 0.000001")
  endif()

  qmc_run_and_check(
    deterministic-LiH_pp-vmc_hf_noj
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_pp"
    det_vmc_hf_noj
    det_vmc_hf_noj.in.xml
    1
    1
    TRUE
    0
    DET_LIH_PP_SCALARS # VMC
  )

  #  LiH molecule pseudopotential tests with Jastrow from XML
  if(QMC_MIXED_PRECISION)
    list(APPEND DET_LIH_PP_XML_SCALARS "kinetic" "0.82744207 0.00027385")
    list(APPEND DET_LIH_PP_XML_SCALARS "totenergy" "-0.78495167 0.00027387")
    list(APPEND DET_LIH_PP_XML_SCALARS "eeenergy" "0.40036895 0.000001")
    list(APPEND DET_LIH_PP_XML_SCALARS "potential" "-1.61249693 0.00000192")
    list(APPEND DET_LIH_PP_XML_SCALARS "ionion" "0.33179338 0.000001")
    list(APPEND DET_LIH_PP_XML_SCALARS "localecp" "-2.39561420 0.00000213")
    list(APPEND DET_LIH_PP_XML_SCALARS "nonlocalecp" "0.05095442 0.000001")
    list(APPEND DET_LIH_PP_XML_SCALARS "samples" "9.00000000 0.0")
    list(APPEND DET_LIH_PP_XML_SCALARS "flux" "-0.44733048 0.00054106")
  else()
    list(APPEND DET_LIH_PP_XML_SCALARS "kinetic" "0.82737098 0.000001")
    list(APPEND DET_LIH_PP_XML_SCALARS "totenergy" "-0.78512613 0.000001")
    list(APPEND DET_LIH_PP_XML_SCALARS "eeenergy" "0.40036884 0.000001")
    list(APPEND DET_LIH_PP_XML_SCALARS "potential" "-1.61249710 0.000001")
    list(APPEND DET_LIH_PP_XML_SCALARS "ionion" "0.33179337 0.000001")
    list(APPEND DET_LIH_PP_XML_SCALARS "localecp" "-2.39561374 0.000001")
    list(APPEND DET_LIH_PP_XML_SCALARS "nonlocalecp" "0.05095443 0.000001")
    list(APPEND DET_LIH_PP_XML_SCALARS "samples" "9.00000000 0.0")
    list(APPEND DET_LIH_PP_XML_SCALARS "flux" "-0.44718758 0.000001")
  endif()

  qmc_run_and_check(
    deterministic-LiH_pp-vmc_hf_sdj_xml
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_pp"
    det_vmc_hf_sdj_xml
    det_vmc_hf_sdj_xml.xml
    1
    1
    TRUE
    0
    DET_LIH_PP_XML_SCALARS # VMC
  )

  # LiH molecule pseudopotential tests with Jastrow from HDF5
  if(QMC_MIXED_PRECISION)
    list(APPEND DET_LIH_PP_HDF5_SCALARS "kinetic" "0.77796918 0.00031073")
    list(APPEND DET_LIH_PP_HDF5_SCALARS "totenergy" "-0.79725042 0.00031023")
    list(APPEND DET_LIH_PP_HDF5_SCALARS "eeenergy" "0.22789440 0.000001")
    list(APPEND DET_LIH_PP_HDF5_SCALARS "potential" "-1.57520943 0.0000022")
    list(APPEND DET_LIH_PP_HDF5_SCALARS "ionion" "0.33179338 0.000001")
    list(APPEND DET_LIH_PP_HDF5_SCALARS "localecp" "-2.13498816 0.00000215")
    list(APPEND DET_LIH_PP_HDF5_SCALARS "nonlocalecp" "0.00009136 0.000001")
    list(APPEND DET_LIH_PP_HDF5_SCALARS "samples" "9.00000000 0.0")
    list(APPEND DET_LIH_PP_HDF5_SCALARS "flux" "-0.31480187 0.0006186")
  else()
    list(APPEND DET_LIH_PP_HDF5_SCALARS "kinetic" "0.77797717 0.000001")
    list(APPEND DET_LIH_PP_HDF5_SCALARS "totenergy" "-0.79723250 0.000001")
    list(APPEND DET_LIH_PP_HDF5_SCALARS "eeenergy" "0.22789450 0.000001")
    list(APPEND DET_LIH_PP_HDF5_SCALARS "potential" "-1.57520957 0.000001")
    list(APPEND DET_LIH_PP_HDF5_SCALARS "ionion" "0.33179337 0.000001")
    list(APPEND DET_LIH_PP_HDF5_SCALARS "localecp" "-2.13498879 0.000001")
    list(APPEND DET_LIH_PP_HDF5_SCALARS "nonlocalecp" "0.00009136 0.000001")
    list(APPEND DET_LIH_PP_HDF5_SCALARS "samples" "9.00000000 0.0")
    list(APPEND DET_LIH_PP_HDF5_SCALARS "flux" "-0.31481816 0.000001")
  endif()

  qmc_run_and_check(
    deterministic-LiH_pp-vmc_hf_sdj_hdf5
    "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_pp"
    det_vmc_hf_sdj_hdf5
    det_vmc_hf_sdj_hdf5.xml
    1
    1
    TRUE
    0
    DET_LIH_PP_HDF5_SCALARS # VMC
  )

  # estimator tests
  include("${qmcpack_SOURCE_DIR}/CMake/python.cmake")
  check_python_reqs("numpy;h5py" diamond-estimator add_tests)
  if(add_tests)
    simple_run_and_check(
      short-LiH_pp-vmc_hf_sdj_estimator_spindensity
      "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_pp"
      vmc_hf_sdj_spindens_short.in.xml
      4
      4
      check_stats.py
      -s
      0
      -q
      spindensity
      -e
      2
      -c
      8
      -p
      vmc_hf_sdj_spindens_short
      -r
      qmc_ref/vmc_hf_sdj_spindens_short.s000.stat_ref_spindensity.dat)

    simple_run_and_check(
      long-LiH_pp-vmc_hf_sdj_estimator_spindensity
      "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_pp"
      vmc_hf_sdj_spindens_long.in.xml
      4
      4
      check_stats.py
      -s
      0
      -q
      spindensity
      -e
      20
      -c
      8
      -p
      vmc_hf_sdj_spindens_long
      -r
      qmc_ref/vmc_hf_sdj_spindens_long.s000.stat_ref_spindensity.dat)

    simple_run_and_check(
      short-LiH_pp-vmc_hf_sdj_estimator_energydensity_voronoi
      "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_pp"
      vmc_hf_sdj_edens_vor_short.in.xml
      4
      4
      check_stats.py
      -s
      0
      -q
      'energydensity,EDvoronoi'
      -e
      2
      -c
      2
      -p
      vmc_hf_sdj_edens_vor_short
      -r
      qmc_ref/vmc_hf_sdj_edens_vor_short.s000.stat_ref_energydensity.dat)

    simple_run_and_check(
      long-LiH_pp-vmc_hf_sdj_estimator_energydensity_voronoi
      "${qmcpack_SOURCE_DIR}/tests/molecules/LiH_pp"
      vmc_hf_sdj_edens_vor_long.in.xml
      4
      4
      check_stats.py
      -s
      0
      -q
      'energydensity,EDvoronoi'
      -e
      2
      -c
      2
      -p
      vmc_hf_sdj_edens_vor_long
      -r
      qmc_ref/vmc_hf_sdj_edens_vor_long.s000.stat_ref_energydensity.dat)
  endif()

else()
  message(VERBOSE "Skipping LiH_pp test because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
endif()
