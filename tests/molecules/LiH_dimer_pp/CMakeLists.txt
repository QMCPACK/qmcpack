
IF (NOT QMC_CUDA)
 IF (NOT QMC_COMPLEX)
#
# LiH molecular dimer pseudopotential tests
#
  LIST(APPEND LIH_PP_SCALARS "kinetic" "0.6354491487 0.0013")
  LIST(APPEND LIH_PP_SCALARS "totenergy" "-0.7507264125 0.00043")
  LIST(APPEND LIH_PP_SCALARS "eeenergy" "0.4814819907 0.00061")
  LIST(APPEND LIH_PP_SCALARS "potential" "-1.3861755612 0.0013")
  LIST(APPEND LIH_PP_SCALARS "samples" "1600000 0.0")
#  LIST(APPEND LIH_PP_SCALARS "flux" "0.0 0.0025")

  QMC_RUN_AND_CHECK(short-LiH_dimer_pp-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp"
                    vmc_hf_noj
                    vmc_hf_noj.in.xml
                    1 16
                    TRUE
                    0 LIH_PP_SCALARS # VMC
                    )

#
# LiH molecular dimer pseudopotential tests with Jastrow from XML
#
  LIST(APPEND LIH_PP_XML_SCALARS "kinetic" "0.6714 0.0022")
  LIST(APPEND LIH_PP_XML_SCALARS "totenergy" "-0.78492 0.00015")
  LIST(APPEND LIH_PP_XML_SCALARS "eeenergy" "0.42605 0.00087")
  LIST(APPEND LIH_PP_XML_SCALARS "potential" "-1.4563 0.0022")
  LIST(APPEND LIH_PP_XML_SCALARS "samples" "320000 0.0")
#  LIST(APPEND LIH_PP_SCALARS "flux" "-0.0082 0.0051")

  QMC_RUN_AND_CHECK(short-LiH_dimer_pp-vmc_hf_sdj_xml
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp"
                    vmc_hf_sdj_xml
                    vmc_hf_sdj_xml.xml
                    1 16
                    TRUE
                    0 LIH_PP_XML_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-LiH_dimer_pp-vmc_hf_sdj_xml
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp"
                    vmc_hf_sdj_xml
                    vmc_hf_sdj_xml.xml
                    16 1
                    TRUE
                    0 LIH_PP_XML_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-LiH_dimer_pp-vmc_hf_sdj_xml
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp"
                    vmc_hf_sdj_xml
                    vmc_hf_sdj_xml.xml
                    4 4
                    TRUE
                    0 LIH_PP_XML_SCALARS # VMC
                    )


#
# LiH molecular dimer pseudopotential tests with Jastrow from HDF5
#

  QMC_RUN_AND_CHECK(short-LiH_dimer_pp-vmc_hf_sdj_hdf5
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp"
                    vmc_hf_sdj_hdf5
                    vmc_hf_sdj_hdf5.xml
                    1 16
                    TRUE
                    0 LIH_PP_XML_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-LiH_dimer_pp-vmc_hf_sdj_hdf5
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp"
                    vmc_hf_sdj_hdf5
                    vmc_hf_sdj_hdf5.xml
                    16 1
                    TRUE
                    0 LIH_PP_XML_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(short-LiH_dimer_pp-vmc_hf_sdj_hdf5
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp"
                    vmc_hf_sdj_hdf5
                    vmc_hf_sdj_hdf5.xml
                    4 4
                    TRUE
                    0 LIH_PP_XML_SCALARS # VMC
                    )


#
# Deterministic tests
#

# LiH molecular dimer pseudopotential tests
IF(QMC_MIXED_PRECISION)
  LIST(APPEND DET_LIH_PP_SCALARS "kinetic" "0.14529655 0.00012115")
  LIST(APPEND DET_LIH_PP_SCALARS "totenergy" "-0.72276621 0.00012158")
  LIST(APPEND DET_LIH_PP_SCALARS "eeenergy" "0.44908026 0.00000069")
  LIST(APPEND DET_LIH_PP_SCALARS "potential" "-0.86806275 0.00000076")
  LIST(APPEND DET_LIH_PP_SCALARS "ionion" "0.33179335 0.0000002")
  LIST(APPEND DET_LIH_PP_SCALARS "localecp" "-1.64893636 0.0000014")
  LIST(APPEND DET_LIH_PP_SCALARS "nonlocalecp" "0.00000000 0.000001")
  LIST(APPEND DET_LIH_PP_SCALARS "samples" "9.00000000 0.0")
  LIST(APPEND DET_LIH_PP_SCALARS "flux" "1.05201758 0.00024441")
ELSE()
  LIST(APPEND DET_LIH_PP_SCALARS "kinetic" "0.14529655 0.000001")
  LIST(APPEND DET_LIH_PP_SCALARS "totenergy" "-0.72276621 0.000001")
  LIST(APPEND DET_LIH_PP_SCALARS "eeenergy" "0.44908026 0.000001")
  LIST(APPEND DET_LIH_PP_SCALARS "potential" "-0.86806275 0.000001")
  LIST(APPEND DET_LIH_PP_SCALARS "ionion" "0.33179335 0.000001")
  LIST(APPEND DET_LIH_PP_SCALARS "localecp" "-1.64893636 0.000001")
  LIST(APPEND DET_LIH_PP_SCALARS "nonlocalecp" "0.00000000 0.000001")
  LIST(APPEND DET_LIH_PP_SCALARS "samples" "9.00000000 0.0")
  LIST(APPEND DET_LIH_PP_SCALARS "flux" "1.05201758 0.000001")
ENDIF()

  QMC_RUN_AND_CHECK(deterministic-LiH_dimer_pp-vmc_hf_noj
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp"
                    det_vmc_hf_noj
                    det_vmc_hf_noj.in.xml
                    1 1
                    TRUE
                    0 DET_LIH_PP_SCALARS # VMC
                    )

#  LiH molecular dimer pseudopotential tests with Jastrow from XML
IF(QMC_MIXED_PRECISION)
  LIST(APPEND DET_LIH_PP_XML_SCALARS "kinetic" "0.82737098 0.00034316")
  LIST(APPEND DET_LIH_PP_XML_SCALARS "totenergy" "-0.78502242 0.00034224")
  LIST(APPEND DET_LIH_PP_XML_SCALARS "eeenergy" "0.40036884 0.00000032")
  LIST(APPEND DET_LIH_PP_XML_SCALARS "potential" "-1.61239340 0.00000155")
  LIST(APPEND DET_LIH_PP_XML_SCALARS "ionion" "0.33179337 0.0000002")
  LIST(APPEND DET_LIH_PP_XML_SCALARS "localecp" "-2.39561374 0.00000189")
  LIST(APPEND DET_LIH_PP_XML_SCALARS "nonlocalecp" "0.05105813 0.0000002")
  LIST(APPEND DET_LIH_PP_XML_SCALARS "samples" "9.00000000 0.0")
  LIST(APPEND DET_LIH_PP_XML_SCALARS "flux" "-0.44718758 0.00068291")
ELSE()
  LIST(APPEND DET_LIH_PP_XML_SCALARS "kinetic" "0.82737098 0.000001")
  LIST(APPEND DET_LIH_PP_XML_SCALARS "totenergy" "-0.78502242 0.000001")
  LIST(APPEND DET_LIH_PP_XML_SCALARS "eeenergy" "0.40036884 0.000001")
  LIST(APPEND DET_LIH_PP_XML_SCALARS "potential" "-1.61239340 0.000001")
  LIST(APPEND DET_LIH_PP_XML_SCALARS "ionion" "0.33179337 0.000001")
  LIST(APPEND DET_LIH_PP_XML_SCALARS "localecp" "-2.39561374 0.000001")
  LIST(APPEND DET_LIH_PP_XML_SCALARS "nonlocalecp" "0.05105813 0.000001")
  LIST(APPEND DET_LIH_PP_XML_SCALARS "samples" "9.00000000 0.0")
  LIST(APPEND DET_LIH_PP_XML_SCALARS "flux" "-0.44718758 0.000001")
ENDIF()

  QMC_RUN_AND_CHECK(deterministic-LiH_dimer_pp-vmc_hf_sdj_xml
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp"
                    det_vmc_hf_sdj_xml
                    det_vmc_hf_sdj_xml.xml
                    1 1
                    TRUE
                    0 DET_LIH_PP_XML_SCALARS # VMC
                    ) 

# LiH molecular dimer pseudopotential tests with Jastrow from HDF5
IF(QMC_MIXED_PRECISION)
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "kinetic" "0.77797717 0.00025129")
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "totenergy" "-0.79724298 0.00025012")
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "eeenergy" "0.22789450 0.0000002")
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "potential" "-1.57522014 0.00000236")
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "ionion" "0.33179337 0.0000002")
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "localecp" "-2.13498879 0.00000248")
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "nonlocalecp" "0.00008078 0.000001")
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "samples" "9.00000000 0.0")
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "flux" "-0.31481816 0.00050133")
ELSE()
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "kinetic" "0.77797717 0.000001")
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "totenergy" "-0.79724298 0.000001")
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "eeenergy" "0.22789450 0.000001")
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "potential" "-1.57522014 0.000001")
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "ionion" "0.33179337 0.000001")
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "localecp" "-2.13498879 0.000001")
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "nonlocalecp" "0.00008078 0.000001")
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "samples" "9.00000000 0.0")
  LIST(APPEND DET_LIH_PP_HDF5_SCALARS "flux" "-0.31481816 0.000001")
ENDIF()

  QMC_RUN_AND_CHECK(deterministic-LiH_dimer_pp-vmc_hf_sdj_hdf5
                    "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp"
                    det_vmc_hf_sdj_hdf5
                    det_vmc_hf_sdj_hdf5.xml
                    1 1
                    TRUE
                    0 DET_LIH_PP_HDF5_SCALARS # VMC
                    )


  # estimator tests
  INCLUDE("${qmcpack_SOURCE_DIR}/CMake/python.cmake")
  set(stat_python_reqs numpy;h5py)
  CHECK_PYTHON_REQS(stat_python_reqs diamond-estimator add_tests)
  if (add_tests)
    SIMPLE_RUN_AND_CHECK(
      short-LiH_dimer_pp-vmc_hf_sdj_estimator_spindensity
      "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp"
      vmc_hf_sdj_spindens_short.in.xml
      4 4
      check_stats.py -s 0 -q spindensity -e 2 -c 8 -p vmc_hf_sdj_spindens_short -r qmc_ref/vmc_hf_sdj_spindens_short.s000.stat_ref_spindensity.dat
      )

    SIMPLE_RUN_AND_CHECK(
      long-LiH_dimer_pp-vmc_hf_sdj_estimator_spindensity
      "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp"
      vmc_hf_sdj_spindens_long.in.xml
      4 4
      check_stats.py -s 0 -q spindensity -e 20 -c 8 -p vmc_hf_sdj_spindens_long -r qmc_ref/vmc_hf_sdj_spindens_long.s000.stat_ref_spindensity.dat
      )

    SIMPLE_RUN_AND_CHECK(
      short-LiH_dimer_pp-vmc_hf_sdj_estimator_energydensity_voronoi
      "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp"
      vmc_hf_sdj_edens_vor_short.in.xml
      4 4
      check_stats.py -s 0 -q 'energydensity,EDvoronoi' -e 2 -c 2 -p vmc_hf_sdj_edens_vor_short -r qmc_ref/vmc_hf_sdj_edens_vor_short.s000.stat_ref_energydensity.dat
      )

    SIMPLE_RUN_AND_CHECK(
      long-LiH_dimer_pp-vmc_hf_sdj_estimator_energydensity_voronoi
      "${CMAKE_SOURCE_DIR}/tests/molecules/LiH_dimer_pp"
      vmc_hf_sdj_edens_vor_long.in.xml
      4 4
      check_stats.py -s 0 -q 'energydensity,EDvoronoi' -e 2 -c 2 -p vmc_hf_sdj_edens_vor_long -r qmc_ref/vmc_hf_sdj_edens_vor_long.s000.stat_ref_energydensity.dat
      )
  endif()

 ELSE()
    MESSAGE_VERBOSE("Skipping LiH_dimer_pp test because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE_VERBOSE("Skipping LiH_dimer_pp test because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

