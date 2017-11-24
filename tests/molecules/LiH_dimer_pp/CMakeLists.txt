
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
 ELSE()
    MESSAGE("Skipping LiH_dimer_pp test because gaussian basis sets are not supported by complex build (QMC_COMPLEX=1)")
 ENDIF()
ELSE()
    MESSAGE("Skipping LiH_dimer_pp test because gaussian basis sets are not supported by CUDA build (QMC_CUDA=1)")
ENDIF()

