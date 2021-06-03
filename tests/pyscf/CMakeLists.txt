#############################################################
# Add tests to ctest
#############################################################

  INCLUDE( "${qmcpack_SOURCE_DIR}/CMake/run_pyscf.cmake" )
  SET(LAST_TEST_NAME "NONE")
  SET(PYSCF_TEST_NAME "NONE")

  RUN_PYSCF_TEST(pyscf-diamond_1x1x1_pp_LCAO ${qmcpack_SOURCE_DIR}/examples/solids/pyscf-inputs
                 diamondC_1x1x1_pp_LCAO PYSCF_TEST_NAME )

  SET_TESTS_PROPERTIES( ${PYSCF_TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "successfully saved to QMCPACK HDF5" )

  SOFTLINK_H5( ${PYSCF_TEST_NAME} diamondC_1x1x1_pp-vmc_gaussian_sdj-1-16 C_Diamond C_Diamond.h5 LAST_TEST_NAME )

  # LCAO and Bspline test should use the same references.
  LIST(APPEND DIAMOND_SCALARS "totenergy" "-10.495941 0.0065")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "11.51198 0.078")
  LIST(APPEND DIAMOND_SCALARS "potential" "-22.00792 0.079")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-2.68548 0.018")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-12.77566 0.0001")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-7.2011 0.090")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "0.65437 0.024")
  LIST(APPEND DIAMOND_SCALARS "samples" "16000 0.0")

  SET(THIS_TEST_NAME "${PYSCF_TEST_NAME}-diamondC_1x1x1_pp-vmc_gaussian_sdj")
  QMC_RUN_AND_CHECK(${THIS_TEST_NAME}
                    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )

  IF(${THIS_TEST_NAME})
    SET_TESTS_PROPERTIES( ${THIS_TEST_NAME}-1-16
                        PROPERTIES DEPENDS ${LAST_TEST_NAME} )
    SET_PROPERTY( TEST ${THIS_TEST_NAME}-1-16
                APPEND PROPERTY LABELS "converter" )
  ENDIF()
            
  RUN_PYSCF_TEST(pyscf-diamond_1x1x1_pp_Bspline ${qmcpack_SOURCE_DIR}/examples/solids/pyscf-inputs
                 diamondC_1x1x1_pp_Bspline PYSCF_TEST_NAME )

  SET_TESTS_PROPERTIES( ${PYSCF_TEST_NAME} PROPERTIES PASS_REGULAR_EXPRESSION "successfully saved to QMCPACK HDF5" )

  SOFTLINK_H5( ${PYSCF_TEST_NAME} diamondC_1x1x1_pp-vmc_gaussian-bspline_sdj-1-16 C_Diamond C_Diamond-Bspline.h5 LAST_TEST_NAME )

  SET(THIS_TEST_NAME "${PYSCF_TEST_NAME}-diamondC_1x1x1_pp-vmc_gaussian-bspline_sdj")
  QMC_RUN_AND_CHECK(${THIS_TEST_NAME}
                    "${qmcpack_SOURCE_DIR}/tests/solids/diamondC_1x1x1-Gaussian_pp"
                    qmc_short-bspline
                    qmc_short-bspline.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )

  IF(${THIS_TEST_NAME})                  
    SET_TESTS_PROPERTIES( ${THIS_TEST_NAME}-1-16
                        PROPERTIES DEPENDS ${LAST_TEST_NAME} )
    SET_PROPERTY( TEST ${THIS_TEST_NAME}-1-16
                APPEND PROPERTY LABELS "converter" )
  ENDIF()
