
  LIST(APPEND DIAMOND_SCALARS "totenergy" "-10.491445 0.002")
  LIST(APPEND DIAMOND_SCALARS "kinetic" "11.434392 0.018")
  LIST(APPEND DIAMOND_SCALARS "potential" "-21.925822 0.018")
  LIST(APPEND DIAMOND_SCALARS "eeenergy" "-2.700534 0.0043")
  LIST(APPEND DIAMOND_SCALARS "ionion" "-12.77566752 0.001")
  LIST(APPEND DIAMOND_SCALARS "localecp" "-7.046740 0.020")
  LIST(APPEND DIAMOND_SCALARS "nonlocalecp" "0.597119 0.0056")
  LIST(APPEND DIAMOND_SCALARS "samples" "128000 0.0")
#  LIST(APPEND DIAMOND_SCALARS "flux" "0.0 0.4")

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short
                    qmc_short.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )

IF(ENABLE_SOA OR QMC_CUDA)
  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_hybridrep_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short
                    qmc_short_hybridrep.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )
ELSE()
  MESSAGE("Skipping diamondC_1x1x1_hybridrep_pp-vmc_sdj test because neiher SOA nor CUDA is enabled")
ENDIF()

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-vmc_sdj-meshf
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short-meshf
                    qmc_short-meshf.in.xml
                    1 16
                    TRUE
                    0 DIAMOND_SCALARS # VMC
                    )

# Reference OPT run in qmc-ref
  LIST(APPEND DIAMOND_OPT_SCALARS "totenergy"   "-10.49370 0.0024")
  LIST(APPEND DIAMOND_OPT_SCALARS "kinetic"     "11.5504 0.021")
  LIST(APPEND DIAMOND_OPT_SCALARS "localecp"    "-7.1979 0.023")
  LIST(APPEND DIAMOND_OPT_SCALARS "nonlocalecp" "0.61540 0.0068")
  LIST(APPEND DIAMOND_OPT_SCALARS "samples"     "128000 0.0")

  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-opt_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short
                    qmc_short_opt.in.xml
                    1 16
                    ${MP_SUCCESS}
                    3 DIAMOND_OPT_SCALARS # VMC
                    )

# Reference DMC run in qmc-ref "-10.531583 0.000265"
  LIST(APPEND DIAMOND_DMC_SCALARS "totenergy" "-10.531583 0.0024")
  QMC_RUN_AND_CHECK(short-diamondC_1x1x1_pp-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_short_vmc_dmc
                    qmc_short_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 DIAMOND_DMC_SCALARS # DMC
                    )

# Coverage test - shorter version of dmc test
  SET(FULL_NAME "coverage-diamondC_1x1x1_pp-vmc_sdj-1-1")
  SET(TEST_ADDED FALSE)
  RUN_QMC_APP(${FULL_NAME}
              "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
              1 1
              TEST_ADDED
              "qmc_short_vmc_dmc_coverage.in.xml")

  IF (TEST_ADDED)
    SET_PROPERTY(TEST ${FULL_NAME} APPEND PROPERTY LABELS "coverage")
  ENDIF()

# Test time limit
  SET(FULL_NAME "cpu_limit-diamondC_1x1x1_pp-vmc-1-1")
  SET(TEST_ADDED FALSE)
  RUN_QMC_APP(${FULL_NAME}
              "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
              1 1
              TEST_ADDED
              qmc_cpu_limit_vmc.xml)

  IF (TEST_ADDED)
    SET_PROPERTY(TEST ${FULL_NAME} APPEND PROPERTY TIMEOUT 180)
    SET_PROPERTY(TEST ${FULL_NAME} APPEND PROPERTY PASS_REGULAR_EXPRESSION "Time limit reached for")
  ENDIF()

  SET(FULL_NAME "cpu_limit-diamondC_1x1x1_pp-vmc-4-4")
  SET(TEST_ADDED FALSE)
  RUN_QMC_APP(${FULL_NAME}
              "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
              4 4
              TEST_ADDED
              qmc_cpu_limit_vmc.xml)

  IF (TEST_ADDED)
    SET_PROPERTY(TEST ${FULL_NAME} APPEND PROPERTY TIMEOUT 120)
    SET_PROPERTY(TEST ${FULL_NAME} APPEND PROPERTY PASS_REGULAR_EXPRESSION "Time limit reached for")
  ENDIF()

  SET(FULL_NAME "cpu_limit-diamondC_1x1x1_pp-dmc-1-1")
  SET(TEST_ADDED FALSE)
  RUN_QMC_APP(${FULL_NAME}
              "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
              1 1
              TEST_ADDED
              qmc_cpu_limit_dmc.xml)

  IF (TEST_ADDED)
    SET_PROPERTY(TEST ${FULL_NAME} APPEND PROPERTY TIMEOUT 240)
    SET_PROPERTY(TEST ${FULL_NAME} APPEND PROPERTY PASS_REGULAR_EXPRESSION "Time limit reached for")
  ENDIF()

  SET(FULL_NAME "cpu_limit-diamondC_1x1x1_pp-dmc-4-4")
  SET(TEST_ADDED FALSE)
  RUN_QMC_APP(${FULL_NAME}
              "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
              4 4
              TEST_ADDED
              qmc_cpu_limit_dmc.xml)

  IF (TEST_ADDED)
    SET_PROPERTY(TEST ${FULL_NAME} APPEND PROPERTY TIMEOUT 120)
    SET_PROPERTY(TEST ${FULL_NAME} APPEND PROPERTY PASS_REGULAR_EXPRESSION "Time limit reached for")
  ENDIF()


#
# Long tests
#

  LIST(APPEND LONG_DIAMOND_SCALARS "totenergy" "-10.491445 0.000065")
  LIST(APPEND LONG_DIAMOND_SCALARS "samples" "122880000 0.0")
  LIST(APPEND LONG_DIAMOND_SCALARS "flux" "0.0 0.03")

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_long
                    qmc_long.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )

  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-vmc_sdj-meshf
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_long-meshf
                    qmc_long-meshf.in.xml
                    1 16
                    TRUE
                    0 LONG_DIAMOND_SCALARS # VMC
                    )

# Reference DMC run in qmc-ref "-10.531583 0.000265"
  LIST(APPEND LONG_DIAMOND_DMC_SCALARS "totenergy" "-10.531583 0.000815")
  QMC_RUN_AND_CHECK(long-diamondC_1x1x1_pp-dmc_sdj
                    "${CMAKE_SOURCE_DIR}/tests/solids/diamondC_1x1x1_pp"
                    qmc_long_vmc_dmc
                    qmc_long_vmc_dmc.in.xml
                    1 16
                    TRUE
                    1 LONG_DIAMOND_DMC_SCALARS # DMC
                    )

