if(NOT MIXED_PRECISION)
  simple_run_and_check(
    BASE_NAME deterministic-heg2d-4_ae-deriv
    BASE_DIR "${qmcpack_SOURCE_DIR}/tests/heg/heg2d_4"
    INPUT_FILE deriv.xml
    CHECK_SCRIPT check_deriv.py)

  simple_run_and_check(
    BASE_NAME deterministic-heg2d-4_ae-grad_lap
    BASE_DIR "${qmcpack_SOURCE_DIR}/tests/heg/heg2d_4"
    INPUT_FILE grad_lap.xml
    CHECK_SCRIPT check_grad_lap.py)

else()
  message(VERBOSE "Skipping eg2d derivative tests in mixed precision (QMC_MIXED_PRECISION=1)")
endif(NOT MIXED_PRECISION)

get_property(CURRENT_TESTS DIRECTORY PROPERTY TESTS)
add_labels_for_tests("${CURRENT_TESTS}")
