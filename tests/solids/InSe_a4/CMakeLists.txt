set(TEST_ADDED FALSE)

run_qmc_app(
  deterministic-InSe_a4_slab
  "${qmcpack_SOURCE_DIR}/tests/solids/InSe_a4"
  1
  1
  TEST_ADDED
  InSe-S1.xml)

get_property(CURRENT_TESTS DIRECTORY PROPERTY TESTS)
add_labels_for_tests("${CURRENT_TESTS}")
