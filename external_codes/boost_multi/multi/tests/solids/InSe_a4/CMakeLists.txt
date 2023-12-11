set(TEST_ADDED FALSE)
set(TEST_LABELS "")

run_qmc_app(
  deterministic-InSe_a4_slab
  "${qmcpack_SOURCE_DIR}/tests/solids/InSe_a4"
  1
  1
  TEST_ADDED
  TEST_LABELS
  InSe-S1.xml)
