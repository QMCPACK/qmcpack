#############################################################
# Add tests to ctest
#############################################################

include("${qmcpack_SOURCE_DIR}/CMake/python.cmake")

message("Adding estimator tests for QMCPACK")

if(NOT TEST_MAX_PROCS)
  set(TEST_MAX_PROCS 16)
  message("TEST_MAX_PROCS was unset. Set to 16")
endif()

set(QMCAPP_ERR "")
if((NOT ${TEST_MAX_PROCS}) OR (${TEST_MAX_PROCS} STREQUAL ""))
  set(QMCAPP_ERR "${QMCAPP_ERR}  TEST_MAX_PROCS not set: '${TEST_MAX_PROCS}'\n")
endif()
#
#
if(NOT ${QMCAPP_ERR} STREQUAL "")
  message("${QMCAPP_ERR}  skipping tests")
endif()

check_python_reqs("numpy;h5py" estimator-skinetic add_tests)

if(add_tests)
  simple_run_and_check(
    BASE_NAME deterministic-estimator-skinetic
    BASE_DIR "${qmcpack_SOURCE_DIR}/tests/estimator/skinetic"
    INPUT_FILE vmc.xml
    THREADS 16
    CHECK_SCRIPT skinetic_check.py)
endif()

check_python_reqs("numpy;h5py" estimator-latdev add_tests)

if(add_tests)
  simple_run_and_check(
    BASE_NAME deterministic-estimator-latdev
    BASE_DIR "${qmcpack_SOURCE_DIR}/tests/estimator/latdev"
    INPUT_FILE vmc.xml
    THREADS 16
    CHECK_SCRIPT latdev_check.py)
endif()

check_python_reqs("numpy;h5py" estimator-latdev-free add_tests)

if(add_tests)
  simple_run_and_check(
    BASE_NAME estimator-latdev-free
    BASE_DIR "${qmcpack_SOURCE_DIR}/tests/estimator/latdev/free"
    INPUT_FILE two.xml
    THREADS 16
    CHECK_SCRIPT flatdev.py)
endif()

check_python_reqs("numpy;h5py" estimator-sofk add_tests)

if(add_tests)
  simple_run_and_check(
    BASE_NAME deterministic-estimator-sofk_pbyp-properties
    BASE_DIR "${qmcpack_SOURCE_DIR}/tests/estimator/sofk"
    INPUT_FILE pbyp_dat-h5.xml
    THREADS 16
    CHECK_SCRIPT check_properties_h5dat.py
    CHECK_ARGS dat-h5_pbyp)
  simple_run_and_check(
    BASE_NAME deterministic-estimator-sofk_pbyp-collectables
    BASE_DIR "${qmcpack_SOURCE_DIR}/tests/estimator/sofk"
    INPUT_FILE pbyp_dat-h5.xml
    THREADS 16
    CHECK_SCRIPT check_collectables_h5dat.py
    CHECK_ARGS dat-h5_pbyp)
  simple_run_and_check(
    BASE_NAME deterministic-estimator-sofk_allp-properties
    BASE_DIR "${qmcpack_SOURCE_DIR}/tests/estimator/sofk"
    INPUT_FILE allp_dat-h5.xml
    THREADS 16
    CHECK_SCRIPT check_properties_h5dat.py
    CHECK_ARGS dat-h5_allp)
  simple_run_and_check(
    BASE_NAME deterministic-estimator-sofk_allp-collectables
    BASE_DIR "${qmcpack_SOURCE_DIR}/tests/estimator/sofk"
    INPUT_FILE allp_dat-h5.xml
    THREADS 16
    CHECK_SCRIPT check_collectables_h5dat.py
    CHECK_ARGS dat-h5_allp)
  simple_run_and_check(
    BASE_NAME estimator-acforce-legacy
    BASE_DIR "${qmcpack_SOURCE_DIR}/tests/estimator/acforce"
    INPUT_FILE vmc.legacy.in.xml
    THREADS 16
    CHECK_SCRIPT check_forces.py)
  simple_run_and_check(
    BASE_NAME estimator-acforce-fast
    BASE_DIR "${qmcpack_SOURCE_DIR}/tests/estimator/acforce"
    INPUT_FILE vmc.fast.in.xml
    THREADS 16
    CHECK_SCRIPT check_forces.py
    CHECK_ARGS vmc.fast.s000.scalar.dat)
endif()

get_property(CURRENT_TESTS DIRECTORY PROPERTY TESTS)
add_labels_for_tests("${CURRENT_TESTS}")
