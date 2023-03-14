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

check_python_reqs("numpy;h5py" estimator-skinetic add_test)

if(add_test)
  simple_run_and_check(estimator-skinetic "${qmcpack_SOURCE_DIR}/tests/estimator/skinetic" vmc.xml 1 16
                       skinetic_check.py)
endif()

check_python_reqs("numpy;pandas;h5py" estimator-latdev add_test)

if(add_test)
  simple_run_and_check(estimator-latdev "${qmcpack_SOURCE_DIR}/tests/estimator/latdev" vmc.xml 1 16 latdev_check.py)
endif()

check_python_reqs("numpy;h5py" estimator-latdev-free add_test)

if(add_test)
  simple_run_and_check(estimator-latdev-free "${qmcpack_SOURCE_DIR}/tests/estimator/latdev/free" two.xml 1 16
                       flatdev.py)
endif()

check_python_reqs("numpy;pandas;h5py" estimator-sofk add_test)

if(add_test)
  simple_run_and_check(
    estimator-sofk_pbyp-properties
    "${qmcpack_SOURCE_DIR}/tests/estimator/sofk"
    pbyp_dat-h5.xml
    1
    16
    check_properties_h5dat.py
    dat-h5_pbyp)
  simple_run_and_check(
    estimator-sofk_pbyp-collectables
    "${qmcpack_SOURCE_DIR}/tests/estimator/sofk"
    pbyp_dat-h5.xml
    1
    16
    check_collectables_h5dat.py
    dat-h5_pbyp)
  simple_run_and_check(
    estimator-sofk_allp-properties
    "${qmcpack_SOURCE_DIR}/tests/estimator/sofk"
    allp_dat-h5.xml
    1
    16
    check_properties_h5dat.py
    dat-h5_allp)
  simple_run_and_check(
    estimator-sofk_allp-collectables
    "${qmcpack_SOURCE_DIR}/tests/estimator/sofk"
    allp_dat-h5.xml
    1
    16
    check_collectables_h5dat.py
    dat-h5_allp)
  simple_run_and_check(estimator-acforce-legacy "${qmcpack_SOURCE_DIR}/tests/estimator/acforce" vmc.legacy.in.xml 1 16
                       check_forces.py)
  simple_run_and_check(
    estimator-acforce-fast
    "${qmcpack_SOURCE_DIR}/tests/estimator/acforce"
    vmc.fast.in.xml
    1
    16
    check_forces.py
    vmc.fast.s000.scalar.dat)
endif()
