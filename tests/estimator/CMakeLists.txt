#############################################################
# Add tests to ctest
#############################################################

INCLUDE("${qmcpack_SOURCE_DIR}/CMake/python.cmake")

message("Adding estimator tests for QMCPACK")

IF (NOT TEST_MAX_PROCS )
  SET ( TEST_MAX_PROCS 16 )
  message ("TEST_MAX_PROCS was unset. Set to 16")
ENDIF ()

set(QMCAPP_ERR "")
if((NOT ${TEST_MAX_PROCS}) OR (${TEST_MAX_PROCS} STREQUAL ""))
    set(QMCAPP_ERR "${QMCAPP_ERR}  TEST_MAX_PROCS not set: '${TEST_MAX_PROCS}'\n")
endif()
#
#
if(NOT ${QMCAPP_ERR} STREQUAL "")
  message("${QMCAPP_ERR}  skipping tests")
endif()

set(skinetic_python_reqs numpy;h5py)
CHECK_PYTHON_REQS(skinetic_python_reqs estimator-skinetic add_test)

if (add_test)
  SIMPLE_RUN_AND_CHECK(estimator-skinetic
    "${qmcpack_SOURCE_DIR}/tests/estimator/skinetic"
    vmc.xml
    1 16
    skinetic_check.py
  )
endif()

set(latdev_python_reqs numpy;pandas;h5py)
CHECK_PYTHON_REQS(latdev_python_reqs estimator-latdev add_test)

if (add_test)
  SIMPLE_RUN_AND_CHECK(estimator-latdev
    "${qmcpack_SOURCE_DIR}/tests/estimator/latdev"
    vmc.xml
    1 16
    latdev_check.py
  )
endif()

set(latdev_free_python_reqs numpy;h5py)
CHECK_PYTHON_REQS(latdev_python_reqs estimator-latdev-free add_test)

if (add_test)
  SIMPLE_RUN_AND_CHECK(estimator-latdev-free
    "${qmcpack_SOURCE_DIR}/tests/estimator/latdev/free"
    two.xml
    1 16
    flatdev.py
  )
endif()

set(sofk_python_reqs numpy;pandas;h5py)
CHECK_PYTHON_REQS(sofk_python_reqs estimator-sofk add_test)

if (add_test)
  SIMPLE_RUN_AND_CHECK(estimator-sofk_pbyp-properties
    "${qmcpack_SOURCE_DIR}/tests/estimator/sofk"
    pbyp_dat-h5.xml
    1 16
    check_properties_h5dat.py dat-h5_pbyp
  )
  SIMPLE_RUN_AND_CHECK(estimator-sofk_pbyp-collectables
    "${qmcpack_SOURCE_DIR}/tests/estimator/sofk"
    pbyp_dat-h5.xml
    1 16
    check_collectables_h5dat.py dat-h5_pbyp
  )
  SIMPLE_RUN_AND_CHECK(estimator-sofk_allp-properties
    "${qmcpack_SOURCE_DIR}/tests/estimator/sofk"
    allp_dat-h5.xml
    1 16
    check_properties_h5dat.py dat-h5_allp
  )
  SIMPLE_RUN_AND_CHECK(estimator-sofk_allp-collectables
    "${qmcpack_SOURCE_DIR}/tests/estimator/sofk"
    allp_dat-h5.xml
    1 16
    check_collectables_h5dat.py dat-h5_allp
  )
endif()
