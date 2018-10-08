
# Add tests for Nexus 

INCLUDE("${qmcpack_SOURCE_DIR}/CMake/python.cmake")

SET(NEXUS_PYTHON_REQS numpy)
CHECK_PYTHON_REQS(NEXUS_PYTHON_REQS "" ADD_TEST)

if (ADD_TEST)
  MESSAGE("Adding Nexus tests")
  SET(TESTLIST "")
  EXECUTE_PROCESS(COMMAND ${PROJECT_SOURCE_DIR}/nexus/executables/ntest --ctestlist OUTPUT_VARIABLE TESTLIST)
  #MESSAGE("  Nexus tests: ${TESTLIST}")
  FOREACH(TESTNAME ${TESTLIST})
    message("Adding test ntest_nexus_${TESTNAME}")
    set(NTEST "${qmcpack_SOURCE_DIR}/nexus/executables/ntest")
    add_test(NAME ntest_nexus_${TESTNAME} COMMAND ${NTEST} -R ${TESTNAME} --ctest --pythonpath=${PROJECT_SOURCE_DIR}/nexus/library)
    set_property(TEST ntest_nexus_${TESTNAME} APPEND PROPERTY LABELS "nexus")
    SET(TEST_LABELS "")
    EXECUTE_PROCESS(COMMAND ${PROJECT_SOURCE_DIR}/tests/scripts/test_labels.py ntest_nexus_${TESTNAME} ${QMC_CUDA} ${ENABLE_SOA} ${QMC_COMPLEX} ${QMC_MIXED_PRECISION} OUTPUT_VARIABLE TEST_LABELS RESULT_VARIABLE SUCCESS)
    SET_PROPERTY(TEST ntest_nexus_${TESTNAME} APPEND PROPERTY LABELS ${TEST_LABELS})
  ENDFOREACH()
ELSE()
  MESSAGE("Skipping Nexus tests because numpy is not present in python installation")
ENDIF()
