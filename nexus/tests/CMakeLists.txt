
# Add tests for Nexus 
INCLUDE("${PROJECT_SOURCE_DIR}/CMake/test_labels.cmake")
INCLUDE("${PROJECT_SOURCE_DIR}/CMake/python.cmake")

SET(NEXUS_PYTHON_REQS numpy)
CHECK_PYTHON_REQS(NEXUS_PYTHON_REQS "" ADD_TEST)

if (ADD_TEST)
  MESSAGE("Adding Nexus tests")
  SET(TESTLIST "")
  EXECUTE_PROCESS(COMMAND ${PROJECT_SOURCE_DIR}/nexus/bin/nxs-test --ctestlist OUTPUT_VARIABLE TESTLIST)
  #MESSAGE("  Nexus tests: ${TESTLIST}")
  FOREACH(TESTNAME ${TESTLIST})
    #message("Adding test ntest_nexus_${TESTNAME}")
    set(NTEST "${qmcpack_SOURCE_DIR}/nexus/bin/nxs-test")
    add_test(NAME ntest_nexus_${TESTNAME} COMMAND ${NTEST} -R ${TESTNAME} --ctest --pythonpath=${PROJECT_SOURCE_DIR}/nexus/lib)
    set_property(TEST ntest_nexus_${TESTNAME} APPEND PROPERTY LABELS "nexus;deterministic")
    SET(TEST_LABELS "")
    ADD_TEST_LABELS( ntest_nexus_${TESTNAME} TEST_LABELS )
  ENDFOREACH()
ELSE()
  MESSAGE("Skipping Nexus tests because numpy is not present in python installation")
ENDIF()
