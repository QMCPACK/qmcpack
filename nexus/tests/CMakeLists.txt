
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
    message("Adding ntest_nexus_${TESTNAME}")
    set(NTEST "${qmcpack_SOURCE_DIR}/nexus/executables/ntest")
    add_test(NAME ntest_nexus_${TESTNAME} COMMAND ${NTEST} -R ${TESTNAME} --ctest)
    set_tests_properties(ntest_nexus_${TESTNAME} PROPERTIES LABELS "nexus")
  ENDFOREACH()
ELSE()
  MESSAGE("Skipping Nexus tests because numpy is not present in python installation")
ENDIF()
