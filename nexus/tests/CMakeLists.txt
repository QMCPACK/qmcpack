# Add tests for Nexus
include("${PROJECT_SOURCE_DIR}/CMake/test_labels.cmake")
include("${PROJECT_SOURCE_DIR}/CMake/python.cmake")

set(NEXUS_PYTHON_REQS numpy)
check_python_reqs(NEXUS_PYTHON_REQS "" ADD_TEST)

if(ADD_TEST)
  message("Adding Nexus tests")
  set(TESTLIST "")
  execute_process(COMMAND ${PROJECT_SOURCE_DIR}/nexus/bin/nxs-test --ctestlist OUTPUT_VARIABLE TESTLIST)
  #MESSAGE("  Nexus tests: ${TESTLIST}")
  foreach(TESTNAME ${TESTLIST})
    #message("Adding test ntest_nexus_${TESTNAME}")
    set(NTEST "${qmcpack_SOURCE_DIR}/nexus/bin/nxs-test")
    add_test(NAME ntest_nexus_${TESTNAME} COMMAND ${NTEST} -R ${TESTNAME} --ctest
                                                  --pythonpath=${PROJECT_SOURCE_DIR}/nexus/lib)
    set_property(
      TEST ntest_nexus_${TESTNAME}
      APPEND
      PROPERTY LABELS "nexus;deterministic")
    set(TEST_LABELS "")
    add_test_labels(ntest_nexus_${TESTNAME} TEST_LABELS)
    # Use a resource lock to avoid potential race condition in copying example inputs
    if (${TESTNAME} MATCHES ".*example.*")
      set_tests_properties(ntest_nexus_${TESTNAME} PROPERTIES RESOURCE_LOCK nexus_examples_resource)
    endif()
  endforeach()
else()
  message("Skipping Nexus tests because numpy is not present in python installation")
endif()
