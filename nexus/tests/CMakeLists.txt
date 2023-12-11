# Add tests for Nexus
include("${PROJECT_SOURCE_DIR}/CMake/test_labels.cmake")
include("${PROJECT_SOURCE_DIR}/CMake/python.cmake")

check_python_reqs(numpy "nexus base" ADD_TEST)

if(ADD_TEST)
  message("Adding Nexus tests")
  file(TIMESTAMP ${qmcpack_SOURCE_DIR}/nexus/bin/nxs-test NEXUS_TESTLIST_TIMESTAMP)
  if(NOT "${NEXUS_TESTLIST_TIMESTAMP}" STREQUAL "${CACHE_NEXUS_TESTLIST_TIMESTAMP}")
    execute_process(COMMAND ${Python3_EXECUTABLE} ${qmcpack_SOURCE_DIR}/nexus/bin/nxs-test --ctestlist
                    OUTPUT_VARIABLE NEXUS_TESTLIST)
    set(CACHE_NEXUS_TESTLIST
        ${NEXUS_TESTLIST}
        CACHE INTERNAL "NEXUS_TESTLIST cache variable" FORCE)
    set(CACHE_NEXUS_TESTLIST_TIMESTAMP
        ${NEXUS_TESTLIST_TIMESTAMP}
        CACHE INTERNAL "Timestamp used to validate NEXUS_TESTLIST cache variables" FORCE)
  else()
    set(NEXUS_TESTLIST ${CACHE_NEXUS_TESTLIST})
  endif()
  foreach(TESTNAME ${NEXUS_TESTLIST})
    #message("Adding test ntest_nexus_${TESTNAME}")
    set(NTEST "${qmcpack_SOURCE_DIR}/nexus/bin/nxs-test")
    add_test(NAME ntest_nexus_${TESTNAME} COMMAND ${Python3_EXECUTABLE} ${NTEST} -R ${TESTNAME}\$ --ctest
                                                  --pythonpath=${PROJECT_SOURCE_DIR}/nexus/lib:$ENV{PYTHONPATH})
    set_tests_properties(ntest_nexus_${TESTNAME} PROPERTIES ENVIRONMENT PYTHONPATH=$ENV{PYTHONPATH})
    set_property(
      TEST ntest_nexus_${TESTNAME}
      APPEND
      PROPERTY LABELS "nexus;deterministic")
    set(TEST_LABELS "")
    add_test_labels(ntest_nexus_${TESTNAME} TEST_LABELS)
    # Use a resource lock to avoid potential race condition in copying example inputs
    if(${TESTNAME} MATCHES ".*example.*")
      set_tests_properties(ntest_nexus_${TESTNAME} PROPERTIES RESOURCE_LOCK nexus_examples_resource)
    endif()
  endforeach()
else()
  message("Skipping Nexus tests because numpy is not present in python installation")
endif()
