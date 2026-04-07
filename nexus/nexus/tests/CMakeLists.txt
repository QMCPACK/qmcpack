# Add tests for Nexus
include("${PROJECT_SOURCE_DIR}/CMake/test_labels.cmake")
include("${PROJECT_SOURCE_DIR}/CMake/python.cmake")

check_python_reqs("numpy;pytest" "nexus base" ADD_TEST)

if(ADD_TEST)
  message("Generating Nexus test list")
  execute_process(COMMAND ${Python3_EXECUTABLE} -m pytest --collect-only -qq ${qmcpack_SOURCE_DIR}/nexus/nexus/tests
                  # Remove test count per module, join lines, remove file path prefix
                  COMMAND sed -E ":a \$!N ; s/: \[0-9\]+\\n/\\;/ ; s/[nexus\\/]{0,6}nexus\\/tests\\/// ; ta"
                  # Remove trailing newline from sed
                  COMMAND tr -d "\\n"
                  RESULT_VARIABLE NEXUS_TESTLIST_RESULT OUTPUT_VARIABLE NEXUS_TESTLIST)
  if(NOT NEXUS_TESTLIST_RESULT EQUAL 0)
    message(FATAL_ERROR "Command line '${Python3_EXECUTABLE} -m pytest --collect-only -qq ${qmcpack_SOURCE_DIR}/nexus/nexus/tests | sed -E ':a $!N ; s/: [0-9]+\\n\\;// ; s/nexus\\/tests\\/test_// ; ta' | tr -d '\\n'' failed!")
  endif()

  set(CACHE_NEXUS_TESTLIST
      ${NEXUS_TESTLIST}
      CACHE INTERNAL "NEXUS_TESTLIST cache variable" FORCE)
  if(NEXUS_TESTLIST)
    message("Adding Nexus tests")
  else()
    message("Nexus test list is empty. No Nexus tests added.")
  endif()
  foreach(TESTFILE ${NEXUS_TESTLIST})
    string(REPLACE "test_" "" TESTNAME ${TESTFILE})
    string(REPLACE ".py" "" TESTNAME ${TESTNAME})
    #message("Adding test ntest_nexus_${TESTNAME}")
    set(NTEST "pytest")
    IF(ENABLE_PYCOV)
      add_test(NAME ntest_nexus_${TESTNAME}
               COMMAND ${Python3_EXECUTABLE} -m ${NTEST} ${qmcpack_SOURCE_DIR}/nexus/nexus/tests/${TESTFILE} -c ${qmcpack_SOURCE_DIR}/nexus/pyproject.toml
                                                                                   --cov=nexus
                                                                                   --cov-report=xml
                                                                                   --cov-append
                                                                                   --cov-config=${qmcpack_SOURCE_DIR}/nexus/pyproject.toml
              )

    ELSE()
      add_test(NAME ntest_nexus_${TESTNAME}
               COMMAND ${Python3_EXECUTABLE} -m ${NTEST} ${qmcpack_SOURCE_DIR}/nexus/nexus/tests/${TESTFILE} -c ${qmcpack_SOURCE_DIR}/nexus/pyproject.toml
              )
    ENDIF()
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
  message("Skipping Nexus tests because numpy and/or pytest is not present in python installation")
endif()
