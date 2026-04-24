# Add tests for Nexus
include("${PROJECT_SOURCE_DIR}/CMake/test_labels.cmake")
include("${PROJECT_SOURCE_DIR}/CMake/python.cmake")

check_python_reqs("numpy;pytest" "nexus base" ADD_TEST)

if(ADD_TEST)
  message("Generating Nexus test list")
  execute_process(COMMAND ${Python3_EXECUTABLE} -m pytest --collect-only -qq ${qmcpack_SOURCE_DIR}/nexus/nexus/tests
                  RESULT_VARIABLE NEXUS_TESTLIST_RESULT OUTPUT_VARIABLE NEXUS_TESTLIST)
  if(NOT NEXUS_TESTLIST_RESULT EQUAL 0)
    message(FATAL_ERROR "Command line ' ${Python3_EXECUTABLE} -m pytest --collect-only -qq ${qmcpack_SOURCE_DIR}/nexus/nexus/tests ' failed!")
  endif()

  string(REPLACE "\n" ";" NEXUS_TESTLIST "${NEXUS_TESTLIST}") # Join into a single line with semicolon separators
  list(REMOVE_ITEM NEXUS_TESTLIST "") # remove empty entries
  set(CACHE_NEXUS_TESTLIST
      ${NEXUS_TESTLIST}
      CACHE INTERNAL "NEXUS_TESTLIST cache variable" FORCE)
  if(NEXUS_TESTLIST)
    message("Adding Nexus tests")
  else()
    message("Nexus test list is empty. No Nexus tests added.")
  endif()
  set(NTEST "pytest")
  if(ENABLE_PYCOV)
    add_test(NAME ntest_nexus_coverage_all
             COMMAND ${Python3_EXECUTABLE} -m ${NTEST} ${qmcpack_SOURCE_DIR}/nexus/nexus/tests -c ${qmcpack_SOURCE_DIR}/nexus/pyproject.toml
                                                                                   --cov=nexus
                                                                                   --cov-report=xml
                                                                                   --cov-config=${qmcpack_SOURCE_DIR}/nexus/pyproject.toml
    )
    set_tests_properties(ntest_nexus_coverage_all PROPERTIES ENVIRONMENT PYTHONPATH=$ENV{PYTHONPATH})
    set_property(
      TEST ntest_nexus_coverage_all
      APPEND
      PROPERTY LABELS "nexus;deterministic")
    set(TEST_LABELS "")
    add_test_labels(ntest_nexus_coverage_all TEST_LABELS)
  else()
    foreach(TESTFILE IN LISTS NEXUS_TESTLIST)
      if("${TESTFILE}" STREQUAL "") # Empty lines should be skipped. They usually precede warnings.
        break()
      else()
        # If the last attempt to cut off the collection early didn't catch the warnings, this will.
        # The warning line looks like
        # =================================== warnings summary ===================================
        string(FIND "${TESTFILE}" "=============================" WARNLINE)
        if(NOT WARNLINE EQUAL -1)
          break()
        endif()
      endif()
      string(REGEX REPLACE "\.py:.*[0-9]" "" TESTFILE "${TESTFILE}")
      string(REPLACE "nexus/tests/test_" "" TESTNAME "${TESTFILE}")
      #message("Adding test ntest_nexus_${TESTNAME}")

      add_test(NAME ntest_nexus_${TESTNAME}
               COMMAND ${Python3_EXECUTABLE} -m ${NTEST} ${qmcpack_SOURCE_DIR}/nexus/${TESTFILE}.py -c ${qmcpack_SOURCE_DIR}/nexus/pyproject.toml
              )
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
  endif()
else()
  message(WARNING "Skipping Nexus tests because python numpy package and/or pytest were not found")
endif()
