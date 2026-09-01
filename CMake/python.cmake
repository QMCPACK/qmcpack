# Support functions for handling python scripts

# Test whether a python modules is present
#   MODULE_NAME - input, name of module to test for
#   MODULE_PRESENT - output - True/False based on success of the import
function(TEST_PYTHON_MODULE MODULE_NAME MODULE_PRESENT)
  message(VERBOSE "Checking import python module ${MODULE_NAME}")
  execute_process(
    COMMAND ${Python3_EXECUTABLE} ${qmcpack_SOURCE_DIR}/tests/scripts/test_import.py ${MODULE_NAME}
    OUTPUT_VARIABLE TMP_OUTPUT_VAR
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(${MODULE_PRESENT}
      ${TMP_OUTPUT_VAR}
      CACHE BOOL "" FORCE)
endfunction()

# Ensure we have a compatible version of Pytest for Nexus's testing.
function(CHECK_PYTEST_VERSION pytest_version_ok)
  execute_process(COMMAND ${Python3_EXECUTABLE} -m pytest --version
                  RESULT_VARIABLE PYTEST_VERSION_RESULT
                  OUTPUT_VARIABLE PYTEST_VERSION
                  ERROR_VARIABLE PYTEST_VERSION_ERROR
                  OUTPUT_STRIP_TRAILING_WHITESPACE
                  ERROR_STRIP_TRAILING_WHITESPACE)
  if(PYTEST_VERSION)
    string(REPLACE "pytest " "" PYTEST_VERSION ${PYTEST_VERSION})
  else()
    # Old versions of pytest put the version output in stderr for some reason
    string(REPLACE "pytest " "" PYTEST_VERSION ${PYTEST_VERSION_ERROR})
  endif()
  if(PYTEST_VERSION VERSION_GREATER_EQUAL "6.2.4")
    message(STATUS "Found compatible Pytest version ${PYTEST_VERSION}")
    set(${pytest_version_ok}
        TRUE
        CACHE BOOL "" FORCE)
  else()
    message(STATUS "Incompatible Pytest version ${PYTEST_VERSION} found, tests will not be added. Required version >= 6.2.4 ")
    set(${pytest_version_ok}
        FALSE
        CACHE BOOL "" FORCE)
  endif()
endfunction()

# Test python module prerequisites for a particular test script
#   module_list - input - list of module names, for example "numpy;h5py" (including the quotation marks)
#   test_name - input - name of test (used for missing module message)
#                     - use empty string to silence output
#   add_test - output - true if all modules are present, false otherwise
function(CHECK_PYTHON_REQS module_list test_name add_test)
  if(NOT ${ARGC} EQUAL 3)
    message(FATAL_ERROR "Expecting 3 arguments in CHECK_PYTHON_REQS. Was called with ${ARGC} arguments. "
                        "Most likely quotation marks are missing from the first argument.")
  endif()
  set(${add_test}
      true
      PARENT_SCOPE)
  foreach(python_module IN LISTS module_list)
    string(TOUPPER ${python_module} module_uppercase)
    set(cached_variable_name HAS_${module_uppercase}_PYTHON_MODULE)
    if(NOT DEFINED ${cached_variable_name})
      test_python_module(${python_module} ${cached_variable_name})
    endif()
    if(NOT ${cached_variable_name})
      if(test_name)
        message("Skipping ${test_name} tests because valid module ${python_module} is not found.")
      endif()
      set(${add_test}
          false
          PARENT_SCOPE)
    else()
      if(${python_module} STREQUAL "pytest")
        check_pytest_version(pytest_version_ok)
        if(NOT ${pytest_version_ok})
          message("Pytest version < 6.2.4 found, will not add Nexus tests")
          set(${add_test}
              false
              PARENT_SCOPE)
        endif()
      endif()
    endif()
  endforeach()
endfunction()

