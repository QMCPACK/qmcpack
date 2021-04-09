# Support functions for handling python scripts

# Test whether a python modules is present
#   MODULE_NAME - input, name of module to test for
#   MODULE_PRESENT - output - True/False based on success of the import
FUNCTION (TEST_PYTHON_MODULE MODULE_NAME MODULE_PRESENT)
  MESSAGE_VERBOSE("Checking import python module ${MODULE_NAME}")
  EXECUTE_PROCESS(
    COMMAND ${qmcpack_SOURCE_DIR}/tests/scripts/test_import.py ${MODULE_NAME}
    OUTPUT_VARIABLE TMP_OUTPUT_VAR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  SET(${MODULE_PRESENT} ${TMP_OUTPUT_VAR} CACHE BOOL "" FORCE)
ENDFUNCTION()

# Test python module prerequisites for a particular test script
#   module_list - input - list of module names
#   test_name - input - name of test (used for missing module message)
#                     - use empty string to silence output
#   add_test - output - true if all modules are present, false otherwise
FUNCTION(CHECK_PYTHON_REQS module_list test_name add_test)
  set(${add_test} true PARENT_SCOPE)
  foreach(python_module IN LISTS ${module_list})
    string(TOUPPER ${python_module} module_uppercase)
    set(cached_result HAS_${module_uppercase}_PYTHON_MODULE)
    if (NOT DEFINED ${cached_result})
      TEST_PYTHON_MODULE(${python_module} ${cached_result})
    endif ()
    if (NOT(${cached_result}))
      if (test_name)
        MESSAGE("Missing python module ${python_module}, not adding test ${test_name}")
      endif()
      set(${add_test} false PARENT_SCOPE)
    endif()
  endforeach()
ENDFUNCTION()

