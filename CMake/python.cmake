# Support functions for handling python scripts

# Test whether a python modules is present
#   MODULE_NAME - input, name of module to test for
#   MODULE_PRESENT - output - True/False based on success of the import
FUNCTION (TEST_PYTHON_MODULE MODULE_NAME MODULE_PRESENT)
  EXECUTE_PROCESS(
    COMMAND python ${qmcpack_SOURCE_DIR}/utils/test_import.py ${MODULE_NAME}
    OUTPUT_VARIABLE TMP_OUTPUT_VAR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  SET(${MODULE_PRESENT} ${TMP_OUTPUT_VAR} PARENT_SCOPE)
ENDFUNCTION()

# Test python module prerequisites for a particular test script
#   module_list - input - list of module names
#   test_name - input - name of test (used for missing module message)
#                     - use empty string to silence output
#   add_test - output - true if all modules are present, false otherwise
FUNCTION(CHECK_PYTHON_REQS module_list test_name add_test)
  set(${add_test} true PARENT_SCOPE)
  foreach(python_module IN LISTS ${module_list})
    TEST_PYTHON_MODULE(${python_module} has_python_module)
    if (NOT(has_python_module))
      if (test_name)
        MESSAGE("Missing python module ${python_module}, not adding test ${test_name}")
      endif()
      set(${add_test} false PARENT_SCOPE)
    endif()
  endforeach()
ENDFUNCTION()

