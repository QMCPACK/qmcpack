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

# Test python module prerequisites for a particular test script
#   module_list - input - list of module names
#   test_name - input - name of test (used for missing module message)
#                     - use empty string to silence output
#   add_test - output - true if all modules are present, false otherwise
function(CHECK_PYTHON_REQS module_list test_name add_test)
  set(${add_test}
      true
      PARENT_SCOPE)
  foreach(python_module IN LISTS ${module_list})
    string(TOUPPER ${python_module} module_uppercase)
    set(cached_variable_name HAS_${module_uppercase}_PYTHON_MODULE)
    if(NOT DEFINED ${cached_variable_name})
      test_python_module(${python_module} ${cached_variable_name})
    endif()
    if(NOT ${cached_variable_name})
      if(test_name)
        message("Missing python module ${python_module}, not adding test ${test_name}")
      endif()
      set(${add_test}
          false
          PARENT_SCOPE)
    endif()
  endforeach()
endfunction()
