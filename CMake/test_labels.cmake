function(ADD_LABELS_FOR_TESTS TEST_NAMES)
  list(LENGTH TEST_NAMES num_tests)
  if(num_tests EQUAL 0)
    return()
  endif()

  set(TEST_LABELS_TEMP "")
  set(SUCCESS FALSE)
  execute_process(
    COMMAND ${Python3_EXECUTABLE} ${qmcpack_SOURCE_DIR}/tests/scripts/test_labels.py "${TEST_NAMES}"
            ${QMC_COMPLEX} ${QMC_MIXED_PRECISION}
    OUTPUT_VARIABLE TEST_LABELS_TEMP
    RESULT_VARIABLE SUCCESS
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(NOT ${SUCCESS} STREQUAL "0")
    message("Warning: test labeling failed.  Test labeling error output:\n${TEST_LABELS_TEMP}")
    return()
  endif()

  string(REPLACE ";" "\\;" TEST_LABELS_TEMP_ESCAPED "${TEST_LABELS_TEMP}")
  string(REPLACE "\n" ";" TEST_LABELS_LINES "${TEST_LABELS_TEMP_ESCAPED}")

  list(LENGTH TEST_LABELS_LINES num_lines)
  if(NOT num_tests EQUAL num_lines)
    message("Warning: test labeling mismatch. Expected ${num_tests} labels, got ${num_lines}.")
    return()
  endif()

  set(i 0)
  foreach(TEST_NAME IN LISTS TEST_NAMES)
    list(GET TEST_LABELS_LINES ${i} TEST_LABELS_LOCAL)
    if(TEST_NAME MATCHES "-r[0-9][0-9]?-t[0-9][0-9]?$")
      list(REMOVE_ITEM TEST_LABELS_LOCAL unstable)
    endif()
    if(TEST_LABELS_LOCAL)
      set_property(
        TEST ${TEST_NAME}
        APPEND
        PROPERTY LABELS ${TEST_LABELS_LOCAL})
    endif()
    math(EXPR i "${i} + 1")
  endforeach()
endfunction()
