function(ADD_TEST_LABELS TEST_NAME TEST_LABELS)
  set(TEST_LABELS_TEMP "")
  set(TEST_LABELS_UNIQUE_NAME TEST_LABELS_${TEST_NAME}_${QMC_COMPLEX}_${QMC_MIXED_PRECISION})
  if(DEFINED ${TEST_LABELS_UNIQUE_NAME})
    set(TEST_LABELS_TEMP ${${TEST_LABELS_UNIQUE_NAME}})
  else()
    set(SUCCESS FALSE)
    execute_process(
      COMMAND ${Python3_EXECUTABLE} ${qmcpack_SOURCE_DIR}/tests/scripts/test_labels.py ${TEST_NAME}
              ${QMC_COMPLEX} ${QMC_MIXED_PRECISION}
      OUTPUT_VARIABLE TEST_LABELS_TEMP
      RESULT_VARIABLE SUCCESS
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(${SUCCESS} STREQUAL "0")
      set(${TEST_LABELS_UNIQUE_NAME}
          ${TEST_LABELS_TEMP}
          CACHE INTERNAL "for internal use only; do not modify")
    else()
      message("Warning: test labeling failed.  Test labeling error output:\n${TEST_LABELS_TEMP}")
      set(TEST_LABELS_TEMP "")
    endif()
  endif()
  # Remove unstable label from direct execution.
  # It will still be added to statistical child tests.
  set(TEST_LABELS_LOCAL ${TEST_LABELS_TEMP})
  list(REMOVE_ITEM TEST_LABELS_LOCAL unstable)
  set_property(
    TEST ${TEST_NAME}
    APPEND
    PROPERTY LABELS ${TEST_LABELS_LOCAL})
  set(${TEST_LABELS}
      ${TEST_LABELS_TEMP}
      PARENT_SCOPE)
endfunction()

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

  math(EXPR loop_max "${num_tests} - 1")
  foreach(i RANGE ${loop_max})
    list(GET TEST_NAMES ${i} TEST_NAME)
    list(GET TEST_LABELS_LINES ${i} TEST_LABELS_LOCAL)
    list(REMOVE_ITEM TEST_LABELS_LOCAL unstable)
    if(TEST_LABELS_LOCAL)
      set_property(
        TEST ${TEST_NAME}
        APPEND
        PROPERTY LABELS ${TEST_LABELS_LOCAL})
    endif()
  endforeach()
endfunction()
