function(ADD_TEST_LABELS TEST_NAME TEST_LABELS)
  set(TEST_LABELS_TEMP "")
  if (DEFINED TEST_LABELS_${TEST_NAME}_${QMC_CUDA}_${QMC_COMPLEX}_${QMC_MIXED_PRECISION})
    set(TEST_LABELS_TEMP TEST_LABELS_${${TEST_NAME}_${QMC_CUDA}_${QMC_COMPLEX}_${QMC_MIXED_PRECISION}})
  else()
    set(SUCCESS FALSE)
    execute_process(
      COMMAND ${qmcpack_SOURCE_DIR}/tests/scripts/test_labels.py ${TEST_NAME} ${QMC_CUDA} ${QMC_COMPLEX}
              ${QMC_MIXED_PRECISION}
      OUTPUT_VARIABLE TEST_LABELS_TEMP
      RESULT_VARIABLE SUCCESS)
    if(${SUCCESS} STREQUAL "0")
      set(TEST_LABELS_${TEST_NAME}_${QMC_CUDA}_${QMC_COMPLEX}_${QMC_MIXED_PRECISION} ${TEST_LABELS_TEMP} CACHE INTERNAL "for internal use only; do not modify")
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
