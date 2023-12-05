# functions for pyscf workflow test

function(ADD_PYSCF_TEST TESTNAME TEST_BINARY WORKDIR TEST_INPUT)
  add_test(NAME ${TESTNAME} COMMAND ${TEST_BINARY} ${TEST_INPUT})

  set_tests_properties(${TESTNAME} PROPERTIES WORKING_DIRECTORY ${WORKDIR})
  set_property(
    TEST ${TESTNAME}
    APPEND
    PROPERTY LABELS "converter")
endfunction()

function(RUN_PYSCF_TEST BASE_NAME SRC_DIR TEST_INPUT_PREFIX TEST_NAME)
  set(${TEST_NAME}
      ${BASE_NAME}
      PARENT_SCOPE)
  set(MY_WORKDIR ${CMAKE_CURRENT_BINARY_DIR}/${BASE_NAME})
  message(VERBOSE "Adding test ${BASE_NAME}")
  copy_directory("${SRC_DIR}" "${MY_WORKDIR}")
  add_pyscf_test(${BASE_NAME} python ${MY_WORKDIR} ${TEST_INPUT_PREFIX}.py)
endfunction()

function(SOFTLINK_H5 SOURCE TARGET PREFIX FILENAME TEST_NAME)
  set(${TEST_NAME}
      "LINK_${SOURCE}_TO_${TARGET}"
      PARENT_SCOPE)
  add_test(NAME LINK_${SOURCE}_TO_${TARGET} COMMAND ${qmcpack_SOURCE_DIR}/tests/scripts/clean_and_link_h5.sh
                                                    ${SOURCE}/${PREFIX}.h5 ${SOURCE}-${TARGET}/${FILENAME})
  set_tests_properties(LINK_${SOURCE}_TO_${TARGET} PROPERTIES DEPENDS ${SOURCE})
  set_property(TEST LINK_${SOURCE}_TO_${TARGET} APPEND PROPERTY LABELS "converter")
endfunction()
