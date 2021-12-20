# add tests for GPAW4QMCPACK converter
function(ADD_GPAW4QMCPACK_TEST test_name)
      set(EXE_NAME ${qmcpack_BINARY_DIR}/bin/gpaw4qmcpack)
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/gpaw4qmcpack_test.py"
                          "${CMAKE_CURRENT_BINARY_DIR}")
  execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}"
                          "${CMAKE_CURRENT_BINARY_DIR}/${test_name}")
  add_test(NAME gpaw4qmcpack_test COMMAND "${CMAKE_CURRENT_BINARY_DIR}/gpaw4qmcpack_test.py"
      "${CMAKE_CURRENT_BINARY_DIR}/${test_name}"
       --exe "${EXE_NAME}" --h5diff "${HDF5_DIFF_EXECUTABLE}")
  set_tests_properties(gpaw4qmcpack_test PROPERTIES
                         FAIL_REGULAR_EXPRESSION "  FAIL"
                         PASS_REGULAR_EXPRESSION "  pass"
                         TIMEOUT 120 LABELS "converter;deterministic")
endfunction()
add_gpaw4qmcpack_test(test_Si_diamond)
