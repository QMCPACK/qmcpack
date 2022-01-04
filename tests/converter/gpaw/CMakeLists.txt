# add tests for GPAW4QMCPACK converter
file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/gpaw4qmcpack_test.py ${CMAKE_CURRENT_BINARY_DIR}/gpaw4qmcpack_test.py SYMBOLIC)
function(ADD_GPAW4QMCPACK_TEST test_name)
  set(EXE_NAME "${PROJECT_SOURCE_DIR}/src/QMCTools/gpaw4qmcpack.py")
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${test_name}")
  file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/restart.gpw"
	  "${CMAKE_CURRENT_BINARY_DIR}/${test_name}/restart.gpw" SYMBOLIC)
  file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/${test_name}/gold.orbs.h5"
	  "${CMAKE_CURRENT_BINARY_DIR}/${test_name}/gold.orbs.h5" SYMBOLIC)
  add_test(NAME gpaw4qmcpack_test COMMAND "${CMAKE_CURRENT_BINARY_DIR}/gpaw4qmcpack_test.py"
      "${CMAKE_CURRENT_BINARY_DIR}/${test_name}"
       --exe "${EXE_NAME}" --h5diff "${HDF5_DIFF_EXECUTABLE}")
  set_tests_properties(gpaw4qmcpack_test PROPERTIES
                         FAIL_REGULAR_EXPRESSION "  FAIL"
                         PASS_REGULAR_EXPRESSION "  pass"
                         TIMEOUT 120 LABELS "converter;deterministic")
endfunction()
add_gpaw4qmcpack_test(test_Si_diamond)
