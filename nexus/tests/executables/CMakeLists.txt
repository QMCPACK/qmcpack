
# Tests for nexus/executables

set(QMCA_EXE "${qmcpack_SOURCE_DIR}/nexus/executables/qmca")
message("Adding qmca tests")

# ADD_QMCA_TEST(<test_name>
#               CMDARGS <arg> [<arg> ...]
#               PASS_REGULAR_EXPRESSION <regexp>
#               FAIL_REGULAR_EXPRESSION <regexp>
#              )
#
#  The test name is prefixed with 'qmca_'
#  CMDARGS are the arguments passed to qmca
#  PASS/FAIL_REGULAR_EXPRESSION get passed to test properties

function(ADD_QMCA_TEST test_name ${ARGN})

  set(oneValueArgs PASS_REGULAR_EXPRESSION FAIL_REGULAR_EXPRESSION)
  set(multiValueArgs CMDARGS)
  cmake_parse_arguments(QMCA "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  add_test(NAME qmca_${test_name} COMMAND ${QMCA_EXE} ${QMCA_CMDARGS})
  set_tests_properties(qmca_${test_name} PROPERTIES LABELS "qmca")

  if (QMCA_PASS_REGULAR_EXPRESSION)
    set_tests_properties(qmca_${test_name} PROPERTIES PASS_REGULAR_EXPRESSION ${QMCA_PASS_REGULAR_EXPRESSION}) 
  endif()

  if (QMCA_FAIL_REGULAR_EXPRESSION)
    set_tests_properties(qmca_${test_name} PROPERTIES FAIL_REGULAR_EXPRESSION ${QMCA_FAIL_REGULAR_EXPRESSION}) 
  endif()

endfunction()


execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}")

ADD_QMCA_TEST(bare PASS_REGULAR_EXPRESSION "Usage")

ADD_QMCA_TEST(file_does_not_exist CMDARGS not_here.s000.scalar.dat
              PASS_REGULAR_EXPRESSION "QMCA Error")

ADD_QMCA_TEST(empty_file CMDARGS empty.s000.scalar.dat
              PASS_REGULAR_EXPRESSION "Empty input file")

ADD_QMCA_TEST(label_only CMDARGS label_only.s000.scalar.dat
              PASS_REGULAR_EXPRESSION "Empty input file")

ADD_QMCA_TEST(one_line CMDARGS test1.s000.scalar.dat
              PASS_REGULAR_EXPRESSION "test1  series 0")

ADD_QMCA_TEST(two_files CMDARGS test1.s000.scalar.dat test1.s001.scalar.dat
              PASS_REGULAR_EXPRESSION "test1  series 0;test1  series 1")

ADD_QMCA_TEST(join_two CMDARGS -j 0,1 test1.s000.scalar.dat test1.s001.scalar.dat
              PASS_REGULAR_EXPRESSION "test1  series 0"
              FAIL_REGULAR_EXPRESSION "test1  series 1")
