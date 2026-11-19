set(THISTESTNAME ppconvert_runs)

execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/O_test.sh ${CMAKE_CURRENT_BINARY_DIR})
execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/O.BFD.gamess ${CMAKE_CURRENT_BINARY_DIR})
#set(ARGS "\-\-gamess_pot O.BFD.gamess") #--s_ref;'1s(2)2p(4)';--p_ref;'1s(2)2p(4)';--d_ref;'1s(2)2p(4)';--xml O.BFD.xml)
add_test(
  NAME ${THISTESTNAME}
  COMMAND bash ${CMAKE_CURRENT_BINARY_DIR}/O_test.sh $<TARGET_FILE:ppconvert>
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(${THISTESTNAME} PROPERTIES LABELS "deterministic")

add_test(
  NAME ppconvert_o_diff
  COMMAND
    ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../numeric_diff.py --absolute-tolerance 0.001
    ${CMAKE_CURRENT_BINARY_DIR}/O.BFD.xml ${qmcpack_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp/O.BFD.xml
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(ppconvert_o_diff PROPERTIES DEPENDS ppconvert_runs LABELS "deterministic")

add_test(
  NAME ppconvert_invalid_option
  COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../expect_failure.py --contains
          "Unrecognized argument invalid-option" -- $<TARGET_FILE:ppconvert> --invalid-option)
set_tests_properties(ppconvert_invalid_option PROPERTIES LABELS "deterministic")

add_test(
  NAME ppconvert_invalid_channel
  COMMAND
    ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/../expect_failure.py --contains "Invalid local channel 99" --
    $<TARGET_FILE:ppconvert> --gamess_pot ${CMAKE_CURRENT_SOURCE_DIR}/O.BFD.gamess --local_channel 99)
set_tests_properties(ppconvert_invalid_channel PROPERTIES LABELS "deterministic")

add_test(
  NAME ppconvert_additional_outputs
  COMMAND
    $<TARGET_FILE:ppconvert> --gamess_pot ${CMAKE_CURRENT_SOURCE_DIR}/O.BFD.gamess --s_ref "1s(2)2p(4)" --p_ref
    "1s(2)2p(4)" --tm ${CMAKE_CURRENT_BINARY_DIR}/O.BFD.psp --upf ${CMAKE_CURRENT_BINARY_DIR}/O.BFD.upf)
set_tests_properties(ppconvert_additional_outputs PROPERTIES LABELS "deterministic")
