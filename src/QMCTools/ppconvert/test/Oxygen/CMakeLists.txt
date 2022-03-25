set(THISTESTNAME ppconvert_runs)

execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/O_test.sh ${CMAKE_CURRENT_BINARY_DIR})
execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/O.BFD.gamess ${CMAKE_CURRENT_BINARY_DIR})
#set(ARGS "\-\-gamess_pot O.BFD.gamess") #--s_ref;'1s(2)2p(4)';--p_ref;'1s(2)2p(4)';--d_ref;'1s(2)2p(4)';--xml O.BFD.xml)
add_test(
  NAME ${THISTESTNAME}
  COMMAND bash ${CMAKE_CURRENT_BINARY_DIR}/O_test.sh $<TARGET_FILE:ppconvert>
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(${THISTESTNAME} PROPERTIES LABELS "deterministic")

find_program(HAVE_NUMDIFF numdiff)
if(HAVE_NUMDIFF)
  add_test(
    NAME ppconvert_o_diff
    COMMAND numdiff -a 0.001 ${CMAKE_CURRENT_BINARY_DIR}/O.BFD.xml
            ${qmcpack_SOURCE_DIR}/tests/solids/monoO_1x1x1_pp/O.BFD.xml
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_tests_properties(ppconvert_o_diff PROPERTIES DEPENDS ppconvert_runs LABELS "deterministic")
endif(HAVE_NUMDIFF)
