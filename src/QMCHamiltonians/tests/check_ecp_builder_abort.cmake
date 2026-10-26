# MPI_Abort and serial aborts have platform-dependent exit statuses. Require the
# pool's fatal diagnostic as well as failure, so an unrelated crash cannot pass.
execute_process(COMMAND ${TEST_COMMAND} "[.ecp_abort]" --turn-on-printout
                RESULT_VARIABLE result OUTPUT_VARIABLE output ERROR_VARIABLE error TIMEOUT 30)
set(log "${output}\n${error}")
if("${result}" STREQUAL "0" OR "${result}" MATCHES "[Tt]imeout")
  message(FATAL_ERROR "Expected ECP application abort; got ${result}\n${log}")
endif()
foreach(expected "Ion species C charge 3 pseudopotential charge 4 mismatch!"
                 "Failed to add pseudopotential for element C"
                 "Fatal Error. Aborting at ECPotentialBuilder::useXmlFormat failed")
  string(FIND "${log}" "${expected}" found)
  if(found EQUAL -1)
    message(FATAL_ERROR "Missing abort diagnostic: ${expected}\nResult: ${result}\n${log}")
  endif()
endforeach()
