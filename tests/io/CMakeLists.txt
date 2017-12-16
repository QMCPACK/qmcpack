#############################################################
# Add I/O tests to ctest
#############################################################

message("Adding I/O tests for QMCPACK")

function(RUN_RESTART_AND_CHECK base_name base_dir input_file procs threads check_script use_for_coverage)

  # "run restart and check" function does 3 things:
  #  1. run qmcpack executable on $input_file.in.xml
  #  1. run qmcpack executable on $input_file.restart.xml in the same folder
  #  2. run $check_script located in the same folder

  # build test name
  set(full_name "${base_name}-${procs}-${threads}")
  message("Adding test ${full_name}")

  # add run (task 1)
  set (test_added false)
  RUN_QMC_APP(${full_name} ${base_dir} ${procs} ${threads} test_added ${input_file}.in.xml)
  if ( test_added )

    # add restart (task 2)
    set(restart_name "${base_name}-${procs}-${threads}-restart")
    #set (test_added false)
    RUN_QMC_APP_NO_COPY( ${restart_name} ${CMAKE_CURRENT_BINARY_DIR}/${full_name} ${procs} ${threads} test_added ${input_file}.restart.xml )

    # make restart depend on the initial run
    set_property(TEST ${restart_name} APPEND PROPERTY DEPENDS ${full_name})

    if ( use_for_coverage )
      set_property(TEST ${full_name} APPEND PROPERTY LABELS "coverage")
      set_property(TEST ${restart_name} APPEND PROPERTY LABELS "coverage")
    endif()

    # set up command to run check, assume check_script is in the same folder as input
    set(check_cmd ${CMAKE_CURRENT_BINARY_DIR}/${full_name}/${check_script})
    #message(${check_cmd})

    # add test (task 3)
    set(test_name "${full_name}-check") # hard-code for single test
    set(work_dir "${CMAKE_CURRENT_BINARY_DIR}/${full_name}")
    #message(${work_dir})
    add_test(NAME "${test_name}"
      COMMAND "${check_cmd}"
      WORKING_DIRECTORY "${work_dir}"
    )

    # make test depend on the restart run
    set_property(TEST ${test_name} APPEND PROPERTY DEPENDS ${full_name}-restart)

  endif()

endfunction()

RUN_RESTART_AND_CHECK( restart "${CMAKE_SOURCE_DIR}/tests/io/restart" qmc_short 8 2 check.sh false)
RUN_RESTART_AND_CHECK( restart "${CMAKE_SOURCE_DIR}/tests/io/restart" qmc_short 1 16 check.sh true)
