#############################################################
# Add I/O tests to ctest
#############################################################

message("Adding I/O tests for QMCPACK")

function(
  RUN_RESTART_AND_CHECK
  base_name
  base_dir
  series
  input_file
  procs
  threads)

  # "run restart and check" function does 3 things:
  #  1. run qmcpack executable on $input_file.in.xml
  #  1. run qmcpack executable on $input_file.restart.xml in the same folder
  #  2. run $check_cmd located under tests/scripts folder

  # build test name
  set(full_name "${base_name}-r${procs}-t${threads}")
  message(VERBOSE "Adding test ${full_name}")

  # add run (task 1)
  set(test_added false)
  run_qmc_app(
    ${full_name}
    ${base_dir}
    ${procs}
    ${threads}
    test_added
    ${input_file}.in.xml)
  if(test_added)

    # add restart (task 2)
    set(restart_name "${base_name}-r${procs}-t${threads}-restart")
    #set (test_added false)
    run_qmc_app_no_copy(
      ${restart_name}
      ${CMAKE_CURRENT_BINARY_DIR}/${full_name}
      ${procs}
      ${threads}
      test_added
      ${input_file}.restart.xml)

    # make restart depend on the initial run
    set_property(
      TEST ${restart_name}
      APPEND
      PROPERTY DEPENDS ${full_name})

    math(EXPR series_plus_1 "${series} + 1")
    set(XXX "00${series}")
    set(XXY "00${series_plus_1}")
    string(LENGTH "${XXX}" len_XXX)
    math(EXPR start_XXX "${len_XXX} - 3")
    string(SUBSTRING "${XXX}" ${start_XXX} 3 XXX)
    string(LENGTH "${XXY}" len_XXY)
    math(EXPR start_XXY "${len_XXY} - 3")
    string(SUBSTRING "${XXY}" ${start_XXY} 3 XXY)
    set(check_cmd ${CMAKE_CURRENT_SOURCE_DIR}/../scripts/check_restart.sh ${input_file}.s${XXX} ${input_file}.s${XXY})
    #message(${check_cmd})

    # add test (task 3)
    set(test_name "${full_name}-check") # hard-code for single test
    set(work_dir "${CMAKE_CURRENT_BINARY_DIR}/${full_name}")
    #message(${work_dir})
    add_test(
      NAME "${test_name}"
      COMMAND ${check_cmd}
      WORKING_DIRECTORY "${work_dir}")

    # make test depend on the restart run
    set_property(TEST ${test_name} APPEND PROPERTY DEPENDS ${full_name}-restart)

  endif()

endfunction()

run_restart_and_check(
  deterministic-restart
  "${qmcpack_SOURCE_DIR}/tests/io/restart"
  1
  qmc
  2
  2)

run_restart_and_check(
  deterministic-restart
  "${qmcpack_SOURCE_DIR}/tests/io/restart"
  1
  qmc
  1
  4)

run_restart_and_check(
  deterministic-restart_batch
  "${qmcpack_SOURCE_DIR}/tests/io/restart_batch"
  1
  qmc_batch
  2
  2)

add_test_check_file_existence(deterministic-restart_batch-r2-t2 qmc_batch.s000.config.h5 TRUE)
add_test_check_file_existence(deterministic-restart_batch-r2-t2 qmc_batch.s000.random.h5 TRUE)
add_test_check_file_existence(deterministic-restart_batch-r2-t2 qmc_batch.s001.config.h5 FALSE)
add_test_check_file_existence(deterministic-restart_batch-r2-t2 qmc_batch.s001.random.h5 FALSE)

run_restart_and_check(
  deterministic-restart_batch
  "${qmcpack_SOURCE_DIR}/tests/io/restart_batch"
  1
  qmc_batch
  1
  4)

if(NOT QMC_MIXED_PRECISION)  
  run_restart_and_check(
    deterministic-restart_batch_opt
    "${qmcpack_SOURCE_DIR}/tests/io/restart_batch_opt"
    1
    qmc_batch
    2
    2)
  
  add_test_check_file_existence(deterministic-restart_batch_opt-r2-t2 qmc_batch.s000.config.h5 TRUE)
  add_test_check_file_existence(deterministic-restart_batch_opt-r2-t2 qmc_batch.s000.random.h5 TRUE)
  add_test_check_file_existence(deterministic-restart_batch_opt-r2-t2 qmc_batch.s001.config.h5 FALSE)
  add_test_check_file_existence(deterministic-restart_batch_opt-r2-t2 qmc_batch.s001.random.h5 FALSE)
  
  run_restart_and_check(
    deterministic-restart_batch_opt
    "${qmcpack_SOURCE_DIR}/tests/io/restart_batch_opt"
    1
    qmc_batch
    1
    4)
  
  add_test_check_file_existence(deterministic-restart_batch_opt-r1-t4 qmc_batch.s000.config.h5 TRUE)
  add_test_check_file_existence(deterministic-restart_batch_opt-r1-t4 qmc_batch.s000.random.h5 TRUE)
  add_test_check_file_existence(deterministic-restart_batch_opt-r1-t4 qmc_batch.s001.config.h5 FALSE)
  add_test_check_file_existence(deterministic-restart_batch_opt-r1-t4 qmc_batch.s001.random.h5 FALSE)
endif()

run_restart_and_check(
  deterministic-restart_dmc
  "${qmcpack_SOURCE_DIR}/tests/io/restart_dmc"
  2
  det_qmc_vmcbatch_dmcbatch_tm
  2
  2)

run_restart_and_check(
  deterministic-restart_dmc
  "${qmcpack_SOURCE_DIR}/tests/io/restart_dmc"
  2
  det_qmc_vmcbatch_dmcbatch_tm
  1
  4)

run_restart_and_check(
  deterministic-restart_dmc_disable_branching
  "${qmcpack_SOURCE_DIR}/tests/io/restart_dmc_disable_branching"
  2
  det_qmc_vmcbatch_dmcbatch_tm
  2
  2)

run_restart_and_check(
  deterministic-restart_dmc_disable_branching
  "${qmcpack_SOURCE_DIR}/tests/io/restart_dmc_disable_branching"
  2
  det_qmc_vmcbatch_dmcbatch_tm
  1
  4)

run_restart_and_check(
  deterministic-save_spline_coefs-batched
  "${qmcpack_SOURCE_DIR}/tests/io/save_spline_coefs"
  1
  qmc_batched
  2
  2)

run_restart_and_check(
  deterministic-save_spline_coefs-batched
  "${qmcpack_SOURCE_DIR}/tests/io/save_spline_coefs"
  1
  qmc_batched
  1
  4)

run_restart_and_check(
  deterministic-save_spline_coefs_hybridrep-batched
  "${qmcpack_SOURCE_DIR}/tests/io/save_spline_coefs_hybridrep"
  1
  qmc_batched
  2
  2)

run_restart_and_check(
  deterministic-save_spline_coefs_hybridrep-batched
  "${qmcpack_SOURCE_DIR}/tests/io/save_spline_coefs_hybridrep"
  1
  qmc_batched
  1
  3)

get_property(CURRENT_TESTS DIRECTORY PROPERTY TESTS)
add_labels_for_tests("${CURRENT_TESTS}")
